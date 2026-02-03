/*
 * readNSxMex.cpp
 * 
 * Reads NSx data directly from file, de-interleaves channels, and applies 
 * clock drift alignment (gap insertion/deletion) in a single pass.
 *
 * Usage:
 * Data = readNSxMex(filename, dataStartOffset, numChannels, totalPackets, ...
 *                   gapIndex, addedSamples, skipFactor, precisionClass)
 *
 * precisionClass: 0 = int16, 1 = double
 */

#include "mex.h"
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cmath>

// 5MB Chunk size for reading
const size_t CHUNK_SIZE_BYTES = 1024 * 1024 * 5; 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check inputs
    if (nrhs < 8) {
        mexErrMsgIdAndTxt("readNSxMex:invalidNumInputs", "8 inputs required (9 optional).");
    }

    // 1. Filename
    char* filename = mxArrayToString(prhs[0]);
    if (filename == NULL) mexErrMsgIdAndTxt("readNSxMex:invalidInput", "Filename must be string.");

    // 2. Data Start Offset (bytes)
    double dataStartOffset = mxGetScalar(prhs[1]);

    // 3. Physical Channels (Total in File)
    int physicalChannels = (int)mxGetScalar(prhs[2]);

    // 4. Total Packets to Read (Input file packets)
    size_t totalPacketsInput = (size_t)mxGetScalar(prhs[3]);

    // 5. Gap Index (Alignment)
    size_t gapIndex = (size_t)mxGetScalar(prhs[4]);

    // 6. Added Samples (Alignment)
    int addedSamples = (int)mxGetScalar(prhs[5]);

    // 7. Skip Factor
    int skipFactor = (int)mxGetScalar(prhs[6]);

    // 8. Output Precision (0=int16, 1=double)
    int precisionClass = (int)mxGetScalar(prhs[7]);

    // 9. Target Indices (Optional)
    std::vector<int> targetIndices;
    if (nrhs > 8 && !mxIsEmpty(prhs[8])) {
        if (!mxIsDouble(prhs[8])) mexErrMsgIdAndTxt("readNSxMex:invalidType", "targetIndices must be double array.");
        double* idxPtr = mxGetPr(prhs[8]);
        size_t nT = mxGetNumberOfElements(prhs[8]);
        targetIndices.reserve(nT);
        for (size_t k=0; k<nT; k++) {
            targetIndices.push_back((int)idxPtr[k]);
        }
    } else {
        targetIndices.reserve(physicalChannels);
        for (int k=0; k<physicalChannels; k++) targetIndices.push_back(k);
    }
    int outputChannels = targetIndices.size();

    // Packet Size Calculation (Assuming PTP / One Sample Per Packet)
    // Packet: Header(1) + Timestamp(8) + Samples(4) + Data(2*physicalChannels)
    // Total: 13 + 2*physicalChannels
    size_t packetHeaderSize = 13; 
    size_t packetSize = packetHeaderSize + physicalChannels * 2;

    // Calculate Output Size (Final columns)
    size_t decimatedInputPackets = totalPacketsInput / skipFactor;
    size_t outputCols = decimatedInputPackets + addedSamples;
    
    if (outputCols > 2147483647) {
        mexPrintf("MEX WARNING: Output cols > 2^31. Ensure -largeArrayDims used.\n");
    }

    // Prepare Output Array
    if (precisionClass == 0) {
        plhs[0] = mxCreateNumericMatrix(outputChannels, outputCols, mxINT16_CLASS, mxREAL);
    } else {
        plhs[0] = mxCreateDoubleMatrix(outputChannels, outputCols, mxREAL);
    }

    void* outDataPtr = mxGetData(plhs[0]);
    int16_t* outInt16 = (int16_t*)outDataPtr;
    double* outDouble = (double*)outDataPtr;

    // Open File
    FILE* fid = fopen(filename, "rb");
    if (!fid) {
        mxFree(filename);
        mexErrMsgIdAndTxt("readNSxMex:fileError", "Could not open file.");
    }
    
#ifdef _WIN32
    _fseeki64(fid, (long long)dataStartOffset, SEEK_SET);
#else
    fseeko(fid, (off_t)dataStartOffset, SEEK_SET);
#endif

    // Buffer for Raw Read
    size_t packetsPerChunk = CHUNK_SIZE_BYTES / packetSize;
    if (packetsPerChunk == 0) packetsPerChunk = 1;
    size_t bytesPerChunk = packetsPerChunk * packetSize;
    std::vector<uint8_t> buffer(bytesPerChunk);

    size_t packetsReadSoFar = 0;
    
    // Alignment Internal State
    size_t samplesUntilEvent = gapIndex; 
    int eventCount = 0;
    int maxEvents = std::abs(addedSamples);
    bool isAdding = (addedSamples > 0);
    bool isRemoving = (addedSamples < 0);
    
    size_t writtenCols = 0;

    // Loop chunks (Single Pass)
    while (packetsReadSoFar < totalPacketsInput) {
        size_t packetsToRead = std::min(packetsPerChunk, totalPacketsInput - packetsReadSoFar);
        size_t bytesToRead = packetsToRead * packetSize;
        
        size_t r = fread(buffer.data(), 1, bytesToRead, fid);
        if (r < bytesToRead) {
            packetsToRead = r / packetSize; 
        }

        for (size_t p = 0; p < packetsToRead; p++) {
            uint8_t* pktPtr = buffer.data() + p * packetSize;
            
            // Match openNSx: Alignment is applied to THE DECIMATED STREAM
            if (packetsReadSoFar % skipFactor == 0) {
                // Determine actions (Skip/Normal/Duplicate)
                bool skipThisInput = false;
                bool duplicateThisInput = false;

                if (gapIndex > 0) {
                    if (isRemoving && eventCount < maxEvents) {
                        if (samplesUntilEvent == 1) {
                             skipThisInput = true;
                             samplesUntilEvent = gapIndex;
                             eventCount++;
                        } else {
                             samplesUntilEvent--;
                        }
                    } else if (isAdding && eventCount < maxEvents) {
                        if (samplesUntilEvent == 1) {
                             duplicateThisInput = true;
                             samplesUntilEvent = gapIndex;
                             eventCount++;
                        } else {
                             samplesUntilEvent--;
                        }
                    }
                }

                int packetsToProduce = 1;
                if (skipThisInput) packetsToProduce = 0;
                if (duplicateThisInput) packetsToProduce = 2;

                for (int k = 0; k < packetsToProduce; k++) {
                    if (writtenCols < outputCols) {
                        // WRITE DATA (De-interleave)
                        int16_t* chansPtr = (int16_t*)(pktPtr + 13);
                        size_t baseIdx = writtenCols * outputChannels;
                        
                        if (precisionClass == 0) {
                            for (int i=0; i<outputChannels; i++) {
                                outInt16[baseIdx + i] = chansPtr[targetIndices[i]];
                            }
                        } else {
                            for (int i=0; i<outputChannels; i++) {
                                outDouble[baseIdx + i] = (double)chansPtr[targetIndices[i]];
                            }
                        }
                        writtenCols++;
                    }
                }
            }
            packetsReadSoFar++;
        }
        if (ftell(fid) >= 0 && r == 0) break; 
    }

    mxFree(filename);
    fclose(fid);
}
