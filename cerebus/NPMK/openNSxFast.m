function NSx = openNSxFast(filename, varargin)
% openNSxFast - Optimized pure-MATLAB NSx reader
% 
% Usage: NSx = openNSxFast(filename, 'uv')
%
% This is a simplified, optimized reader that's faster than openNSx
% but doesn't require MEX files. Supports basic functionality.

    %% Parse arguments
    flagUV = any(strcmpi(varargin, 'uv'));
    
    %% Open file
    fid = fopen(filename, 'r', 'ieee-le');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    fprintf('[openNSxFast] Reading headers from %s\n', filename);
    
    try
        %% Read File Type ID
        fseek(fid, 0, 'bof');
        FileTypeID = fread(fid, [1,8], '*char');
        
        %% Determine file format and parse accordingly
        if strcmpi(FileTypeID, 'NEURALSG')
            % File Spec 2.1 format
            basicHeaderBytes = fread(fid, 24, '*uint8');
            Label = char(basicHeaderBytes(1:16))';
            Label = Label(Label ~= 0);
            TimeRes = 30000;
            Period = double(typecast(basicHeaderBytes(17:20), 'uint32'));
            ChannelCount = double(typecast(basicHeaderBytes(21:24), 'uint32'));
            
            % Calculate sampling frequency
            if Period > 0
                Fs = TimeRes / Period;
            else
                Fs = 30000;
            end
            
            % Read Extended Header (simplified for fast reading)
            extendedHeaderSize = ChannelCount * 4;
            extendedHeaderBytes = fread(fid, extendedHeaderSize, '*uint8');
            ChannelIDs = double(typecast(extendedHeaderBytes, 'uint32'));  % Column vector
            
            % For NEURALSG, we don't have detailed electrode info, so create minimal structure
            ElectrodesInfo = struct([]);
            for i = 1:ChannelCount
                ElectrodesInfo(i).ChannelID = ChannelIDs(i);
                ElectrodesInfo(i).MaxDigiValue = 32767;  % Default int16 max
                ElectrodesInfo(i).MaxAnalogValue = 8191;  % Default analog max
            end
            MaxAnalogVals = [ElectrodesInfo.MaxAnalogValue];
            MaxDigiVals = [ElectrodesInfo.MaxDigiValue];
            
            BytesInHeaders = 8 + 24 + extendedHeaderSize;
            
        elseif strcmpi(FileTypeID, 'NEURALCD') || strcmpi(FileTypeID, 'BRSMPGRP')
            % File Spec 2.2, 2.3, and 3.0 format
            basicHeaderBytes = fread(fid, 306, '*uint8');
            Label = char(basicHeaderBytes(7:22))';
            Label = Label(Label ~= 0);
            
            % Read TimeRes and Period from correct byte positions
            TimeRes = double(typecast(basicHeaderBytes(283:286), 'uint32'));
            Period = double(typecast(basicHeaderBytes(279:282), 'uint32'));
            ChannelCount = double(typecast(basicHeaderBytes(303:306), 'uint32'));
            
            % Calculate sampling frequency
            if Period > 0
                Fs = 30000 / Period;  % Different formula for NEURALCD!
            else
                Fs = 30000;
            end
            
            % Read Extended Header
            extHeaderEntrySize = 66;
            extendedHeaderSize = ChannelCount * extHeaderEntrySize;
            extendedHeaderBytes = fread(fid, extendedHeaderSize, '*uint8');
            
            ElectrodesInfo = struct([]);
            ChannelIDs = zeros(ChannelCount, 1);  % Column vector
            MaxAnalogVals = zeros(1, ChannelCount);
            MaxDigiVals = zeros(1, ChannelCount);
            
            for i = 1:ChannelCount
                byteOffset = (i-1) * extHeaderEntrySize;
                Type = char(extendedHeaderBytes((1:2)+byteOffset))';
                
                if strcmp(Type, 'CC')
                    ElectrodesInfo(i).Type = 'CC';
                    ElectrodesInfo(i).ChannelID = double(typecast(extendedHeaderBytes((3:4)+byteOffset), 'uint16'));
                    ChannelIDs(i) = ElectrodesInfo(i).ChannelID;
                    
                    ElectrodesInfo(i).Label = char(extendedHeaderBytes((5:20)+byteOffset))';
                    ElectrodesInfo(i).Label = ElectrodesInfo(i).Label(ElectrodesInfo(i).Label ~= 0);
                    
                    ElectrodesInfo(i).MinDigiValue = double(typecast(extendedHeaderBytes((23:24)+byteOffset), 'int16'));
                    ElectrodesInfo(i).MaxDigiValue = double(typecast(extendedHeaderBytes((25:26)+byteOffset), 'int16'));
                    MaxDigiVals(i) = ElectrodesInfo(i).MaxDigiValue;
                    
                    ElectrodesInfo(i).MinAnalogValue = double(typecast(extendedHeaderBytes((27:28)+byteOffset), 'int16'));
                    ElectrodesInfo(i).MaxAnalogValue = double(typecast(extendedHeaderBytes((29:30)+byteOffset), 'int16'));
                    MaxAnalogVals(i) = ElectrodesInfo(i).MaxAnalogValue;
                    
                    ElectrodesInfo(i).Units = char(extendedHeaderBytes((31:46)+byteOffset))';
                    ElectrodesInfo(i).Units = ElectrodesInfo(i).Units(ElectrodesInfo(i).Units ~= 0);
                    
                    ElectrodesInfo(i).HighFreqCorner = double(typecast(extendedHeaderBytes((47:50)+byteOffset), 'uint32'));
                    ElectrodesInfo(i).HighFreqOrder = double(typecast(extendedHeaderBytes((51:54)+byteOffset), 'uint32'));
                    ElectrodesInfo(i).HighFilterType = double(typecast(extendedHeaderBytes((55:56)+byteOffset), 'uint16'));
                    ElectrodesInfo(i).LowFreqCorner = double(typecast(extendedHeaderBytes((57:60)+byteOffset), 'uint32'));
                    ElectrodesInfo(i).LowFreqOrder = double(typecast(extendedHeaderBytes((61:64)+byteOffset), 'uint32'));
                    ElectrodesInfo(i).LowFilterType = double(typecast(extendedHeaderBytes((65:66)+byteOffset), 'uint16'));
                end
            end
            
            BytesInHeaders = 8 + 306 + extendedHeaderSize;
        else
            error('Unsupported file format: %s', FileTypeID);
        end
        
        %% Determine data location
        fseek(fid, BytesInHeaders, 'bof');
        
        % Check for data packet header
        PacketID = fread(fid, 1, 'uint8=>double');
        if PacketID == 1
            % Standard packet format
            Timestamp = fread(fid, 1, 'uint32=>double');
            NumDataPoints = fread(fid, 1, 'uint32=>double');
        else
            % Legacy format - rewind
            fseek(fid, BytesInHeaders, 'bof');
            Timestamp = 0;
            % Determine data points from file size
            fseek(fid, 0, 'eof');
            fileSize = ftell(fid);
            dataBytes = fileSize - BytesInHeaders;
            NumDataPoints = dataBytes / (ChannelCount * 2); % 2 bytes per int16 sample
            fseek(fid, BytesInHeaders, 'bof');
        end
        
        fprintf('[openNSxFast] Reading %d samples from %d channels\n', NumDataPoints, ChannelCount);
        
        %% Read data efficiently
        % Position at data start
        if PacketID == 1
            dataStart = BytesInHeaders + 9; % Header + packet header
        else
            dataStart = BytesInHeaders;
        end
        fseek(fid, dataStart, 'bof');
        
        % Read all data at once - interleaved format
        fprintf('[openNSxFast] Reading data...\n');
        tic;
        rawData = fread(fid, [ChannelCount, NumDataPoints], 'int16=>int16');
        readTime = toc;
        fprintf('[openNSxFast] Data read completed in %.2f seconds\n', readTime);
        
        % Convert to double if UV conversion requested
        if flagUV
            fprintf('[openNSxFast] Converting to microvolts...\n');
            tic;
            % Calculate scale factors
            scale = MaxAnalogVals ./ MaxDigiVals;
            scale = scale(:); % Column vector
            
            % Apply scaling
            NSxData = double(rawData) .* scale;
            convTime = toc;
            fprintf('[openNSxFast] Conversion completed in %.2f seconds\n', convTime);
        else
            NSxData = rawData;
        end
        
        fclose(fid);
        
        %% Build output structure (compatible with openNSx format)
        NSx.MetaTags.FileTypeID = FileTypeID;
        NSx.MetaTags.SamplingLabel = Label;
        NSx.MetaTags.ChannelCount = ChannelCount;
        NSx.MetaTags.SamplingFreq = Fs;
        NSx.MetaTags.TimeRes = TimeRes;
        NSx.MetaTags.ChannelID = ChannelIDs;
        NSx.MetaTags.DateTime = '';  % Not parsed in fast version
        NSx.MetaTags.DataPoints = NumDataPoints;
        NSx.MetaTags.DataDurationSec = NumDataPoints / Fs;
        NSx.MetaTags.Timestamp = Timestamp;
        
        NSx.ElectrodesInfo = ElectrodesInfo;
        NSx.Data = NSxData;
        
        fprintf('[openNSxFast] Complete! Data size: [%d x %d]\n', size(NSxData, 1), size(NSxData, 2));
        
    catch ME
        fclose(fid);
        rethrow(ME);
    end
end