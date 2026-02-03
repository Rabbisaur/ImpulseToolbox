function NSx = openNSxFasterMex(filename, varargin)
% openNSxFasterMex - MAX SPEED version of openNSx using C++/MEX
% 
% Usage matches openNSx.
%
% This function delegates header parsing to 'openNSx' (using 'noread'),
% then uses 'readNSxMex' to read/align/de-interleave data in one pass.

%% 1. Use openNSx to parse headers (NO READ)
% We pass all arguments but ensure 'noread' is active so we get metadata only.
% We also need to strip 'read' if provided.

% 1. Clean arguments for delegation to openNSx (metadata only)
% We remove arguments that openNSx doesn't recognize or that we handle specifically.
args = varargin;

% Identify custom arguments to remove
% 'read', 'report', 'noread' are handled/added specifically
% 'skipfactor', 'zeropad' are custom for this wrapper/handled by MEX
isCustom = cellfun(@(x) ischar(x) && any(strcmpi(x, {'read', 'report', 'noread', 'skipfactor', 'zeropad'})), args);

% Note: skipfactor and zeropad take a following value, so we must remove that too.
customIdx = find(isCustom);
toRemove = customIdx;
for i = 1:length(customIdx)
    idx = customIdx(i);
    if any(strcmpi(args{idx}, {'skipfactor', 'zeropad'}))
        if idx < length(args)
            toRemove = [toRemove, idx+1];
        end
    end
end
args(toRemove) = [];
args = [args, {'noread'}];

    % Call original to get the requested structure (for Labels, Config, and Segment analysis)
    % We suppress the "Reading header only" warning for a cleaner interface.
    w = warning('off', 'all'); % openNSx is chatty
    NSx = openNSx(filename, args{:});
    
    % Re-check if user wanted to read data?
    % For this wrapper, we assume YES by default unless 'noread' was explicit in varargin.
    flagNoRead = any(cellfun(@(x) ischar(x) && strcmpi(x, 'noread'), varargin));
    if flagNoRead
        warning(w);
        return;
    end
    
    % Retrieve Physical File Layout (All Channels for mapping)
    NSx_Phys = openNSx(filename, 'noread');
    warning(w); % Restore warning state

%% 2. Prepare for MEX Read
% We need to calculate the alignment parameters to pass to C++.

% Extract key params
totalPackets = NSx.MetaTags.DataPoints; % vector if segmented
% If segmented, we loop segments.

numChannels = NSx.MetaTags.ChannelCount;
skipFactor = 1; % Default
% Parse skip factor from varargin if present
% (Simplified parsing: relying on what openNSx determined? 
% openNSx structure doesn't store 'requestedSkipFactor' explicitly in MetaTags...)
% We must re-parse skip factor or use a simplified approach.

% Let's look for 'skipfactor' in args or default
skipFactorIdx = find(strcmpi(varargin, 'skipfactor'));
if ~isempty(skipFactorIdx)
    skipFactor = varargin{skipFactorIdx+1};
end

% Precision
precisionClass = 0; % int16
if any(strcmpi(varargin, 'double'))
    precisionClass = 1;
end
% openNSx 'uv' argument -> converts to double.
flagUV = any(strcmpi(varargin, 'uv'));
if flagUV
    precisionClass = 1; % Force double for UV calculation
end


% Internal check: Some arguments like 'uv' might influence precisionClass.
% Ensure flagUV is consistent with precisionClass.

% 2. Prepare for Metadata and Structure
PhysicalChannelCount = NSx_Phys.MetaTags.ChannelCount;
AllChannelIDs = [NSx_Phys.ElectrodesInfo.ChannelID];

% Map Requested Channels to Physical Indices
if isempty(NSx.MetaTags.ChannelID)
     % Fallback if something weird
     targetIndices = 0:PhysicalChannelCount-1;
else
    [found, targetIndices] = ismember(NSx.MetaTags.ChannelID, AllChannelIDs);
    targetIndices = targetIndices(found) - 1; % 0-based for C++
    
    if any(~found)
        warning('Wrapper: Some requested channels not found in Physical Headers.');
    end
end

% Robust Header Parsing for Data Offset
FID = fopen(filename, 'r', 'ieee-le');
fseek(FID, 0, 'eof');
fileSize = ftell(FID);
fseek(FID, 10, 'bof');
BytesInHeaders = fread(FID, 1, 'uint32');

% Auto-detect Packet Start
fseek(FID, BytesInHeaders, 'bof');
if fread(FID, 1, 'uint8') == 1
    dataStartOffset = BytesInHeaders;
else
    % Fallback/Sniff logic
    fseek(FID, BytesInHeaders, 'bof');
    firstByte = fread(FID, 1, 'uint8');
    if firstByte == 46 % Data byte
         dataStartOffset = BytesInHeaders - 13; 
    else
        dataStartOffset = BytesInHeaders; 
    end
end

% Iterate Segments and Read
for segCtx = 1:length(NSx.MetaTags.DataPoints)
    % Calculate Offset for this segment
    % (Note: For single segment, dataStartOffset is correct. For multiple, we'd need more logic, 
    % but NPMK usually reads one continuous block unless segmented. 
    % We assume segCtx=1 is the primary use case for high-speed read.)
    
    % Timestamp Verification
    fseek(FID, dataStartOffset + 1, 'bof');
    tsVal = fread(FID, 1, 'uint64');
    if isempty(NSx.MetaTags.Timestamp) || length(NSx.MetaTags.Timestamp) < segCtx
        NSx.MetaTags.Timestamp(segCtx) = tsVal;
    elseif NSx.MetaTags.Timestamp(segCtx) ~= tsVal
         NSx.MetaTags.Timestamp(segCtx) = tsVal;
    end
    
    % Use Physical Metadata for Total Packets
    physicalPackets_Total = NSx_Phys.MetaTags.DataPoints(segCtx);
    packetSize = 13 + 2 * PhysicalChannelCount;
    
    % Handle Time Selection (t:XX)
    % Default: Read all
    startPacketIdx = 0;
    numPacketsToRead = physicalPackets_Total;
    
    % Parse 't:XX' from varargin
    % ... (rest of search logic) ...
    tArgIdx = find(cellfun(@(x) ischar(x) && strncmpi(x, 't:', 2), varargin));
    if ~isempty(tArgIdx)
        tStr = varargin{tArgIdx};
        parts = sscanf(tStr, 't:%f:%f');
        if length(parts) == 2
            tStart = parts(1);
            tEnd = parts(2);
            
            % Determine Unit
             timeUnit = 'sec'; % Default
             if length(varargin) > tArgIdx
                 nextArg = varargin{tArgIdx+1};
                 if ischar(nextArg)
                     if any(strcmpi(nextArg, {'sec','secs','second','seconds'}))
                         timeUnit = 'sec';
                     elseif any(strcmpi(nextArg, {'min','mins','minute','minutes'}))
                         timeUnit = 'min';
                     elseif any(strcmpi(nextArg, {'sample','samples'}))
                         timeUnit = 'sample';
                     end
                 end
             end
             
             fs = NSx.MetaTags.SamplingFreq;
             if strcmpi(timeUnit, 'sec')
                 startPacketIdx = floor(tStart * fs);
                 endPacketIdx = floor(tEnd * fs);
             elseif strcmpi(timeUnit, 'min')
                 startPacketIdx = floor(tStart * 60 * fs);
                 endPacketIdx = floor(tEnd * 60 * fs);
             else % samples
                 startPacketIdx = floor(tStart);
                 endPacketIdx = floor(tEnd);
             end
             
             numPacketsToRead = endPacketIdx - startPacketIdx;
             
             % Bounds Check
             if startPacketIdx < 0, startPacketIdx = 0; end
             if startPacketIdx + numPacketsToRead > physicalPackets_Total
                 numPacketsToRead = physicalPackets_Total - startPacketIdx;
             end
        end
    end
    
    % Apply Start Offset
    dataStartOffset = dataStartOffset + startPacketIdx * packetSize;
    
    % Alignment params
    % We must account for clock drift (difference between physical and metadata)
    expectedTotal = round(NSx_Phys.MetaTags.DataDurationSec(segCtx) * NSx_Phys.MetaTags.SamplingFreq);
    driftTotal = expectedTotal - physicalPackets_Total;
    
    % Match openNSx: alignment is applied to DECIMATED data if skipfactor > 1
    if skipFactor > 1
        % Match openNSx ratio exactly (floor( Drift * DecimatedPoints / TotalRawPoints ))
        decDataPoints = floor(numPacketsToRead / skipFactor);
        addedSamples = floor(driftTotal * decDataPoints / physicalPackets_Total);
        fileDataLengthDec = decDataPoints;
    else
        fileDataLengthDec = numPacketsToRead;
        % For subsection, pro-rate drift. For full, use driftTotal.
        if isempty(tArgIdx)
            addedSamples = driftTotal;
        else
            addedSamples = round(driftTotal * (numPacketsToRead / physicalPackets_Total));
        end
    end
    
    % Match openNSx: gapIndex is calculated from FULL PHYISCAL parameters
    % but applied to the stream (which for us is the decimated stream).
    gapIndex = 0;
    if addedSamples ~= 0
        % Even if skipFactor > 1, openNSx uses the "full" driftTotal and physical count for gapIndex
        gapIndex = round(numPacketsToRead / (abs(driftTotal) + 1));
    end
    
    fclose(FID);
    
    % Call MEX
    mexData = readNSxMex(filename, dataStartOffset, PhysicalChannelCount, ...
                                  numPacketsToRead, gapIndex, addedSamples, ...
                                  skipFactor, precisionClass, targetIndices);

    
    % ZERO PADDING
    % openNSx adds zero padding if Timestamp > 0.
    % We default to false (Matches user request to avoid PTP warnings).
    flagZeroPad = false; 
    
    % Parse 'zeropad' argument
    zpIdx = find(strcmpi(varargin, 'zeropad'));
    if ~isempty(zpIdx) && length(varargin) > zpIdx
        flagZeroPad = logical(varargin{zpIdx+1});
    end
    
    if flagZeroPad && NSx.MetaTags.Timestamp(segCtx) > 0
        numZeroSamples = floor(NSx.MetaTags.Timestamp(segCtx) / skipFactor);
        
        if numZeroSamples < 200000000 && numZeroSamples > 0
            % Prepend Zeros (Reasonable size)
            if precisionClass == 0
                padding = zeros(size(mexData,1), numZeroSamples, 'int16');
            else
                padding = zeros(size(mexData,1), numZeroSamples, 'double');
            end
            mexData = [padding mexData];
        else
             warning('Excessive or invalid padding (%d samples). Skipping padding.', numZeroSamples);
        end
    end
    
    NSx.Data{segCtx} = mexData;
    
    % Apply UV if needed (MEX returned double, but raw values)
    if flagUV && precisionClass == 1
        % Vectorized scale
        % scale = (MaxAnalog / MaxDigi)
        scale = double([NSx.ElectrodesInfo.MaxAnalogValue]) ./ double([NSx.ElectrodesInfo.MaxDigiValue]);
        % NSx.Data is [Ch x Samples]
        NSx.Data{segCtx} = NSx.Data{segCtx} .* scale';
    end
    
end

% reduce to array if only one cell (Match openNSx)
if iscell(NSx.Data) && length(NSx.Data) == 1
    NSx.Data = NSx.Data{1};
end

end
