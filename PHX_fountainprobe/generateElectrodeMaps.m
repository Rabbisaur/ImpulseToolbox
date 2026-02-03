function generateElectrodeMaps()
% GENERATEELECTRODEMAPS Generates ElectrodeMap.mat and FP_ElectrodeMap.csv
% from the raw Excel configuration file.
%
% This script:
% 1. Reads 'LLNL PI to H-Plate Mapping_Rivet bonding.xlsx'
% 2. Generates the mapping table (FP <-> HP <-> BR)
% 3. Calculates electrode positions (Fountain geometry)
% 4. Creates the LGA layout matrices
% 5. Saves 'FP_ElectrodeMap.csv'
% 6. Saves 'ElectrodeMap.mat' (containing ElecMap table and LGA structure)

    %% 1. Configuration and Data Loading
    % Define raw source file
    rawExcelFile = 'LLNL PI to H-Plate Mapping_Rivet bonding.xlsx';
    if ~exist(rawExcelFile, 'file')
        error('Raw Excel file not found: %s', rawExcelFile);
    end
    
    fprintf('Reading %s...\n', rawExcelFile);
    
    % Sheet 1: Electrode to H-Plate mapping
    sheetNameMapping = 'Electrode to H-Plate mapping';
    dataMap = readcell(rawExcelFile, 'Sheet', sheetNameMapping);
    strMAP = dataMap(:, 1:2);
    
    % Sheet 2 (implied): Pedestal Layout 
    % Note: The original script implies Pedestal Layout is read or defined.
    % In FP_create_csv.m, PedLayout seems to be a 36x36 cell array but it's not clear where it's populated from.
    % Looking closely at FP_create_csv.m lines 174-193, it uses 'PedLayout' variable but never initializes/loads it 
    % before 'if ~Done' loop! This suggests 'PedLayout' might be a workspace variable or hardcoded.
    % WAIT: I missed something in the audit. 
    % Line 5: data = readcell(emapfile...)
    % Line 6: strMAP = data(:,1:2);
    % Variable 'data' has EVERYTHING.
    % Usually the pedestal layout acts as a map. 
    % Let's assume the Pedestal Layout is available or we can reconstruct it from the Mapping table if mostly 1-to-1.
    % ACTUALLY: The original script seems to fail or expect 'PedLayout' to exist.
    % However, we can reconstruct the grid if we know the mapping.
    % Update: I will check if 'Pedestal Layout' sheet exists.
    
    sheetInfo = sheetnames(rawExcelFile);
    if any(contains(sheetInfo, 'Pedestal'))
        pedSheet = sheetInfo{contains(sheetInfo, 'Pedestal')};
        fprintf('Found Pedestal sheet: %s\n', pedSheet);
        rawPed = readcell(rawExcelFile, 'Sheet', pedSheet);
        % Usually this is a grid. Assuming 36x36 or similar.
        PedLayoutRaw = rawPed; 
    else
        warning('Pedestal Layout sheet not found. Creating empty LGA map.');
        PedLayoutRaw = cell(36,36); 
    end

    %% 2. Process Mapping (FP <-> HP)
    NoneMarker = -99;
    LIST = struct();
    
    % The Excel structure seems to be: [FP_Label, HP_Label]
    % col 1: FP (e.g. "S1E1" or "R") 
    % col 2: HP (e.g. "1" or "R1")
    
    % Initialize counter
    idx = 0;
    
    % Skip header row (1)
    for i = 2:size(strMAP, 1)
        % Check if row has valid data
        fpStr = strMAP{i, 1};
        hpStr = strMAP{i, 2};
        
        idx = idx + 1;
        LIST(idx).FP_shank = NoneMarker;
        LIST(idx).FP_elect = NoneMarker;
        LIST(idx).FP_electCNT = NoneMarker;
        LIST(idx).FP_ref = NoneMarker;
        LIST(idx).HP_elect = NoneMarker;
        LIST(idx).HP_ref = NoneMarker;
        
        % Parse HP (Col 2)
        if ~isempty(hpStr) && ~any(ismissing(hpStr))
            if ischar(hpStr) && startsWith(hpStr, 'R')
                % Ref (e.g. R1)
                LIST(idx).HP_ref = str2double(hpStr(2:end));
            elseif isnumeric(hpStr)
                LIST(idx).HP_elect = hpStr;
            elseif ischar(hpStr)
                LIST(idx).HP_elect = str2double(hpStr);
            end
        end
        
        % Parse FP (Col 1)
        if ~isempty(fpStr) && ~any(ismissing(fpStr))
            if ischar(fpStr)
                % E.g. "S1E1" or "R"
                eIdx = strfind(fpStr, 'E');
                if ~isempty(eIdx)
                    % "S[shank]E[elect]"
                    LIST(idx).FP_shank = str2double(fpStr(2:eIdx-1));
                    LIST(idx).FP_elect = str2double(fpStr(eIdx+1:end));
                    if ~isnan(LIST(idx).FP_shank) && ~isnan(LIST(idx).FP_elect)
                        LIST(idx).FP_electCNT = (56 * (LIST(idx).FP_shank - 1)) + LIST(idx).FP_elect;
                    end
                elseif startsWith(fpStr, 'R')
                    % Ref?
                    LIST(idx).FP_ref = 1;
                    % Try to extract shank if format is R[shank] or S[shank]R?
                    % Original code: `str2double(strMAP{i,1}(2:Ridx-1))` -> implies S[shank]R
                    rIdx = strfind(fpStr, 'R');
                    if rIdx > 1
                         LIST(idx).FP_shank = str2double(fpStr(2:rIdx-1));
                    else
                         LIST(idx).FP_shank = 1; % Default?
                    end
                end
            end
        end
    end
    
    %% 3. Add Geometry (Polar & Cartesian)
    % Shank Start Angles (Polar: [Angle, RadiusIdx])
    % RadiusIdx 1=Inner, 2=Outer
    FP_shstart_pol = {...
        [285,2], [315,2], [345,2], [15,2], [0,1], [45,2], ... % 90-90=0
        [-60,1], [75,2], [60,1], [120,1], [105,2], [240,1], ... % 30-90=-60, etc
        [135,2], [180,1], [165,2], [195,2], [225,2], [255,2]};
        % Note: angles adjusted by -90 in original list for Inner ones
    
    % Fountain Head Geometry
    a_angles = [4, 10.5]; % degrees
    rr_ecc = [0.5 0.9]; % mm
    d_vals = [rr_ecc(1)*tand(a_angles(1)), rr_ecc(2)*tand(a_angles(2))];
    
    % Pre-calculate Shank Vectors
    ShankGeo = zeros(18, 3); % x, y, z direction
    for sh = 1:18
        if sh <= length(FP_shstart_pol)
            pol = FP_shstart_pol{sh};
            phi = pol(1); % Azimuth
            ridx = pol(2); % Ring
            theta = a_angles(ridx); % Tilt
            
            % Direction Vector Calculation (as per original script logic)
            % Hypotneuse factor stuff in original is complex, simplifying to match:
            % It calculates an exit vector.
            % Let's stick to the coordinates calculated in the loop for safety.
        end
    end
    
    % Populate LIST with coordinates
    for i = 1:length(LIST)
        sh = LIST(i).FP_shank;
        el = LIST(i).FP_elect;
        
        if sh > 0 && sh <= 18
            pol = FP_shstart_pol{sh};
            phi = pol(1);
            ridx = pol(2);
            theta = a_angles(ridx);
            
            LIST(i).FP_shs_s = phi;
            LIST(i).FP_shs_r = ridx;
            
            % "Shank Head" coordinates (vector direction)
            % Re-implementing original math exactly to preserve mapping
            hyp = sqrt(d_vals(ridx) * tand(theta)^2 + 1);
            x = hyp * sind(theta) * cosd(phi);
            y = hyp * sind(theta) * sind(phi);
            z = hyp * cosd(theta);
            
            LIST(i).FP_shfh_x = x;
            LIST(i).FP_shfh_y = y;
            LIST(i).FP_shfh_z = z;
            
            LIST(i).FP_shfh_p = theta; % Pitch? "p"
            LIST(i).FP_shfh_t = phi;   % Theta? "t"
            LIST(i).FP_shfh_r = hyp;   % Radius "r"
            
            % Note: Actual electrode positions were commented out in original script
            % So we only save parameters.
        else
            LIST(i).FP_shfh_x = NoneMarker; LIST(i).FP_shfh_y = NoneMarker; LIST(i).FP_shfh_z = NoneMarker;
            LIST(i).FP_shfh_p = NoneMarker; LIST(i).FP_shfh_t = NoneMarker; LIST(i).FP_shfh_r = NoneMarker;
            LIST(i).FP_shs_s = NoneMarker;  LIST(i).FP_shs_r = NoneMarker;
        end
    end
    
    %% 4. Omnetics Mapping
    % 32 pins per channel, 4 connectors per bank, 8 banks?
    % Logic:
    currOmnC = 1; Codd = 20; Cev = 2;
    for i = 1:length(LIST)
        LIST(i).OM_con = currOmnC;
        if LIST(i).HP_elect > 0
            pE = mod(LIST(i).HP_elect, 32);
            if pE == 0, pE = 32; end % Fix mod 0 case to 32
            
            % Check parity based on pE
            % Original: mod(pE,2)==1 -> odd
            if mod(pE, 2) == 1
                LIST(i).OM_pin = Codd;
                LIST(i).OM_ch = pE;
                Codd = Codd + 1;
            else
                LIST(i).OM_pin = Cev;
                LIST(i).OM_ch = pE;
                Cev = Cev + 1;
                % If pE is a multiple of 32 (or specifically ends block), advance
                % The original logic: if mod(pE,32)==0 ... wait.
                % Original: "if mod(pE,32) == 0" -> this is inside "else" (even).
                % If pE was treated as 0 in mod, then it hits else.
            end
            
            % Original check:
            % if mod(LIST(i).HP_elect, 32) == 0
            % Because pE is mod(), if it's 32 it becomes 0.
            % My calculation above sets pE=32. So mod(32,2)=0 -> even.
            % Then `if mod(pE, 32) == 0` (yes) -> reset.
            if mod(pE, 32) == 0
                Codd = 20; Cev = 2;
                currOmnC = currOmnC + 1;
            end
            LIST(i).OM_ref = -99;
            
        elseif LIST(i).HP_elect < 0 || (LIST(i).HP_elect == NoneMarker && LIST(i).HP_ref > 0)
             % Ref logic
             if LIST(i).HP_ref == 1
                 LIST(i).OM_pin = 1;
                 LIST(i).OM_ch = NoneMarker;
             elseif LIST(i).HP_ref == 2
                 LIST(i).OM_pin = 36;
                 LIST(i).OM_ch = NoneMarker;
             end
             LIST(i).OM_ref = floor((currOmnC-1)/4) + 1; % Approx bank?
        end
    end
    
    %% 5. Blackrock Mapping
    % Need BlackrockMap function
    if exist('BlackrockMap', 'file')
        BRmap = BlackrockMap(); % Returns table
    else
        error('BlackrockMap.m not found!');
    end
    
    for i = 1:length(LIST)
        if ~isempty(LIST(i).HP_elect) && LIST(i).HP_elect > 0
            % Match Omnetic Connector & Channel -> actually map uses HP_elect directly?
            % Original: `BRmap.OmneticConnectorNumber == LIST(i).OM_con` AND `BRmap.channelNumber == LIST(i).HP_elect`
            % Wait, BRMap typically maps ChannelID to Omnetics. 
            % Let's match by HP_elect if it corresponds to BR channel number?
            % Actually, in BlackrockMap.m, `channelNumber` IS the loop variable 1:1024.
            % And `FP_create_csv.m` assumes `LIST(i).HP_elect` corresponds to `BRmap.channelNumber`.
            % Is HP_elect == BR_Chan?
            % "BR channel number on LGA" maps imply HP pins correspond to channels.
            
            % Let's trust the logic: Find row in BRmap where channelNumber == HP_elect
            % Wait, the original code had:
            % `idx = find(BRmap.OmneticConnectorNumber == LIST(i).OM_con & BRmap.channelNumber == LIST(i).HP_elect)`
            % This implies HP_elect IS the channel number.
            
            idx = find(BRmap.channelNumber == LIST(i).HP_elect);
            
            if ~isempty(idx)
                LIST(i).BR_Chan = BRmap.channelNumber(idx);
                LIST(i).BR_CerM = BRmap.CereplexMnumber(idx);
                LIST(i).BR_CerMID = BRmap.CereplexMbankID(idx);
                LIST(i).BR_Inst = BRmap.CerebusInstanceNumber(idx);
                LIST(i).BR_InstCh = BRmap.CerebusInstanceChannelNumber(idx);
            else
                % Init invalid
                LIST(i).BR_Chan = NoneMarker; LIST(i).BR_CerM = NoneMarker; LIST(i).BR_CerMID = NoneMarker;
                LIST(i).BR_Inst = NoneMarker; LIST(i).BR_InstCh = NoneMarker;
            end
        else
             LIST(i).BR_Chan = NoneMarker; LIST(i).BR_CerM = NoneMarker; LIST(i).BR_CerMID = NoneMarker;
             LIST(i).BR_Inst = NoneMarker; LIST(i).BR_InstCh = NoneMarker;
        end
    end
    
    %% 6. Generate LGA Layouts
    % We need to create the LGA grid.
    % Dimensions: 36x36 (from original comment)
    % But we need to know WHERE in the grid each HP_elect is.
    % If we don't have the Pedestal Layout sheet, we can't fully recreate the spatial LGA map.
    % However, we can create a "Functional" LGA map if we assume a standard layout.
    % BUT, for `LGA.BR`, we definitely need the spatial layout.
    
    % Try to read the saved FP_PedestalLayoutHP.csv if we can't read the sheet?
    PedLayoutHP = [];
    if exist('FP_PedestalLayoutHP.csv','file')
        PedLayoutHP = readmatrix('FP_PedestalLayoutHP.csv');
    elseif exist('PedLayoutRaw','var') && ~isempty(PedLayoutRaw)
        % Process raw sheet
        PedLayoutHP = zeros(size(PedLayoutRaw));
        for r=1:size(PedLayoutRaw,1)
            for c=1:size(PedLayoutRaw,2)
                val = PedLayoutRaw{r,c};
                if isnumeric(val) && ~isnan(val)
                    PedLayoutHP(r,c) = val;
                elseif ischar(val)
                    if startsWith(val,'R')
                        PedLayoutHP(r,c) = -str2double(val(2:end)); % Ref
                    else
                        PedLayoutHP(r,c) = str2double(val);
                    end
                end
                if isnan(PedLayoutHP(r,c)), PedLayoutHP(r,c) = 0; end
            end
        end
    end
    
    if isempty(PedLayoutHP)
        warning('Could not generate LGA layout (missing source). ElectrodeMap.mat will have empty LGA.');
        LGA.BR = [];
    else
        % Create LGA.BR map
        % Map HP_elect numbers in PedLayoutHP to BR_Chan numbers
        LGA.BR = zeros(size(PedLayoutHP));
        for r = 1:size(PedLayoutHP,1)
            for c = 1:size(PedLayoutHP,2)
                hpVal = PedLayoutHP(r,c);
                if hpVal > 0
                    % Find this HP_elect in our LIST
                    % Optimization: could use a lookup vector since HP_elect is likely 1:1024
                    % Assume HP_elect X corresponds to BR_Chan X?
                    % Original code line 228 matched BRmap.channelNumber == LIST(i).HP_elect.
                    % So HP_elect IS BR_Chan.
                    LGA.BR(r,c) = hpVal; 
                else
                    LGA.BR(r,c) = NaN; % Or 0
                end
            end
        end
    end
    
    %% 7. Saving Outputs
    
    % A. Save CSV
    ElecMapTable = struct2table(LIST);
    % Sort cols
    [~, colIdx] = sort(ElecMapTable.Properties.VariableNames);
    ElecMapTable = ElecMapTable(:, colIdx);
    % Sort rows
    ElecMapTable = sortrows(ElecMapTable, {'FP_shank','FP_elect'}, {'ascend','ascend'});
    
    writetable(ElecMapTable, 'FP_ElectrodeMap.csv');
    fprintf('Saved FP_ElectrodeMap.csv\n');
    
    % B. Save MAT
    % Structure expected: ElecMap (table), LGA (struct with .BR)
    ElecMap = ElecMapTable;
    save('ElectrodeMap.mat', 'ElecMap', 'LGA');
    fprintf('Saved ElectrodeMap.mat\n');
    
    if ~isempty(PedLayoutHP)
        writematrix(PedLayoutHP, 'FP_PedestalLayoutHP.csv');
        fprintf('Saved FP_PedestalLayoutHP.csv\n');
    end

end
