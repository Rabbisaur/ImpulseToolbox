function plotFountainProbeData(dataVector, varargin)
% PLOTFOUNTAINPROBEDATA Visualizes 1024-channel data on LGA and Fountain Probe geometry.
%
% Usage:
%   plotFountainProbeData(dataVector)
%   plotFountainProbeData(dataVector, 'Name', 'Value', ...)
%
% Inputs:
%   dataVector - A 1x1024 or 1024x1 numerical vector representing values for
%                Blackrock channels 1-1024.
%
% Optional Parameters (Name-Value pairs):
%   'MapFile'     - Path to ElectrodeMap.mat (default: looks in path)
%   'CSVFile'     - Path to FP_ElectrodeMap.csv (default: looks in path)
%   'FigureTitle' - Title for the figures (default: 'Fountain Probe Data')
%   'CLim'        - Color limits [min max] (default: auto)
%   'PlotLGA'     - Boolean, plot LGA map (default: true)
%   'PlotShanks'  - Boolean, plot Shank map (default: true)
%
% Example:
%   data = rand(1, 1024);
%   plotFountainProbeData(data, 'FigureTitle', 'Noise Levels');

    p = inputParser;
    addRequired(p, 'dataVector', @(x) isnumeric(x) && numel(x) == 1024);
    addParameter(p, 'MapFile', 'ElectrodeMap.mat', @ischar);
    addParameter(p, 'CSVFile', 'FP_ElectrodeMap.csv', @ischar);
    addParameter(p, 'FigureTitle', 'Fountain Probe Data', @ischar);
    addParameter(p, 'CLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    addParameter(p, 'PlotLGA', true, @islogical);
    addParameter(p, 'PlotShanks', true, @islogical);
    parse(p, dataVector, varargin{:});
    
    dataVector = p.Results.dataVector(:); % Ensure column vector
    mapFile = p.Results.MapFile;
    csvFile = p.Results.CSVFile;
    figTitle = p.Results.FigureTitle;
    cLim = p.Results.CLim;
    
    %% 1. Plot LGA Map
    if p.Results.PlotLGA
        plotLGA(dataVector, mapFile, figTitle, cLim);
    end
    
    %% 2. Plot Fountain Probe Shanks
    if p.Results.PlotShanks
        plotShanks(dataVector, csvFile, figTitle, cLim);
    end
end

function plotLGA(data, mapFile, titleStr, cLim)
    % Load LGA Map
    if ~exist(mapFile, 'file')
        % Try to find it
        found = dir(fullfile('**', mapFile));
        if isempty(found)
            warning('ElectrodeMap.mat not found. Skipping LGA plot.');
            return;
        end
        mapFile = fullfile(found(1).folder, found(1).name);
    end
    
    loaded = load(mapFile);
    if ~isfield(loaded, 'LGA') || ~isfield(loaded.LGA, 'BR')
        warning('ElectrodeMap.mat does not contain LGA.BR. Skipping LGA plot.');
        return;
    end
    
    lgaGrid = loaded.LGA.BR; % Matrix where values are BR Channel IDs
    
    % Map data to grid
    plotGrid = nan(size(lgaGrid));
    validChans = lgaGrid > 0 & lgaGrid <= 1024 & ~isnan(lgaGrid);
    plotGrid(validChans) = data(lgaGrid(validChans));
    
    % Plot
    figure('Name', [titleStr ' - LGA'], 'Color', 'w');
    imagesc(plotGrid);
    colorbar;
    if ~isempty(cLim), caxis(cLim); end
    title([titleStr ' - LGA Layout']);
    axis image;
    xlabel('Column'); ylabel('Row');
    
    % Add text for orientation if known, or just leave as grid
end

function plotShanks(data, csvFile, titleStr, cLim)
    % Load CSV Map
    if ~exist(csvFile, 'file')
        % Try to find it
        found = dir(fullfile('**', csvFile));
        if isempty(found)
            warning('FP_ElectrodeMap.csv not found. Skipping Shank plot.');
            return;
        end
        csvFile = fullfile(found(1).folder, found(1).name);
    end
    
    try
        T = readtable(csvFile);
    catch
        warning('Failed to read FP_ElectrodeMap.csv. Skipping Shank plot.');
        return;
    end
    
    % Extract valid rows (where BR_Chan is valid)
    validRows = T.BR_Chan > 0 & T.BR_Chan <= 1024 & ~isnan(T.BR_Chan);
    T = T(validRows, :);
    
    % Geometry: 
    % We want to plot 18 shanks.
    % We can arrange them in a 2D layout:
    % X-axis: Shank ID (1-18)
    % Y-axis: Electrode Depth (FP_elect 1-56)
    
    % Create grid: 56 depths x 18 shanks
    shankGrid = nan(56, 18);
    
    for i = 1:height(T)
        brCh = T.BR_Chan(i);
        shank = T.FP_shank(i);
        elec = T.FP_elect(i);
        
        if shank >= 1 && shank <= 18 && elec >= 1 && elec <= 56
            shankGrid(elec, shank) = data(brCh);
        end
    end
    
    % Plot flattened view
    figure('Name', [titleStr ' - Shanks Flat'], 'Color', 'w');
    imagesc(shankGrid);
    colorbar;
    if ~isempty(cLim), caxis(cLim); end
    title([titleStr ' - Shank Depth View']);
    xlabel('Shank ID (1-18)'); 
    ylabel('Electrode ID (1=Deep, 56=Shallow)');
    set(gca, 'XTick', 1:18);
    axis xy; % Normal Y axis (1 at bottom) -> Wait, 1 is usually deepest? 
    % README: "lower numbers are deeper in the brain" -> So 1 is Deep.
    % If we want Deep at bottom, 'axis xy' puts 1 at bottom. Correct.
    
    
    % Optional: 3D "Fountain" Plot
    % If we have spatial coordinates
    if ismember('FP_shfh_x', T.Properties.VariableNames)
        figure('Name', [titleStr ' - 3D Fountain'], 'Color', 'w');
        scatter3(T.FP_shfh_x, T.FP_shfh_y, T.FP_shfh_z, 20, data(T.BR_Chan), 'filled');
        colorbar;
        if ~isempty(cLim), caxis(cLim); end
        title([titleStr ' - 3D Fountain Geometry']);
        axis equal;
        rotate3d on;
        xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
        view(3);
    end
end
