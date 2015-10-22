function [D_final,velocity_final,plotdata] = RunDiffusionAnalysis(dataObject,trajectoryMethod)
% [D_final,velocity_final] = RunDiffusionAnalysis(dataObject,trajectoryMethod)
% Perform diffusion and velocity analysis on trajectories obtained from TrackNTrace, MOSAIC or other tracking programs according to all options defined in SetDefaultOptions.m
% 
% INTPUT:
%     dataObject: Can be a cell array of filename strings, a filename string, or a trajectory array.
%     If it's a cell array, the routine will try to load each file's trajectory data and combine all loaded trajectories into one array to analyze.
%     If it's a filename string, the same applies but the routine will also save the analysis to this file if possible.
%     If it's a trajectory array, the routine will assume that the array is formatted in the same way as if it was loaded from its respective file and then analyze it.
%     
%     trajectoryMethod: String of tracking method of choice. Can be 'TNT', 'MOSAIC','Icytrack' or 'utrack'. 
% 
% 
% OUTPUT:
%     Length and time scale of the output are [µm] and [s], respectively.
%     
%     D_final: Struct containing result of diffusion analysis. The content
%     will vary depending on the analyis methods chosen in
%     setDefaultOptions.m.
%     
%     velocity_final: Struct containing result of velocity analysis. The
%     content will vary depending on the analyis methods chosen in
%     setDefaultOptions.m.
% 
%     plotdata: Struct of histogram plot curves.

%% Initialize
fullPath = mfilename('fullpath');
[path,~,~] = fileparts(fullPath);
addpath(genpath([path,filesep,'subfun']));
addpath(genpath([path,filesep,'external']));
addpath(genpath([path,filesep,'helper']));

% Set options
[experimParam,trajParam,fitParam,printParam] = setDefaultOptions();

%% Calculate displacements/find trajectories
%no input? then ask user
if isempty(dataObject)
    [filename, pathname, ~] = uigetfile({'*.mat';'*.txt';'*.xls'}, 'Select all files to analyze.', 'MultiSelect','on');
    dataObject = cell(numel(filename),1);
    for iFiles=1:numel(dataObject)
        dataObject(iFiles) = {[pathname,filename{iFiles}]};
    end
    saveToFile = true;
    saveFilename = [pathname,'diffusionAnalysis.mat'];
    saveOverwrite = false;
else
    % otherwise, parse input
    if ischar(dataObject)
        dataObject = {dataObject};
        saveToFile = true;
        saveFilename = dataObject{1};
    else
        if ~iscell(dataObject)
            dataObject = {dataObject};
        end
        saveToFile = false;
    end
    saveOverwrite = true;
end

fprintf('\nProcessing tracks ... \n');
tracks_processed = cell(numel(dataObject),1);
for iFiles=1:numel(dataObject)
    [tracks_temp] = convertTrajectories(dataObject{iFiles},trajParam,trajectoryMethod);
    tracks_processed(iFiles) = {vertcat(tracks_temp{:})};
end

tracks_processed = vertcat(tracks_processed{:});
fprintf(repmat(sprintf('\b'), 1, 5));
fprintf('done. \n');

%% Calculate diffusion coefficients / MSD values
fprintf('Fitting diffusion coefficient ... \n');
switch fitParam.method
    case('histogram')
        [msd_result, velocity_final, plotdata,fitParam.fitVelocity] = diffusionAnalysisHistogram(tracks_processed,fitParam,experimParam);
        [D_final] = diffusionAnalysisFit(msd_result,experimParam,fitParam,printParam);
        
    case('msd')
        [msd_result,fitParam.fitVelocity] = diffusionAnalysisMSD(tracks_processed,fitParam,experimParam);
        [D_final,velocity_final] = diffusionAnalysisMSDFit(msd_result,experimParam,fitParam,printParam);
        
    case('msd-single')
        [D_final,velocity_final] = diffusionAnalysisMSDSingle(tracks_processed,fitParam,experimParam);
        
    otherwise
        fprintf('Fit method unknown, aborting. \n');
        return
end
fprintf('Fitting done. \n');
if ~exist('plotdata','var');
    plotdata = [];
end

if saveToFile
    diffusionAnalysis_result.D_final = D_final;
    diffusionAnalysis_result.velocity_final = velocity_final;
    diffusionAnalysis_result.plotdata = plotdata;
    diffusionAnalysis_result.experimParam = experimParam;
    diffusionAnalysis_result.trajParam = trajParam;
    diffusionAnalysis_result.fitParam = fitParam;
    diffusionAnalysis_result.printParam = printParam; %#ok<STRNU>
    
    if saveOverwrite
        save(saveFilename,'diffusionAnalysis_result','-append');
    else
        save(saveFilename,'diffusionAnalysis_result');
    end
end