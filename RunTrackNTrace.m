% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
%     Copyright (C) 2016  Simon Christoph Stein, scstein@phys.uni-goettingen.de
%     Copyright (C) 2016  Jan Thiart, jthiart@phys.uni-goettingen.de
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% 	CT, 2020:
% 	- support for import plugins
% 	- support for postprocessing plugins
%	- preset can be set as input arguments
%   - GUI-less mode
%
%
function [varargout] = RunTrackNTrace(varargin)
% Run this function to start the TrackNTrace application.

clearGlobals(); % Clear global variables used by TNT

% Global variables accessible by plugins
global globalOptions;
global importOptions;
global candidateOptions;
global refinementOptions;
global trackingOptions;
global postprocOptions;
global movie;
global filename_movie;
global imgCorrection; %#ok<NUSED>

addRequiredPathsTNT();

fprintf('Starting Track''N''Trace.\n')

%% Load default options
[globalOptions_def, TNToptions] = getDefaultOptions();

% check for overwrite arguments, eg. RunTrackNTrace('enableParallelProcessing',false)
TNToptions = cell2struct(cellfun(@(fname)parsevarargs(varargin,fname,TNToptions.(fname)),fieldnames(TNToptions),'UniformOutput',false),fieldnames(TNToptions));
globalOptions_def = cell2struct(cellfun(@(fname)parsevarargs(varargin,fname,globalOptions_def.(fname)),fieldnames(globalOptions_def),'UniformOutput',false),fieldnames(globalOptions_def));
%% Check if parallel processing is available
global parallelProcessingAvailable
parallelProcessingAvailable = false;

MATLAB_2013b_or_newer = false;
% Check MATLAB version
MATLABversion = strsplit(version,'.');
if(str2double(MATLABversion(1))>=8 && str2double(MATLABversion(2))>=2) % matlabpool -> parpool in MATLAB 2013b (8.2.x) and later
    MATLAB_2013b_or_newer = true; 
end


if TNToptions.enableParallelProcessing
    if MATLAB_2013b_or_newer
        try
            p = gcp('nocreate');
            if isempty(p)
                nrRunningWorkers = 0;
            else
                nrRunningWorkers = p.NumWorkers;
            end
            
            if(nrRunningWorkers == 0)
                parpool('local');
                p = gcp('nocreate');
                nrRunningWorkers = p.NumWorkers;
            end
            parallelProcessingAvailable = true;
            fprintf('TNT: Parallel processing available (%i workers).\n', nrRunningWorkers)
        catch
            parallelProcessingAvailable = false;
            fprintf('TNT: Parallel processing unavailable.\n')
        end
    else      
        try
            nrRunningWorkers = matlabpool('size'); %#ok<*DPOOL>
            if(nrRunningWorkers == 0)
                matlabpool('open','local');
            end
            parallelProcessingAvailable = true;
            fprintf('TNT: Parallel processing available (%i workers).\n', matlabpool('size'))
        catch
            parallelProcessingAvailable = false;
            fprintf('TNT: Parallel processing unavailable.\n')
        end
    end
end

if nargin==1 && iscellstr(varargin{1})
    %Load TNT file
    %TODO add checks for filetype and option to load files and show
    %settings GUI first
    list_filenames_TNTdata = varargin{1};
else 
    [import_plugins,accepted_formats] = loadImportPlugins();
    if nargin > 2 %eg passed from TNTvisulizer
        movie_list = cellstr(varargin{1}); % Can be one or multiple files
        globalOptions = parsevarargs(varargin(2:end),'globalOptions');
        importOptions = parsevarargs(varargin(2:end),'importOptions');
        candidateOptions_loaded = parsevarargs(varargin(2:end),'candidateOptions');
        refinementOptions_loaded = parsevarargs(varargin(2:end),'refinementOptions');
        trackingOptions_loaded = parsevarargs(varargin(2:end),'trackingOptions');
        postprocOptions_loaded = parsevarargs(varargin(2:end),'postprocOptions');
        
        candidateData_loaded = parsevarargs(varargin(2:end),'candidateData');
        refinementData_loaded = parsevarargs(varargin(2:end),'refinementData');
        trackingData_loaded = parsevarargs(varargin(2:end),'trackingData');
        postprocData_loaded = parsevarargs(varargin(2:end),'postprocData');
        movieSize_loaded = parsevarargs(varargin(2:end),'movieSize');
        firstFrame_lastFrame_loaded = parsevarargs(varargin(2:end),'firstFrame_lastFrame');
        outputPath_loaded = parsevarargs(varargin(2:end),'outputPath',-1);
        
        showSettingsGUI = parsevarargs(varargin(2:end),'showSettingsGUI',true);
        showStartupGUI = parsevarargs(varargin(2:end),'showStartupGUI',false); % This overrides all passed options beside the movie_list
    else% No inputs
        showSettingsGUI = true;
        showStartupGUI = true;
        movie_list = {};
    end
    if showStartupGUI
        %% Startup (loading input data information)
        [movie_list, globalOptions,importOptions,candidateOptions_loaded,refinementOptions_loaded,trackingOptions_loaded,postprocOptions_loaded, ...
            candidateData_loaded, refinementData_loaded, trackingData_loaded, postprocData_loaded, movieSize_loaded, firstFrame_lastFrame_loaded, metadata_loaded, outputPath_loaded, GUIreturns] = startupGUI(accepted_formats,movie_list);
        if GUIreturns.userExit
            userExitFunc();
            return;
        end
    end
    candidateOptions = candidateOptions_loaded;
    refinementOptions = refinementOptions_loaded;
    trackingOptions = trackingOptions_loaded;
    postprocOptions = postprocOptions_loaded;

    if(isempty(globalOptions))
        globalOptions = globalOptions_def;
    else
        % Ensure that all options exist. This enables compatibility with older
        % versions.
        globalOptions_temp = globalOptions;
        globalOptions = globalOptions_def;
        for j=fieldnames(globalOptions_temp)'
            globalOptions.(j{1}) = globalOptions_temp.(j{1});
        end
    end
    %% Main GUI: Adjust options for each movie and test settings if desired
    GUIinputs.TNToptions = TNToptions;
    GUIinputs.showStartupInformation = true;
    if (isfloat(outputPath_loaded))
        GUIinputs.outputFolderSameAsMovie = true;
        GUIinputs.outputFolder = '';
    else
        GUIinputs.outputFolderSameAsMovie = false;
        GUIinputs.outputFolder = outputPath_loaded;
    end
    GUIinputs.candidateOptions_loaded = candidateOptions_loaded;
    GUIinputs.refinementOptions_loaded = refinementOptions_loaded;
    GUIinputs.show_candidateData_fromFile_cbx  =                                                 ~isempty(candidateData_loaded);
    GUIinputs.show_refinementData_fromFile_cbx = (GUIinputs.show_candidateData_fromFile_cbx   && ~isempty(refinementData_loaded)); % candidateData must always be there
    GUIinputs.show_trackingData_fromFile_cbx   = (GUIinputs.show_refinementData_fromFile_cbx  && ~isempty(trackingData_loaded));   % candidateData + refinementData must always be there
    GUIinputs.show_postprocData_fromFile_cbx   = (GUIinputs.show_trackingData_fromFile_cbx    && ~isempty(postprocData_loaded));   % candidateData + refinementData + trackingData must always be there
    GUIinputs.use_loaded_candidateData = true;
    GUIinputs.use_loaded_refinementData = true;
    GUIinputs.use_loaded_trackingData = true;
    GUIinputs.use_loaded_postprocData = true;
    GUIinputs.firstFrame_lastFrame_loaded = firstFrame_lastFrame_loaded;

    % Calculate default dark image if given in default options
    dark_img_def = [];
    if(~isempty(globalOptions_def.filename_dark_movie))
        dark_img_def = CalculateDark(read_tiff(globalOptions_def.filename_dark_movie));
    end

    % Create timestamp for output files
    time = clock;
    timestamp = sprintf('%i-m%02i-d%02i-%02ih%02i',time (1),time(2),time(3),time(4),time(5));

    % Iterate through all movies in the list
    if showSettingsGUI
        GUIreturns.useSettingsForAll = false;
    else
        GUIreturns = GUIinputs;
        GUIreturns.useSettingsForAll = true;
    end
    
    list_filenames_TNTdata = cell(0);
    for iMovie=1:numel(movie_list)
        % Is this movie the last one to analyze?
        if(iMovie == numel(movie_list))
           GUIinputs.lastMovieInList = true;
        else
           GUIinputs.lastMovieInList = false;
        end

        filename_movie = movie_list{iMovie};
        [pathToMovie,filename,~] = fileparts(filename_movie);

        % Check if movie can be read
        selected_import_plugin = selectImportPlugin(filename_movie,import_plugins);

        GUIinputs.hasTCSPC = isfield(import_plugins(selected_import_plugin).info,'hasTCSPC')&&import_plugins(selected_import_plugin).info.hasTCSPC;
        if(~isempty(filename_movie)) || selected_import_plugin == 0
            % Test file import
            try
                import_plugins(selected_import_plugin).setOptions(importOptions);   % This ensures the options are valid options.
                if 0==selectImportPlugin(filename_movie,import_plugins(selected_import_plugin))
                    % This should never error. Calling selectImportPlugin again
                    % is necessary to call a plugins setFile function with the
                    % new options.
                    error('Invalid import plugin options.');
                end
                importOptions = import_plugins(selected_import_plugin).getOptions();% If importOptions was empty it is now initilsied with the defaults.
                [movie, metadata] = import_plugins(selected_import_plugin).mainFunc(importOptions,filename_movie,[1 0],1); % Read 0 frames. Movie is supposed to have the right dimensions
            catch err
                warning('Could not read movie ''%s''.\n  Error: %s',filename_movie,err.message);
                continue;
            end
        else
            continue;
        end

        % Set dark image to default for this movie
        dark_img = dark_img_def;

        % Does the user want to adjust the settings per movie?
        if not(  GUIreturns.useSettingsForAll )
            % Give file to adjust settings for to GUI
            GUIinputs.filename_movie = filename_movie;

            % Enable photon conversion for PTU files
            if GUIinputs.hasTCSPC
                globalOptions.photonBias = 0;
                globalOptions.photonGain = 1;
                globalOptions.photonSensitivity = 1;
                globalOptions.usePhotonConversion = true;
            elseif (globalOptions.photonBias == 0 && globalOptions.photonGain == 1 && globalOptions.photonSensitivity == 1)
                globalOptions.photonBias = globalOptions_def.photonBias;
                globalOptions.photonGain = globalOptions_def.photonGain;
                globalOptions.photonSensitivity = globalOptions_def.photonSensitivity;
                globalOptions.usePhotonConversion = globalOptions_def.usePhotonConversion;
            end

            [globalOptions,importOptions, candidateOptions,refinementOptions,trackingOptions,postprocOptions, GUIreturns] = settingsGUI(globalOptions,import_plugins(selected_import_plugin),candidateOptions,refinementOptions,trackingOptions,postprocOptions, GUIinputs);
            % Save output folder settings
            GUIinputs.outputFolderSameAsMovie = GUIreturns.outputFolderSameAsMovie;
            GUIinputs.outputFolder = GUIreturns.outputFolder;

            % Save use loaded data settings
            GUIinputs.use_loaded_candidateData = GUIreturns.use_loaded_candidateData;
            GUIinputs.use_loaded_refinementData = GUIreturns.use_loaded_refinementData;
            GUIinputs.use_loaded_trackingData = GUIreturns.use_loaded_trackingData;
            GUIinputs.use_loaded_postprocData = GUIreturns.use_loaded_postprocData;

            GUIinputs.showStartupInformation = false; % Only show this on first startup
            if GUIreturns.userExit
                userExitFunc(); % Cleanup
                return;
            end

            % Check if different dark movie was given
            if(~strcmp(globalOptions_def.filename_dark_movie, globalOptions.filename_dark_movie))
                if(~isempty(globalOptions.filename_dark_movie))
                    try
                        dark_img = double(globalOptions.binFrame) .* CalculateDark(read_tiff(globalOptions.filename_dark_movie));
                    catch err
                        error('Error when calculating dark image from movie ''%s''.\n  Error: %s',globalOptions.filename_dark_movie,err.message);
                    end
                end
            end

            % Set same settings for all remaining movies if user said so
            if GUIreturns.useSettingsForAll; GUIreturns.previewMode = false; end

            % If preview mode is enabled, analyze first X frames and show GUI
            if GUIreturns.previewMode
                first_run = true;
                filename_dark_movie = globalOptions.filename_dark_movie;
                while true
                    if not(first_run)
                        [globalOptions,importOptions, candidateOptions,refinementOptions,trackingOptions,postprocOptions, GUIreturns] = settingsGUI(globalOptions,import_plugins(selected_import_plugin), candidateOptions,refinementOptions,trackingOptions,postprocOptions, GUIinputs);
                        % Save output folder settings
                        GUIinputs.outputFolderSameAsMovie = GUIreturns.outputFolderSameAsMovie;
                        GUIinputs.outputFolder = GUIreturns.outputFolder;

                        % Save use loaded data settings
                        GUIinputs.use_loaded_candidateData = GUIreturns.use_loaded_candidateData;
                        GUIinputs.use_loaded_refinementData = GUIreturns.use_loaded_refinementData;
                        GUIinputs.use_loaded_trackingData = GUIreturns.use_loaded_trackingData;
                        GUIinputs.use_loaded_postprocData = GUIreturns.use_loaded_postprocData;

                        if GUIreturns.userExit
                            userExitFunc();
                            return;
                        end
                        if GUIreturns.useSettingsForAll; GUIreturns.previewMode = false; end %dont go through other movies anymore
                    end

                    if not(GUIreturns.previewMode); break; end % If preview mode was disabled by user in the settingsGUI
                    % Check if requested frame interval has changed -> re-read movie if neccessary
                    % If possible read FLIM. movie is than a cell array
                    % containing {counts,lifetimes}
                    if first_run || GUIreturns.previewIntervalChanged || GUIreturns.timegateChanged || GUIreturns.importOptionsChanged
                        movieArgs = {filename_movie, [globalOptions.firstFrameTesting, globalOptions.lastFrameTesting], globalOptions.binFrame,true};
                        if globalOptions.useTimegate
                            movieArgs = [movieArgs,{[globalOptions.tgStart globalOptions.tgEnd]}]; %#ok<AGROW>
                        end
                        [movie, metadata] = importOptions.mainFunc(importOptions,movieArgs{:});
                    end
                    % Check if different dark movie was given
                    if(~strcmp(filename_dark_movie, globalOptions.filename_dark_movie))
                        if(~isempty(globalOptions.filename_dark_movie))
                            try
                                dark_img = double(globalOptions.binFrame) .* CalculateDark(read_tiff(globalOptions.filename_dark_movie));
                            catch err
                                error('Error when calculating dark image from movie ''%s''.\n  Error: %s',globalOptions.filename_dark_movie,err.message);
                            end
                        end
                        filename_dark_movie = globalOptions.filename_dark_movie;
                    end

                    % IF: this is the first run perform all steps. ELSE: Reuse unchanged data from the last run
                    if first_run
                        [candidateData_preview,refinementData_preview, trackingData_preview, postprocData_preview, previewOptions] = runPreview(movie,dark_img,metadata);
                    else
                        [candidateData_preview,refinementData_preview, trackingData_preview, postprocData_preview, previewOptions] = runPreview(movie,dark_img,metadata, candidateData_preview, refinementData_preview, trackingData_preview, postprocData_preview, previewOptions, GUIreturns);
                    end
                    first_run = false;
                end
            end

        end %not( GUIreturns.useSettingsForAll )


        % --  Save options   --

        % Build filename for TNT file
        if GUIreturns.outputFolderSameAsMovie
            TNToutputPath = [pathToMovie,filesep];
        else
            if isempty(GUIreturns.outputFolder)
                TNToutputPath = [];
            else            
                TNToutputPath = [GUIreturns.outputFolder,filesep];
            end
        end

        filename_TNTdata = [TNToutputPath,filename,'_',timestamp,'_TNT.mat'];
        % prevent overriding filenames, e.g. when one movie is chosen multiple times
        fileID = 1;
        while(true)
           if exist(filename_TNTdata,'file') % Test if output file already exists
               filename_TNTdata = [TNToutputPath,filename,'_',timestamp,'_',int2str(fileID),'_TNT.mat'];
               fileID = fileID+1;
           else
               break;
           end
        end
        list_filenames_TNTdata = [list_filenames_TNTdata; {filename_TNTdata}]; %#ok<AGROW> % Append name of datafile to list

        % Should loaded xData be used?
        if(GUIreturns.use_loaded_candidateData)
            candidateOptions =  candidateOptions_loaded;
            struct_helper.candidateData = candidateData_loaded;
        end

        if(GUIreturns.use_loaded_refinementData)
            refinementOptions =  refinementOptions_loaded;
            struct_helper.refinementData = refinementData_loaded;
        end

        if(GUIreturns.use_loaded_trackingData)
            trackingOptions =  trackingOptions_loaded;
            struct_helper.trackingData = trackingData_loaded;
        end

        if(GUIreturns.use_loaded_postprocData)
            postprocOptions =  postprocOptions_loaded;
            struct_helper.postprocData = postprocData_loaded;
        end

        % Turn preview mode off. Toggles frame range in plugin_fitLT.
        globalOptions.previewMode = false;

        % Save the globalOptions and Options of the enabled steps
        savevars = {'filename_movie','globalOptions','importOptions','dark_img'};
        if(globalOptions.enableCandidate)
            savevars = [savevars,'candidateOptions']; %#ok<AGROW>
        end
        if(globalOptions.enableRefinement)
            savevars = [savevars,'refinementOptions']; %#ok<AGROW>
        end
        if(globalOptions.enableTracking)
            savevars = [savevars,'trackingOptions']; %#ok<AGROW>
        end
        if(globalOptions.enablePostproc)
            savevars = [savevars,'postprocOptions']; %#ok<AGROW>
        end
        save(filename_TNTdata,savevars{:},'-v7.3');

        % Append loaded data and movieSize + firstFrame_lastFrame so that movie
        % must not be loaded to do tracking
        if(GUIreturns.use_loaded_candidateData || GUIreturns.use_loaded_refinementData || GUIreturns.use_loaded_trackingData || GUIreturns.use_loaded_postprocData)
            struct_helper.movieSize = movieSize_loaded;
            struct_helper.firstFrame_lastFrame = firstFrame_lastFrame_loaded;
            metadata =  metadata_loaded;
            struct_helper.metadata = metadata_loaded;
            save(filename_TNTdata,'-append','-struct','struct_helper');
        end    

        if not(GUIreturns.useSettingsForAll || TNToptions.rememberSettingsForNextMovie)
            globalOptions    = [];
            candidateOptions = [];
            refinementOptions   = [];
            trackingOptions  = [];
            postprocOptions  = [];
        else
            dark_img_def = dark_img;
        end
    end
end
%% Clear all vars which are not needed any more.
clearvars -except list_filenames_TNTdata parallelProcessingAvailable TNToptions % This preserves the global attribute, eventhough the vars are cleared.
%% Candidate detection and refinement for every movie
for iMovie=1:numel(list_filenames_TNTdata)
    filename_TNTdata = list_filenames_TNTdata{iMovie};
    load(filename_TNTdata,'-mat'); %#ok<LOAD>
    
    if not(exist('candidateData','var') && exist('refinementData','var')) || isempty(candidateData) || isempty(refinementData)
        % Read movie (intensity only)
        if iMovie==1 || ~isequaln(parameters_movie_last_loop,{filename_movie,importOptions,[globalOptions.firstFrame,globalOptions.lastFrame],globalOptions.binFrame,globalOptions.useTimegate,[globalOptions.tgStart globalOptions.tgEnd]})
            movieArgs = {filename_movie, [globalOptions.firstFrame, globalOptions.lastFrame], globalOptions.binFrame,true};
            if globalOptions.useTimegate
                movieArgs = [movieArgs,{[globalOptions.tgStart globalOptions.tgEnd]}]; %#ok<AGROW>
            end
            clear movie metadata;%Free up memory
            [movie,metadata] = importOptions.mainFunc(importOptions,movieArgs{:});
        end
        parameters_movie_last_loop = {filename_movie,importOptions,[globalOptions.firstFrame,globalOptions.lastFrame],globalOptions.binFrame,globalOptions.useTimegate,[globalOptions.tgStart globalOptions.tgEnd]};

        % Compute the positions
        if globalOptions.enableCandidate
            if not(exist('candidateData','var')) || isempty(candidateData)
                fprintf('######\nTNT: Locating candidates in movie %s.\n',filename_movie);
                [candidateData, candidateOptions] = findCandidateParticles(movie{1}, dark_img, globalOptions, candidateOptions);
            else
                fprintf('######\nTNT: Using loaded candidateData for processing.\n');
            end
        else
            candidateData = [];
        end
        if globalOptions.enableRefinement
            if not(exist('refinementData','var')) || isempty(refinementData)
                fprintf('######\nTNT: Refining positions in movie %s.\n',filename_movie);
                [refinementData, refinementOptions] = fitParticles(movie{1}, dark_img, globalOptions, refinementOptions, candidateData);
            else
                fprintf('######\nTNT: Using loaded refinementData for processing.\n');
            end
        else
            refinementData = [];
        end

        % Save positions, movieSize, and index of first and last frame processed
        firstFrame_lastFrame = [globalOptions.firstFrame,  globalOptions.firstFrame + size(movie{1},3)-1];  % % Note: lastFrame could have been set to 'inf', now we synchronize with the correct number
        movieSize = size(movie{1}); % % Save size of movie (nice to have)
        
        save(filename_TNTdata,'metadata','-append'); % For backwards compatibility. Metadata is only added if it is generated from the importPlugin
    else
        fprintf('######\nTNT: Using loaded candidateData & refinementData for processing.\n');
    end
    
    save(filename_TNTdata,'candidateData','refinementData','globalOptions','candidateOptions','refinementOptions','movieSize','firstFrame_lastFrame','-append');
    candidateData = []; refinementData = [];
end

%% Compute trajectories for every movie
for iMovie=1:numel(list_filenames_TNTdata)
    
    % Load global options to check if tracking is desired (skip movie if it is not)
    load(list_filenames_TNTdata{iMovie},'globalOptions','filename_movie');
    if (~globalOptions.enableTracking)
        continue
    end
    
    
    hasTrackingData = whos('-file',list_filenames_TNTdata{iMovie},'trackingData');
    hasTrackingData = ~isempty(hasTrackingData) && prod(hasTrackingData.size)>0;
    hasPostprocData = whos('-file',list_filenames_TNTdata{iMovie},'postprocData');
    hasPostprocData = ~isempty(hasPostprocData) && prod(hasPostprocData.size)>0;
    if ~hasTrackingData
        % Load options and data needed for processing
        load(list_filenames_TNTdata{iMovie},'trackingOptions','refinementData');
        % Compute trajectories
        fprintf('######\nTNT: Tracking particles in movie %s.\n',filename_movie);
        [trackingData, trackingOptions] = trackParticles(refinementData,trackingOptions);
        
        clear refinementData;        
        %Save trajectories
        save(list_filenames_TNTdata{iMovie},'trackingData','trackingOptions','-append');
    elseif globalOptions.enablePostproc && ~hasPostprocData
        load(list_filenames_TNTdata{iMovie},'trackingData');        
    end
    if (~globalOptions.enablePostproc) || hasPostprocData
        continue
    end
    
    % Load options and data needed for processing
    load(list_filenames_TNTdata{iMovie},'postprocOptions');
    
    % Compute trajectories
    fprintf('######\nTNT: Postproccessing tracks in movie %s.\n',filename_movie);
    [postprocData, postprocOptions] = postprocTracks(trackingData,globalOptions,importOptions,postprocOptions);
    
    %Save trajectories
    save(list_filenames_TNTdata{iMovie},'postprocData','postprocOptions','-append');
end

% Deactivate parallel processing if requested
if parallelProcessingAvailable && TNToptions.closeMatlabpoolOnExit
    matlabpool('close');
end

% Return TNT files if requsted
if nargout>0
    varargout{1} = list_filenames_TNTdata;
end

% Clear globals
clearGlobals();

%% Add required folders and subfolders to path
    function addRequiredPathsTNT()
        fullPathToThisFile = mfilename('fullpath');
        [path,~,~] = fileparts(fullPathToThisFile);
        addpath(genpath([path,filesep,'external']));
        addpath(genpath([path,filesep,'helper']));
        addpath(genpath([path,filesep,'plugins']));
        addpath(genpath([path,filesep,'subfun']));
        addpath(genpath([path,filesep,'analysis']));
    end

%% Cleanup function if something goes wrong
    function userExitFunc()
        % Remove TNTdata files
        if exist('list_filenames_TNTdata','var')
            fprintf('TNT: User exit.\n')
%             warning off backtrace
%             warning('User exit. Stopping TrackNTrace. Deleting settings files that have been saved already.');
%             warning on backtrace
            for iTNTfile = 1:numel(list_filenames_TNTdata)
                if exist(list_filenames_TNTdata{iTNTfile},'file')
                    delete(list_filenames_TNTdata{iTNTfile});
                end
            end
            if nargout>0
                varargout{1} = {};
            end
        end
        
        if parallelProcessingAvailable && TNToptions.closeMatlabpoolOnExit
            matlabpool('close');
        end
        clearGlobals();
    end

% Clear all global variables
    function clearGlobals()
        clear global globalOptions importOptions candidateOptions refinementOptions trackingOptions postrprocOptions movie imgCorrection parallelProcessingAvailable filename_movie;
    end
    
    function [value,isfallback] = parsevarargs(args,arg,fallback)
        idx = find(strcmpi(args(1:end-1),arg),1,'last');
        isfallback = isempty(idx);
        if isfallback
            if nargin<3
                value = [];
            else
                value = fallback;
            end
        else
            value = args{idx+1};
        end
    end

end

