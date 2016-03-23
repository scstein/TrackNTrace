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
function RunTrackNTrace()
% Run this function to start the TrackNTrace application.

clearGlobals(); % Clear global variables used by TNT

% Global variables accessible by plugins
global globalOptions;
global candidateOptions
global refinementOptions
global trackingOptions
global movie;
global filename_movie;
global imgCorrection; %#ok<NUSED>

addRequiredPathsTNT();

fprintf('Starting Track''N''Trace.\n')

%% Load default options
[globalOptions_def, TNToptions] = getDefaultOptions();

%% Check if parallel processing is available
global parallelProcessingAvailable
parallelProcessingAvailable = false;

if TNToptions.enableParallelProcessing
    try
        nrRunningWorkers = matlabpool('size');
        if(nrRunningWorkers == 0);
            matlabpool('open','local');
        end
        parallelProcessingAvailable = true;
        fprintf('TNT: Parallel processing available (%i workers).\n', matlabpool('size'))
    catch
        parallelProcessingAvailable = false;
        fprintf('TNT: Parallel processing unavailable.\n')
    end
end

%% Startup (loading input data information)
[movie_list, globalOptions, candidateOptions_loaded,refinementOptions_loaded,trackingOptions_loaded, ...
    candidateData_loaded, refinementData_loaded, movieSize_loaded, firstFrame_lastFrame_loaded, outputPath_loaded, GUIreturns] = startupGUI();
candidateOptions = candidateOptions_loaded;
refinementOptions = refinementOptions_loaded;
trackingOptions = trackingOptions_loaded;

if(isempty(globalOptions))
    globalOptions = globalOptions_def;
end;
if GUIreturns.userExit;
    userExitFunc();
    return;
end;

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
GUIinputs.show_candidateData_fromFile_cbx = ~isempty(candidateData_loaded);
GUIinputs.show_refinementData_fromFile_cbx = (GUIinputs.show_candidateData_fromFile_cbx && ~isempty(refinementData_loaded)); % candidateData must always be there
GUIinputs.use_loaded_candidateData = true;
GUIinputs.use_loaded_refinementData = true;
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
GUIreturns.useSettingsForAll = false;

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
    if(~isempty(filename_movie))
        try
            [~,movie] = evalc(['read_tiff(''',filename_movie,''',false,[1,2])']); % Read 2 frames. note: evalc suppresses output
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
        
        [globalOptions, candidateOptions,refinementOptions,trackingOptions, GUIreturns] = settingsGUI(globalOptions, candidateOptions,refinementOptions,trackingOptions, GUIinputs);
        % Save output folder settings
        GUIinputs.outputFolderSameAsMovie = GUIreturns.outputFolderSameAsMovie;
        GUIinputs.outputFolder = GUIreturns.outputFolder;
        
        % Save use loaded data settings
        GUIinputs.use_loaded_candidateData = GUIreturns.use_loaded_candidateData;
        GUIinputs.use_loaded_refinementData = GUIreturns.use_loaded_refinementData;
        
        GUIinputs.showStartupInformation = false; % Only show this on first startup
        if GUIreturns.userExit;
            userExitFunc(); % Cleanup
            return;
        end;
        
        % Check if different dark movie was given
        if(~strcmp(globalOptions_def.filename_dark_movie, globalOptions.filename_dark_movie))
            if(~isempty(globalOptions.filename_dark_movie))
                try
                    dark_img = CalculateDark(read_tiff(globalOptions.filename_dark_movie));
                catch err
                    error('Error when calculating dark image from movie ''%s''.\n  Error: %s',globalOptions.filename_dark_movie,err.message);
                end
            end
        end
        
        % Set same settings for all remaining movies if user said so
        if GUIreturns.useSettingsForAll; GUIreturns.previewMode = false; end;
        
        % If preview mode is enabled, analyze first X frames and show GUI
        if GUIreturns.previewMode
            first_run = true;
            filename_dark_movie = globalOptions.filename_dark_movie;
            while true
                if not(first_run)
                    [globalOptions, candidateOptions,refinementOptions,trackingOptions, GUIreturns] = settingsGUI(globalOptions, candidateOptions,refinementOptions,trackingOptions, GUIinputs);
                    % Save output folder settings
                    GUIinputs.outputFolderSameAsMovie = GUIreturns.outputFolderSameAsMovie;
                    GUIinputs.outputFolder = GUIreturns.outputFolder;
                    
                    % Save use loaded data settings
                    GUIinputs.use_loaded_candidateData = GUIreturns.use_loaded_candidateData;
                    GUIinputs.use_loaded_refinementData = GUIreturns.use_loaded_refinementData;
                    
                    if GUIreturns.userExit;
                        userExitFunc();
                        return;
                    end;
                    if GUIreturns.useSettingsForAll; GUIreturns.previewMode = false; end; %dont go through other movies anymore
                end
                
                if not(GUIreturns.previewMode); break; end; % If preview mode was disabled by user in the settingsGUI
                % Check if requested frame interval has changed -> re-read movie if neccessary
                if first_run || GUIreturns.previewIntervalChanged
                    movie = read_tiff(filename_movie, false, [globalOptions.firstFrameTesting, globalOptions.lastFrameTesting]);
                end
                % Check if different dark movie was given
                if(~strcmp(filename_dark_movie, globalOptions.filename_dark_movie))
                    if(~isempty(globalOptions.filename_dark_movie))
                        try
                            dark_img = CalculateDark(read_tiff(globalOptions.filename_dark_movie));
                        catch err
                            error('Error when calculating dark image from movie ''%s''.\n  Error: %s',globalOptions.filename_dark_movie,err.message);
                        end
                    end
                    filename_dark_movie = globalOptions.filename_dark_movie;
                end
                
                % IF: this is the first run perform all steps. ELSE: Reuse unchanged data from the last run
                if first_run
                    [candidateData_preview,refinementData_preview, trackingData_preview, previewOptions] = runPreview(movie,dark_img);
                else
                    [candidateData_preview,refinementData_preview, trackingData_preview, previewOptions] = runPreview(movie,dark_img, candidateData_preview, refinementData_preview, trackingData_preview, previewOptions, GUIreturns);
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
    
    % Should loaded candidateData be used?
    if(GUIreturns.use_loaded_candidateData)
        candidateOptions =  candidateOptions_loaded;
        struct_helper.candidateData = candidateData_loaded;
    end
    
    % Should loaded candidateData be used?
    if(GUIreturns.use_loaded_refinementData)
        refinementOptions =  refinementOptions_loaded;
        struct_helper.refinementData = refinementData_loaded;
    end
    
    save(filename_TNTdata,'filename_movie','globalOptions','candidateOptions','refinementOptions','dark_img');
    if(globalOptions.enableTracking) % Save tracking options only if tracking is desired
        save(filename_TNTdata,'trackingOptions','-append');
    end
    
    % Append loaded data and movieSize + firstFrame_lastFrame so that movie
    % must not be loaded to do tracking
    if(GUIreturns.use_loaded_candidateData || GUIreturns.use_loaded_refinementData)
        struct_helper.movieSize = movieSize_loaded;
        struct_helper.firstFrame_lastFrame = firstFrame_lastFrame_loaded;
        save(filename_TNTdata,'-append','-struct','struct_helper');
    end    
    
    if not(GUIreturns.useSettingsForAll || TNToptions.rememberSettingsForNextMovie)
        globalOptions    = [];
        candidateOptions = [];
        refinementOptions   = [];
        trackingOptions  = [];
    else
        dark_img_def = dark_img;
    end
end
clearvars -except list_filenames_TNTdata parallelProcessingAvailable TNToptions

%% Candidate detection and refinement for every movie
for iMovie=1:numel(list_filenames_TNTdata)
    filename_TNTdata = list_filenames_TNTdata{iMovie};
    load(filename_TNTdata,'-mat');
        
    if  not(exist('candidateData','var') && exist('refinementData','var'))
        % Read movie
        if iMovie==1 || ~strcmp(filename_movie,filename_movie_last_loop)
            movie = read_tiff(filename_movie, false, [globalOptions.firstFrame,globalOptions.lastFrame]);
        end
        filename_movie_last_loop = filename_movie;

        % Compute the positions
        if not(exist('candidateData','var'))
            fprintf('######\nTNT: Locating candidates in movie %s.\n',filename_movie);
            [candidateData, candidateOptions] = findCandidateParticles(movie, dark_img, globalOptions, candidateOptions);
        else
            fprintf('######\nTNT: Using loaded candidateData for processing.\n');
        end
        if not(exist('refinementData','var'))
            fprintf('######\nTNT: Refining positions in movie %s.\n',filename_movie);
            [refinementData, refinementOptions] = fitParticles(movie, dark_img, globalOptions, refinementOptions, candidateData);
        else
            fprintf('######\nTNT: Using loaded refinementData for processing.\n');
        end

        % Save positions, movieSize, and index of first and last frame processed
        firstFrame_lastFrame = [globalOptions.firstFrame,  globalOptions.firstFrame + size(movie,3)-1];  %#ok<NASGU> % Note: lastFrame could have been set to 'inf', now we synchronize with the correct number
        movieSize = size(movie); %#ok<NASGU> % Save size of movie (nice to have)
    else
        fprintf('######\nTNT: Using loaded candidateData & refinementData for processing.\n');
    end
    
    save(filename_TNTdata,'candidateData','refinementData','globalOptions','candidateOptions','refinementOptions','movieSize','firstFrame_lastFrame','-append');
    clearvars -except list_filenames_TNTdata parallelProcessingAvailable TNToptions globalOptions filename_movie_last_loop
end

%% Compute trajectories for every movie
for iMovie=1:numel(list_filenames_TNTdata)
    
    % Load global options to check if tracking is desired (skip movie if it is not)
    load(list_filenames_TNTdata{iMovie},'globalOptions');
    if (~globalOptions.enableTracking)
        continue
    end
    
    % Load options and data needed for processing
    load(list_filenames_TNTdata{iMovie},'trackingOptions','refinementData','filename_movie');
    
    % Compute trajectories
    fprintf('######\nTNT: Tracking particles in movie %s.\n',filename_movie);
    [trackingData, trackingOptions] = trackParticles(refinementData,trackingOptions); %#ok<ASGLU>
    
    %Save trajectories
    save(list_filenames_TNTdata{iMovie},'trackingData','trackingOptions','-append');
end

% Deactivate parallel processing if requested
if parallelProcessingAvailable && TNToptions.closeMatlabpoolOnExit
    matlabpool('close');
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
        end
        
        if parallelProcessingAvailable && TNToptions.closeMatlabpoolOnExit
            matlabpool('close');
        end
        clearGlobals();
    end

% Clear all global variables
    function clearGlobals()
        clear global globalOptions candidateOptions refinementOptions trackingOptions movie imgCorrection parallelProcessingAvailable filename_movie;
    end
end


