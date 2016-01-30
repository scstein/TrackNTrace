function RunTrackNTrace()
clearGlobals(); % Clear global variables used by TNT

% Global variables accessible by plugins
global globalOptions;
global candidateOptions
global fittingOptions
global trackingOptions
global movie;
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
            matlabpool('open');
        end
        parallelProcessingAvailable = true;
        fprintf('TNT: Parallel processing available (%i workers).\n', matlabpool('size'))
    catch
        parallelProcessingAvailable = false;
        fprintf('TNT: Parallel processing unavailable.\n')
    end
end

%% Load and adjust the default settings for this batch
[movie_list, globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIreturns] = startupGUI();
if(isempty(globalOptions))
    globalOptions = globalOptions_def;
end;
if GUIreturns.userExit;
    exitFunc();
    return;
end;

%% Adjust options for each movie and test settings if desired
GUIinputs.TNToptions = TNToptions;
GUIinputs.showStartupInformation = true;

% Calculate default dark image if given
dark_img_def = [];
if(~isempty(globalOptions_def.filename_dark_movie))
    dark_img_def = CalculateDark(read_tiff(globalOptions_def.filename_dark_movie));
end

% Create timestamp for output files
time = clock;
timestamp = sprintf('%i-m%02i-d%02i-%02ih%i',time(1),time(2),time(3),time(4),time(5));

% Iterate through all movies in the list
GUIreturns.useSettingsForAll = false;

posFit_list = cell(0);
list_filenames_TNTdata = cell(0);
for iMovie=1:numel(movie_list)
    filename_movie = movie_list{iMovie};
    [path,filename,~] = fileparts(filename_movie);
    filename_TNTdata = [path,filesep,filename,'_',timestamp,'_TNT.mat'];
    list_filenames_TNTdata = [list_filenames_TNTdata; {filename_TNTdata}]; % Append name of datafile to list
    
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
        
        [globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIreturns] = settingsGUI(globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs);
        GUIinputs.showStartupInformation = false; % Only show this on first startup
        if GUIreturns.userExit;
            exitFunc(); % Cleanup
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
                    [globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIreturns] = settingsGUI(globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs);
                    if GUIreturns.userExit;
                        exitFunc();
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
                    [candidateData_preview,fittingData_preview, trackingData_preview, previewOptions] = runPreview(movie,dark_img);
                else
                    [candidateData_preview,fittingData_preview, trackingData_preview, previewOptions] = runPreview(movie,dark_img, candidateData_preview, fittingData_preview, trackingData_preview, previewOptions, GUIreturns);
                end
                first_run = false;
            end
        end
        
    end %not( GUIreturns.useSettingsForAll )
        
    % Save options    
    save(filename_TNTdata,'filename_movie','globalOptions','candidateOptions','fittingOptions','dark_img');
    if(globalOptions.enableTracking) % Save tracking options only if tracking is desired
        save(filename_TNTdata,'trackingOptions','-append');
    end
    posFit_list = [posFit_list;{filename_TNTdata}]; %#ok<AGROW>
    
    
    if not(GUIreturns.useSettingsForAll || TNToptions.rememberSettingsForNextMovie)
        globalOptions    = [];
        candidateOptions = [];
        fittingOptions   = [];
        trackingOptions  = [];
    else
        dark_img_def = dark_img;
    end
end
clearvars -except posFit_list

%% Candidate detection and fitting for every movie
for iMovie=1:numel(posFit_list)
    filename_TNTdata = posFit_list{iMovie};
    
    load(filename_TNTdata,'-mat');
    
    % Read movie
    movie = read_tiff(filename_movie, false, [globalOptions.firstFrame,globalOptions.lastFrame]);
    % Compute the positions
    fprintf('######\nTNT: Locating particles in movie %s.\n',filename_movie);
    [candidateData, candidateOptions] = findCandidateParticles(movie, dark_img, globalOptions, candidateOptions);
    [fittingData, fittingOptions] = fitParticles(movie, dark_img, globalOptions, fittingOptions, candidateData);
    
    % Save positions, movieSize, and index of first and last frame processed
    firstFrame_lastFrame = [globalOptions.firstFrame,  globalOptions.firstFrame + size(movie,3)-1];  %#ok<NASGU> % Note: lastFrame could have been set to 'inf', now we synchronize with the correct number
    movieSize = size(movie); %#ok<NASGU> % Save size of movie (nice to have)
    save(filename_TNTdata,'candidateData','fittingData','globalOptions','candidateOptions','fittingOptions','movieSize','firstFrame_lastFrame','-append');
end
clearvars -except posFit_list

%% Compute trajectories for every movie
for iMovie=1:numel(posFit_list)
    
    % Load global options to check if tracking is desired (skip movie if it is not)
    load(posFit_list{iMovie},'globalOptions');
    if (~globalOptions.enableTracking)
        continue
    end
    
    % Load options and data needed for processing
    load(posFit_list{iMovie},'trackingOptions','fittingData','filename_movie');
    
    % Compute trajectories
    fprintf('######\nTNT: Tracking particles in movie %s.\n',filename_movie);
    [trackingData, trackingOptions] = trackParticles(fittingData,trackingOptions); %#ok<ASGLU>
    
    %Save trajectories
    save(posFit_list{iMovie},'trackingData','trackingOptions','-append');
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
    function exitFunc()
        % Remove TNTdata files
        if exist('list_filenames_TNTdata','var')
            warning off backtrace
            warning('User abort. Stopping TrackNTrace. Deleting settings files that have been saved already.');
            warning on backtrace
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
        clear global globalOptions candidateOptions fittingOptions trackingOptions movie imgCorrection parallelProcessingAvailable;
    end
end


