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
[globalOptions_def] = getDefaultGlobalOptions(); 

%% Check if parallel processing is available
global parallelProcessingAvailable
parallelProcessingAvailable = false;
closeMatlabpoolOnExit = globalOptions_def.closeMatlabpoolOnExit ;

if globalOptions_def.enableParallelProcessing    
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
GUIinputs.titleText = 'Please select a list of movies to process.';
GUIinputs.fileText  = 'Default settings for this batch';
GUIinputs.singleFileMode = false; % false -> movie list can be edited

[globalOptions_def, candidateOptions_def,fittingOptions_def,trackingOptions_def, GUIreturns] = settingsGUI(globalOptions_def, [],[],[], GUIinputs);
if GUIreturns.userExit;
    exitFunc()
    return; 
end;

%% Adjust options for each movie and test settings if desired
GUIinputs.singleFileMode = true; % No editing of movie list possible

movie_list = globalOptions_def.filename_movies;

% Calculate default dark image if given
dark_img_def = [];
if(~isempty(globalOptions_def.filename_dark_movie))
    dark_img_def = CalculateDark(read_tiff(globalOptions_def.filename_dark_movie));
end

% Create timestamp for output files
time = clock;
timestamp = sprintf('%i-m%02i-d%02i-%02ih%i',time(1),time(2),time(3),time(4),time(5));

% Iterate through all movies in the list
posFit_list = cell(0);
for i=1:numel(movie_list)
    filename_movie = movie_list{i};
    [path,filename,~] = fileparts(filename_movie);
    filename_fittingData = [path,filesep,filename,'_',timestamp,'_TNT.mat'];
    
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
    
    % Set options to default for this movie
    globalOptions = globalOptions_def;
    candidateOptions = candidateOptions_def;
    fittingOptions = fittingOptions_def;
    trackingOptions = trackingOptions_def;
    dark_img = dark_img_def;
    
    % Does the user want to adjust the settings per movie?
    if not(  GUIreturns.useSettingsForAll )
        % Show filename in GUI
        GUIinputs.fileText = filename_movie;
        
        GUIinputs.titleText = 'Adjust movie specific options.';
        [globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIreturns] = settingsGUI(globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs);
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
        if GUIreturns.useSettingsForAll; globalOptions.previewMode = false; end;
        
        % If test mode is enabled, analyze first X frames and show GUI
        if globalOptions.previewMode
            run_again = true;
            first_run = true;
            filename_dark_movie = globalOptions.filename_dark_movie;
            while run_again
                if not(first_run)
                    [globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIreturns] = settingsGUI(globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs);
                    if GUIreturns.userExit;
                        exitFunc();
                        return;
                    end;
                    if GUIreturns.useSettingsForAll; globalOptions.previewMode = false; end; %dont go through other movies anymore
                end
                
                if not(globalOptions.previewMode); break; end; % If test mode was disabled by user in the settingsGUI
                % Check if requested frame interval has changed -> re-read movie if neccessary
                if first_run || GUIreturns.testWindowChanged
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
                
                % IF: this is the first run, the preview window changed or the fitting/candidate options changed locate and
                % track particles and save fittingData. ELSE: reuse fittingData acquired in the last run without re-fitting
                if first_run || GUIreturns.globalOptionsChanged || GUIreturns.fittingOptionsChanged || GUIreturns.candidateOptionsChanged
                    [run_again, candidateData_test,fittingData_test] = runPreview(movie,dark_img);
                else
                    [run_again] = runPreview(movie,dark_img, candidateData_test, fittingData_test);
                end
                first_run = false;
            end
        end
        
    end %not( GUIreturns.useSettingsForAll )
    
    if GUIreturns.useSettingsForAll
        globalOptions_def = globalOptions;
        candidateOptions_def = candidateOptions;
        fittingOptions_def = fittingOptions;
        trackingOptions_def = trackingOptions;
        dark_img_def = dark_img;
    end
    
    % Save options
    globalOptions.filename_movies = {filename_movie}; % Save only name of this file in its settings (important when loading options)
    save(filename_fittingData,'filename_movie','globalOptions','candidateOptions','fittingOptions','trackingOptions','dark_img');
    posFit_list = [posFit_list;{filename_fittingData}]; %#ok<AGROW>
end
clearvars -except posFit_list

%% Candidate detection and fitting for every movie
for i=1:numel(posFit_list)
    filename_fittingData = posFit_list{i};
    
    load(filename_fittingData,'-mat');
    
    % Read movie
    movie = read_tiff(filename_movie, false, [globalOptions.firstFrame,globalOptions.lastFrame]);
    % Compute the positions
    fprintf('######\nTNT: Locating particles in movie %s.\n',filename_movie);
    [candidateData, candidateOptions] = findCandidateParticles(movie, dark_img, globalOptions, candidateOptions);
    [fittingData, fittingOptions] = fitParticles(movie, dark_img, globalOptions, fittingOptions, candidateData);
    
    % Save positions and movieSize, update globalOptions.lastFrame
    globalOptions.lastFrame = globalOptions.firstFrame + size(movie,3)-1; % lastFrame could have been set to 'inf', now we synchronize with the correct number
    movieSize = size(movie); %#ok<NASGU> % Save size of movie (nice to have)
    save(filename_fittingData,'candidateData','fittingData','globalOptions','candidateOptions','fittingOptions','movieSize','-append');
end
clearvars -except posFit_list

%% Compute trajectories for every movie
for i=1:numel(posFit_list)
    load(posFit_list{i},'globalOptions','trackingOptions','fittingData','filename_movie');
    
    % If no tracking is desired for this movie, continue
    if (~globalOptions.enableTracking)
        continue
    end
    
    % Compute trajectories
    fprintf('######\nTNT: Tracking particles in movie %s.\n',filename_movie);
    [trackingData, trackingOptions] = trackParticles(fittingData,trackingOptions); %#ok<ASGLU>
    
    %Save trajectories
    save(posFit_list{i},'trackingData','trackingOptions','-append');
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
    warning off backtrace
    warning(sprintf('User abort. Stopping TrackNTrace.\nDelete unwanted settings files that might have been saved already.'));
    warning on backtrace
    if parallelProcessingAvailable && closeMatlabpoolOnExit
        matlabpool('close');
    end
    clearGlobals();
end

% Clear all global variables
function clearGlobals()
    clear global globalOptions candidateOptions fittingOptions trackingOptions movie imgCorrection parallelProcessingAvailable;
end
end


