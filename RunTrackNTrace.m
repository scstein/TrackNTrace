%% Add required folders and subfolders to path
addpath(genpath('external'));
addpath(genpath('helper'));
addpath(genpath('subfun'));
addpath(genpath('testing'));

%% Load and adjust the default settings for this batch
GUIinputs.titleText = 'Please select a list of movies to process.';
GUIinputs.fileText  = 'Default settings for this batch';
GUIinputs.singleFileMode = false;
[generalOptions_def, candidateOptions_def,fittingOptions_def,trackingOptions_def] = setDefaultOptions();
[generalOptions_def, candidateOptions_def,fittingOptions_def,trackingOptions_def, GUIreturns] = settingsGUI(generalOptions_def, candidateOptions_def,fittingOptions_def,trackingOptions_def, GUIinputs);

%% Adjust options for each movie and test settings if desired
GUIinputs.singleFileMode = true; % No editing of movie list possible

% [movie_list,dark_stack] = getMovieFilenames(generalOptions.filename_movies, generalOptions.filename_dark_movie);
movie_list = generalOptions_def.filename_movies;
% Calculate default dark image if given
dark_img_def = [];
if(~isempty(generalOptions_def.filename_dark_movie))
    dark_img_def = CalculateDark(read_tiff(generalOptions_def.filename_dark_movie));
end

% Get timestamp for output files
time = clock;
timestamp = sprintf('%i-m%i-d%i-%ih%i',time(1),time(2),time(3),time(4),time(5));

posFit_list = cell(0);
for i=1:numel(movie_list)
    filename_movie = movie_list{i};
    [path,filename,~] = fileparts(filename_movie);
    filename_fitData = [path,filesep,filename,'_',timestamp,'_TNT.mat'];
    
    
    %     Skip nonsensical input
    [T,movie] = evalc(['read_tiff(''',filename_movie,''',false,[1,2])']); % Read 2 frames. note: evalc suppresses output
    if size(movie,3)<=1
       continue; 
    end
    
    % Set options to default for this batch
    generalOptions = generalOptions_def;
    candidateOptions = candidateOptions_def;
    fittingOptions = fittingOptions_def;
    trackingOptions = trackingOptions_def;
    dark_img = dark_img_def;
    
    % Does the user want to adjust the settings per movie?
    if not(GUIreturns.useSettingsForAll)
        % Show filename in GUI
        GUIinputs.fileText = filename_movie;
        
        GUIinputs.titleText = 'Adjust movie specific options.';
        [generalOptions, candidateOptions,fittingOptions,trackingOptions] = settingsGUI(generalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs);
        
        
        % Check if different dark movie was given
        if(~strcmp(generalOptions_def.filename_dark_movie, generalOptions.filename_dark_movie))
            if(~isempty(generalOptions.filename_dark_movie))
                dark_img   = CalculateDark(read_tiff(generalOptions.filename_dark_movie));
            end
        end
        
        % If test mode is enabled, analyze first X frames and show GUI
        if generalOptions.previewMode
            run_again = true;
            first_run = true;
            firstFrameTesting = 0;
            lastFrameTesting  = 0;
            filename_dark_movie = generalOptions.filename_dark_movie;
            while run_again
                if not(first_run); [generalOptions, candidateOptions,fittingOptions,trackingOptions] = settingsGUI(generalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs); end;
                if not(generalOptions.previewMode); break; end; % If test mode was disabled by user in the settingsGUI
                % Check if requested frame interval has changed -> re-read movie if neccessary
                if (firstFrameTesting ~= generalOptions.firstFrameTesting) || (lastFrameTesting ~= generalOptions.lastFrameTesting)
                    firstFrameTesting = generalOptions.firstFrameTesting;
                    lastFrameTesting  = generalOptions.lastFrameTesting;
                    movie = read_tiff(filename_movie, false, [firstFrameTesting, lastFrameTesting]);
                end
                % Check if different dark movie was given
                if(~strcmp(filename_dark_movie, generalOptions.filename_dark_movie))
                    if(~isempty(generalOptions.filename_dark_movie))
                        dark_img   = CalculateDark(read_tiff(generalOptions.filename_dark_movie));
                    end
                    filename_dark_movie = generalOptions.filename_dark_movie;
                end
                [run_again] = testTrackerSettings(movie,dark_img,candidateOptions,fittingOptions,trackingOptions);
                first_run = false;
            end
        end
        
    end
    % Save options
    save(filename_fitData,'filename_movie','generalOptions','candidateOptions','fittingOptions','trackingOptions','dark_img');
    posFit_list = [posFit_list;{filename_fitData}]; %#ok<AGROW>
end

%% Compute positions
clearvars -except posFit_list
for i=1:numel(posFit_list)
    filename_fitData = posFit_list{i};
    load(filename_fitData,'-mat');
    
    % Read movie  
    movie = read_tiff(filename_movie, false, [generalOptions.firstFrame,generalOptions.lastFrame]);
    
    % Compute the positions
    fprintf('\n######\nLocating particles in movie %s.\n',filename_movie);
    fitData = locateParticles(movie, dark_img, candidateOptions, fittingOptions); %#ok<NASGU>
    
    % Save positions
    save(filename_fitData,'fitData','-append');
end

clear fitData movie
%% TODO
% correct for half pixel shift?
% convert amplitude to photons
% keep amplitude even if tracker cant manage
% replace Ptest by something better? --> amp>mean+k*std, k by user
% automatically determine if memory large enough

%% Compute trajectories
for i=1:numel(posFit_list)
    load(posFit_list{i},'trackingOptions','fitData');
    
    % If no tracking is desired for this movie, continue
    if (~trackingOptions.enableTracking)
        continue
    end
    
    % Compute trajectories
    fprintf('\n######\nTracking particles in movie %s.\n',filename_movie);
    trajectoryData = trackParticles(fitData,trackingOptions);
    
    %Save trajectories
    save(posFit_list{i},'trajectoryData','-append');
end

