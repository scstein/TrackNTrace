function [] = TrackNTraceBatch(movie_list,dark_img,generalOptions,candidateOptions,fittingOptions,trackingOptions,filterOptions)
% TrackNTraceBatch(movie_list,dark_img,generalOptions,candidateOptions,fittingOptions,trackingOptions,filterOptions)
% The function complements RunTrackNTrace in that it can be used to handle
% a large batch of movies by directly giving all relevant settings and
% movie names. See RunTrackNTrace for details.
% 
% INPUT:
%     movie_list: cell array of full movie path strings
%     
%     dark_img: dark stack correction image
%     
%     <default>Options: settings structs
%     
%     filterOptions: struct for filtering
%         .fun: filter function name (must be in search path) .args: cell
%         array of arguments
    
addpath(genpath('external'));
addpath(genpath('helper'));
addpath(genpath('subfun'));
addpath(genpath('testing'));

% Get timestamp for output files
time = clock;
timestamp = sprintf('%i-m%02i-d%02i-%ih%i',time(1),time(2),time(3),time(4),time(5));

posFit_list = cell(0);
for i=1:numel(movie_list)
    filename_movie = movie_list{i};
    [path,filename,~] = fileparts(filename_movie);
    filename_fitData = [path,filesep,filename,'_',timestamp,'_TNT.mat'];
    
    % Skip nonsensical input
    [~,movie] = evalc(['read_tiff(''',filename_movie,''',false,[1,2])']); % Read 2 frames. note: evalc suppresses output
    if size(movie,3)<=1
       continue; 
    end
    
    % Save options
    save(filename_fitData,'filename_movie','generalOptions','candidateOptions','fittingOptions','trackingOptions','dark_img','filterOptions');
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
    fprintf('######\nLocating particles in movie %s.\n',filename_movie);
    fitData = locateParticles(movie, dark_img, candidateOptions, fittingOptions);
    
    if ~isempty(filterOptions)
        filter_eval_string = [filterOptions.fun,'(fitData,'];
        for iArgs=1:numel(filterOptions.args)-1
            filter_eval_string = [filter_eval_string,filterOptions.args{iArgs},',']; %#ok<AGROW>
        end
        filter_eval_string = [filter_eval_string,num2str(filterOptions.args{end}),')']; %#ok<AGROW>
        
        fitData = eval(filter_eval_string);
    end
    
    % Save positions
    save(filename_fitData,'fitData','-append');
end

clear fitData movie


%% Compute trajectories
for i=1:numel(posFit_list)
    load(posFit_list{i},'trackingOptions','fitData','filename_movie');
    
    % If no tracking is desired for this movie, continue
    if (~trackingOptions.enableTracking)
        continue
    end
    
    % Compute trajectories
    fprintf('######\nTracking particles in movie %s.\n',filename_movie);
    trajectoryData = trackParticles(fitData,trackingOptions); %#ok<NASGU>
    
    %Save trajectories
    save(posFit_list{i},'trajectoryData','-append');
end

