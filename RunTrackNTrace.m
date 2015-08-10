%% Get filenames
filename.movies = 'U:\Jan Thiart\dicty\2015-07-17\';
filename.dark_movies = [];
[movie_list,dark_stack] = readMovies_TrackNTrace(filename);

%% Configure options
testMode = [true,1]; %enable/disable testMode by setting demoMode(1) true/false. A small portion of the movie at position demoMode(2) of movie_list is run through by the tracking algorithm to
[candidateOptions,fittingOptions,trackingOptions] = setOptionsTrackNTrace();

%% Test settings if desired
if testMode(1)
    run_again = true;
    movie = read_tiff(movie_list{testMode(2)}, false, 50);
    if ~isempty(dark_stack)
        if sum([size(movie,1),size(movie,2)]==[size(dark_stack,1),size(dark_stack,2)])~=2 || size(movie,3)<=1
            dark_stack = [];
        end
    end
    while run_again
        [run_again] = testTrackerSettings(movie,dark_stack,candidateOptions,fittingOptions,trackingOptions);
        rehash;
        [candidateOptions,fittingOptions,trackingOptions] = setOptionsTrackNTrace();
        
    end
end

%% Compute positions
posFit_list = cell(0);
for i=1:numel(movie_list)
    movie_name = movie_list{i};
    
    % Read movie
    movie = read_tiff(movie_list{i}, false);
    movie = movie(:,:,candidateOptions.firstFrame:min(candidateOptions.lastFrame,size(movie,3)));
    
    % Skip nonsensical input
    if ~isempty(dark_stack)
        if sum([size(movie,1),size(movie,2)]==[size(dark_stack,1),size(dark_stack,2)])~=2 || size(movie,3)<=1
            continue
        end
    end
    
    % Test candidate search
%     testCandidateSearch(movie, dark_stack(:,:,1), candidateOptions,1000); %last value: frame where seach is tested
    
    % Compute the positions
    fitData = locateParticles(movie, dark_stack, candidateOptions, fittingOptions);
    
    % Save positions
    filename_fitData = [movie_name(1:end-4),'_TNT'];
    save(filename_fitData,'candidateOptions','fittingOptions','fitData');
    posFit_list = [posFit_list;{filename_fitData}]; %#ok<AGROW>
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
    % Compute trajectories
    trajectoryData = trackParticles(posFit_list{i},trackingOptions);
    
    %Save trajectories
    save(posFit_list{i},'trajectoryData','trackingOptions','-append');
end

