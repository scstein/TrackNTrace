function [run_again] = testTrackerSettings(movie_filename,dark_stack,candidateOptions,fittingOptions,trackingOptions)

%read first 50 frames, hopefully this is enough for testing. If a movie is
%smaller than that, we're screwed as read_tiff_DEMO doesn't check for movie
%length to be fast enough
movie = read_tiff_DEMO(movie_filename, false, [1,50]);
%don't correct for dark image counts if size doesn'T match
if ~isempty(dark_stack)
    if sum([size(movie,1),size(movie,2)]==[size(dark_stack,1),size(dark_stack,2)])~=2 || size(movie,3)<=1
        dark_stack = [];
    end
end

%slicing up a small movie diesn't make sense, so don't
trackingOptions.splitMovieParts = 1;

%get particle positions, track them
fitData = locateParticles(movie, dark_stack, candidateOptions, fittingOptions);
trajectoryData = trackParticles(fitData,trackingOptions);

%visualize all trajectories
[hGUI, selected_option] = visualizeTracksGUI(movie,trajectoryData,5,[],[],[],true);

run_again = strcmp(selected_option,'Yes');

