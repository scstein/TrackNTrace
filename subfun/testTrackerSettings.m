function [run_again, candidateData, fitData] = testTrackerSettings(movie,darkImage,globalOptions,candidateOptions,fittingOptions,trackingOptions, candidateData, fitData)
% 
% %read first 50 frames, hopefully this is enough for testing. If a movie is
% %smaller than that, we're screwed as read_tiff_DEMO doesn't check for movie
% %length to be fast enough
% movie = read_tiff_DEMO(movie_filename, false, [1,50]);
% %don't correct for dark image counts if size doesn'T match
% if ~isempty(dark_stack)
%     if sum([size(movie,1),size(movie,2)]==[size(dark_stack,1),size(dark_stack,2)])~=2 || size(movie,3)<=1
%         dark_stack = [];
%     end
% end

% if not supplied by the user, get the particle positions
if nargin < 7 || isempty(vertcat(fitData{:}))
    [candidateData, ~] = findCandidateParticles(movie, darkImage, globalOptions, candidateOptions);
    [fitData, ~] = fitParticles(movie, darkImage, globalOptions, fittingOptions, candidateData);
end
% track the particles
if globalOptions.enableTracking
    trajectoryData = trackParticles(fitData,trackingOptions);
end

%visualize all trajectories
if globalOptions.enableTracking
    [hGUI, run_again] = visualizeTracksGUI(movie,trajectoryData,5,[],[],[],[],true);
else
    [hGUI, run_again] = visualizeFitDataGUI(movie,fitData,5,[],true);
end


