function [run_again, candidateData, fitData] = runPreview(movie,darkImage, candidateData, fitData)
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

% We need the global access here, as plugins might want to access/change
% the values in their options or options of other plugins.
% This is somewhat shitty and can only be changed by recoding Track'N'Trace
% in a class based way;
global globalOptions
global candidateOptions
global fittingOptions
global trackingOptions

% Make copy of original options which should be left untouched by the preview
cpy_globalOptions = globalOptions;
cpy_candidateOptions = candidateOptions;
cpy_fittingOptions = fittingOptions;
cpy_trackingOptions = trackingOptions;

% if not supplied by the user, get the particle positions
if nargin < 7 || isempty(vertcat(fitData{:}))
    [candidateData, candidateOptions] = findCandidateParticles(movie, darkImage, globalOptions, candidateOptions);
    [fitData, fittingOptions] = fitParticles(movie, darkImage, globalOptions, fittingOptions, candidateData);
end

% track the particles
if globalOptions.enableTracking
    [trajectoryData,trackingOptions] = trackParticles(fitData,trackingOptions);
end

%visualize all trajectories
if globalOptions.enableTracking
    [hGUI, run_again] = visualizeTracksGUI(movie,trajectoryData,5,[],[],[],[],true);
else
    [hGUI, run_again] = visualizeFitDataGUI(movie,fitData, fittingOptions.outParamDescription, 5,[],true);
end

% Restore original options
globalOptions = cpy_globalOptions;
candidateOptions = cpy_candidateOptions;
fittingOptions = cpy_fittingOptions;
trackingOptions = cpy_trackingOptions;
    
end



