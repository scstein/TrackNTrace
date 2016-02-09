function [candidateData, fittingData, trackingData, previewOptions] = runPreview(movie,darkImage, candidateData, fittingData,trackingData, previewOptions, GUIreturns)
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
% This is somewhat shitty and can probably only be changed by recoding
% Track'N'Trace in a class based way;
global globalOptions
global candidateOptions
global fittingOptions
global trackingOptions

% Make copy of original globals which should be left untouched by the preview
cpy_globalOptions = globalOptions;
cpy_candidateOptions = candidateOptions;
cpy_fittingOptions = fittingOptions;
cpy_trackingOptions = trackingOptions;

% if no data is supplied, just run every step
if nargin < 3 || GUIreturns.globalOptionsChanged_ExcludingEnableTracking
    [candidateData, candidateOptions] = findCandidateParticles(movie, darkImage, globalOptions, candidateOptions);
    [fittingData, fittingOptions] = fitParticles(movie, darkImage, globalOptions, fittingOptions, candidateData);
    
    % track the particles
    if globalOptions.enableTracking
        [trackingData,trackingOptions] = trackParticles(fittingData,trackingOptions);
    else
        trackingData = [];
    end
else % Reuse data where possible
    % Perform candidate search if required
    if GUIreturns.candidateOptionsChanged
        [candidateData, candidateOptions] = findCandidateParticles(movie, darkImage, globalOptions, candidateOptions);
    else % or use the options from the last preview run
        candidateOptions = previewOptions.candidateOptions;
    end
    
    % Perform fitting if required
    if GUIreturns.candidateOptionsChanged || GUIreturns.fittingOptionsChanged
        [fittingData, fittingOptions] = fitParticles(movie, darkImage, globalOptions, fittingOptions, candidateData);
    else % or use the options from the last preview run
        fittingOptions = previewOptions.fittingOptions;
    end
    
    
    % By invalidating tracking data we remember that the tracking options
    % were changed by the user, even if he set enableTracking = false in at the same time
    if GUIreturns.trackingOptionsChanged
        trackingData = [];
    end
    
    % Perform tracking if required
    if GUIreturns.candidateOptionsChanged || GUIreturns.fittingOptionsChanged || GUIreturns.trackingOptionsChanged || isempty(trackingData)
        if globalOptions.enableTracking
            % Track the particles
            [trackingData,trackingOptions] = trackParticles(fittingData,trackingOptions);
        end
    else % or use the options from the last preview run
        trackingOptions = previewOptions.trackingOptions;
    end    
end


%visualize all trajectories
if globalOptions.enableTracking
%     [hGUI] = visualizeTracksGUI(movie,trackingData,5,[],[],[],[],true);
    [hGUI] = TNTvisualizer(movie, candidateData, candidateOptions.outParamDescription, fittingData, fittingOptions.outParamDescription, trackingData, trackingOptions.outParamDescription, 5, true);
else
%     [hGUI] = visualizeFitDataGUI(movie,fittingData, fittingOptions.outParamDescription, 5,[],true);
    [hGUI] = TNTvisualizer(movie, candidateData, candidateOptions.outParamDescription, fittingData, fittingOptions.outParamDescription, [], [], 5, true);
end

% Store options from this preview
previewOptions.candidateOptions = candidateOptions;
previewOptions.fittingOptions = fittingOptions;
previewOptions.trackingOptions = trackingOptions;

% Restore original globals
globalOptions = cpy_globalOptions;
candidateOptions = cpy_candidateOptions;
fittingOptions = cpy_fittingOptions;
trackingOptions = cpy_trackingOptions;

end



