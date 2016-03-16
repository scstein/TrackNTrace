function [candidateData, fittingData, trackingData, previewOptions] = runPreview(movie,darkImage, candidateData, fittingData,trackingData, previewOptions, GUIreturns)
% Computes candidates, refines the position and tracks if requested.
%
% Note: previewOptions is a struct with the fields candidateOptions, fittingOptions, trackingOptions
%         We need this, as plugins might alter their options (e.g. add
%         information computed in the init function) when running. As
%         runnning the preview should have no side effects, we can not use
%         the global candidateOptions,fittingOptions,trackingOptions directly.
%         
%       GUIreturns is a struct output by the settings GUI, which gives us
%         information on which plugins settings changed so that we do not
%         need to compute steps (e.g. fit again) if parameters are unchanged
%

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
    TNTvisualizer(movie, candidateData, candidateOptions.outParamDescription, fittingData, fittingOptions.outParamDescription, trackingData, trackingOptions.outParamDescription, 5, true, globalOptions.firstFrameTesting);
else
    TNTvisualizer(movie, candidateData, candidateOptions.outParamDescription, fittingData, fittingOptions.outParamDescription, [], [], 5, true, globalOptions.firstFrameTesting);
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



