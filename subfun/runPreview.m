% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
%     Copyright (C) 2016  Simon Christoph Stein, scstein@phys.uni-goettingen.de
%     Copyright (C) 2016  Jan Thiart, jthiart@phys.uni-goettingen.de
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [candidateData, refinementData, trackingData, previewOptions] = runPreview(movie,darkImage, candidateData, refinementData,trackingData, previewOptions, GUIreturns)
% Computes candidates, refines the position and tracks if requested.
%
% Note: previewOptions is a struct with the fields candidateOptions, refinementOptions, trackingOptions
%         We need this, as plugins might alter their options (e.g. add
%         information computed in the init function) when running. As
%         runnning the preview should have no side effects, we can not use
%         the global candidateOptions,refinementOptions,trackingOptions directly.
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
global refinementOptions
global trackingOptions

% Make copy of original globals which should be left untouched by the preview
cpy_globalOptions = globalOptions;
cpy_candidateOptions = candidateOptions;
cpy_refinementOptions = refinementOptions;
cpy_trackingOptions = trackingOptions;

% if no data is supplied, just run every step
if nargin < 3 || GUIreturns.globalOptionsChanged_ExcludingEnableTracking
    [candidateData, candidateOptions] = findCandidateParticles(movie, darkImage, globalOptions, candidateOptions);
    [refinementData, refinementOptions] = fitParticles(movie, darkImage, globalOptions, refinementOptions, candidateData);
    
    % track the particles
    if globalOptions.enableTracking
        [trackingData,trackingOptions] = trackParticles(refinementData,trackingOptions);
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
    
    % Perform refinement if required
    if GUIreturns.candidateOptionsChanged || GUIreturns.refinementOptionsChanged
        [refinementData, refinementOptions] = fitParticles(movie, darkImage, globalOptions, refinementOptions, candidateData);
    else % or use the options from the last preview run
        refinementOptions = previewOptions.refinementOptions;
    end
    
    
    % By invalidating tracking data we remember that the tracking options
    % were changed by the user, even if he set enableTracking = false in at the same time
    if GUIreturns.trackingOptionsChanged
        trackingData = [];
    end
    
    % Perform tracking if required
    if GUIreturns.candidateOptionsChanged || GUIreturns.refinementOptionsChanged || GUIreturns.trackingOptionsChanged || isempty(trackingData)
        if globalOptions.enableTracking
            % Track the particles
            [trackingData,trackingOptions] = trackParticles(refinementData,trackingOptions);
        end
    else % or use the options from the last preview run
        trackingOptions = previewOptions.trackingOptions;
    end    
end


%visualize all trajectories
if globalOptions.enableTracking
    TNTvisualizer(movie, candidateData, candidateOptions.outParamDescription, refinementData, refinementOptions.outParamDescription, trackingData, trackingOptions.outParamDescription, 5, true, globalOptions.firstFrameTesting);
else
    TNTvisualizer(movie, candidateData, candidateOptions.outParamDescription, refinementData, refinementOptions.outParamDescription, [], [], 5, true, globalOptions.firstFrameTesting);
end

% Store options from this preview
previewOptions.candidateOptions = candidateOptions;
previewOptions.refinementOptions = refinementOptions;
previewOptions.trackingOptions = trackingOptions;

% Restore original globals
globalOptions = cpy_globalOptions;
candidateOptions = cpy_candidateOptions;
refinementOptions = cpy_refinementOptions;
trackingOptions = cpy_trackingOptions;

end



