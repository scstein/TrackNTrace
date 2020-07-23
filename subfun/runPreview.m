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
% 	CT, 2020:
% 	- support for postprossing plugins
%	- struct based call of TNTvisualizer
%
function [candidateData, refinementData, trackingData, postprocData, previewOptions] = runPreview(movie,darkImage,metadata, candidateData, refinementData, trackingData, postprocData, previewOptions, GUIreturns)
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
global importOptions
global candidateOptions
global refinementOptions
global trackingOptions
global postprocOptions

% Make copy of original globals which should be left untouched by the preview
% Note: the importOptions should not be changed by the plugins.
cpy_globalOptions = globalOptions;
cpy_candidateOptions = candidateOptions;
cpy_refinementOptions = refinementOptions;
cpy_trackingOptions = trackingOptions;
cpy_postprocOptions = postprocOptions;


% if no data is supplied or globalOptions changed run everything
if nargin < 4 || GUIreturns.globalOptionsChanged_ExcludingEnable || GUIreturns.importOptionsChanged
    % Initilise output data
    [candidateData, refinementData, trackingData, postprocData]=deal([]);
else
    % if the options changed delete the old data and all dependent data.
    % Otherwise restore the old options.
    if GUIreturns.candidateOptionsChanged || ~globalOptions.enableCandidate
        candidateData = [];
    else
        candidateOptions = previewOptions.candidateOptions;
    end
    if GUIreturns.refinementOptionsChanged || ~globalOptions.enableRefinement || isempty(candidateData)
        refinementData = [];
    else
        refinementOptions = previewOptions.refinementOptions;
    end
    if GUIreturns.trackingOptionsChanged || ~globalOptions.enableTracking || isempty(refinementData)
        trackingData = [];
    else
        trackingOptions = previewOptions.trackingOptions;
    end
    if GUIreturns.postprocOptionsChanged || ~globalOptions.enablePostproc || isempty(trackingData)
        postprocData = [];
    else
        postprocOptions = previewOptions.postprocOptions;
    end
end

% Execute step if data is requested and not existing
if globalOptions.enableCandidate && isempty(candidateData)
    [candidateData, candidateOptions] = findCandidateParticles(movie{1}, darkImage, globalOptions, candidateOptions);
end
if globalOptions.enableRefinement && isempty(refinementData) && ~isempty(candidateData)
    [refinementData, refinementOptions] = fitParticles(movie{1}, darkImage, globalOptions, refinementOptions, candidateData);
end
if globalOptions.enableTracking && isempty(trackingData) && ~isempty(refinementData)
    [trackingData,trackingOptions] = trackParticles(refinementData,trackingOptions);
end
if globalOptions.enablePostproc && isempty(postprocData) && ~isempty(trackingData)
    [postprocData,postprocOptions] = postprocTracks(trackingData, globalOptions,importOptions, postprocOptions);
end

fprintf('\n');

% Build TNTdata struct
if isstruct(metadata) && isfield(metadata,'filename')
    [~,titlename] = fileparts(metadata.filename);
    titlename = ['Preview of ', titlename];
else
    titlename = 'Preview';
end

TNTdata = struct(...
    'title',titlename,...
    'is_blocking',true,...
    'firstFrame_lastFrame',[globalOptions.firstFrameTesting, globalOptions.lastFrameTesting],...
    'FPS',5,...
    'metadata',metadata,...
    'globalOptions',globalOptions...
    );

% tntvisoptions = cell(1,8);
if globalOptions.enableCandidate
%     tntvisoptions(1:2) = {candidateData, candidateOptions.outParamDescription};
    TNTdata.candidateData = candidateData;
    TNTdata.candidateOptions = candidateOptions;
end
if globalOptions.enableRefinement
%     tntvisoptions(3:4) = {refinementData, refinementOptions.outParamDescription};
    TNTdata.refinementData = refinementData;
    TNTdata.refinementOptions = refinementOptions;
end
if globalOptions.enableTracking
%     tntvisoptions(5:6) = {trackingData, trackingOptions.outParamDescription};
    TNTdata.trackingData = trackingData;
    TNTdata.trackingOptions = trackingOptions;
end
if globalOptions.enablePostproc
%     tntvisoptions(7:8) = {postprocData, postprocOptions.outParamDescription};
    TNTdata.postprocData = postprocData;
    TNTdata.postprocOptions = postprocOptions;
end
% tntvisoptions = [tntvisoptions,{5, true, globalOptions.firstFrameTesting}];
% TNTvisualizer(movie, tntvisoptions{:});
TNTvisualizer(movie, TNTdata);

% Store options from this preview
previewOptions.candidateOptions = candidateOptions;
previewOptions.refinementOptions = refinementOptions;
previewOptions.trackingOptions = trackingOptions;
previewOptions.postprocOptions = postprocOptions;

% Restore original globals
globalOptions = cpy_globalOptions;
candidateOptions = cpy_candidateOptions;
refinementOptions = cpy_refinementOptions;
trackingOptions = cpy_trackingOptions;
postprocOptions = cpy_postprocOptions;

end