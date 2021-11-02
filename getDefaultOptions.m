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
%     Extended CT, 2019
%
function [globalOptions, TNToptions] = getDefaultOptions()
% [globalOptions] = getDefaultOptions 
% 
% Set default options for TrackNTrace program and return them
% to main script. Users may adapt these settings to their liking.
% 
% OUTPUT:
%     globalOptions: struct containing general options (filename,
%     processing modes etc.)
%     TNToptions: struct containing TNT internal options

% TNT internal options, these have nothing to do with the movie and are not
% saved along the other options
TNToptions.enableParallelProcessing = true;
TNToptions.closeMatlabpoolOnExit = false;
TNToptions.showResultLinks = true;

TNToptions.rememberSettingsForNextMovie = true; % If this is set to true, settings for the current movie are used as default for the next when pressing continue.

TNToptions.defaultCandidatePlugin = 'Cross correlation';
TNToptions.defaultRefinementPlugin = 'TNT Fitter';
TNToptions.defaultTrackingPlugin = 'TNT NearestNeighbor';
TNToptions.defaultPostprocPlugin = 'fit lifetime';

% General Options of the GUI
globalOptions.filename_dark_movie = ''; % This directory is chosen as startup when clicking to select a movie
globalOptions.previewMode = true; %enable/disable previewMode by setting true/false.
globalOptions.firstFrame = 1; %integer, first movie frame to consider. 
globalOptions.lastFrame = inf; %integer, last movie frame to consider. Put inf or 'end' to read to the end.
globalOptions.binFrame = 1; %integer, number of frames which are combined to one frame.
globalOptions.firstFrameTesting = 1; %integer, first movie frame to consider during testMode
globalOptions.lastFrameTesting = 50; %integer, last movie frame to consider during testMode. Put inf to read to the end.

% Options for photon conversion
globalOptions.usePhotonConversion = false; %boolean, enables photon conversion (absolutely necessary for MLE)
globalOptions.photonBias = 100; %integer, camera A/D count floor which is added to the image to remove negativ counts
globalOptions.photonSensitivity = 5; %double, electrons per image count, depends on pre-amp setting and readout mode. Check your camera manual and spec sheet!
globalOptions.photonGain = 100; %integer, EMCCD gain of camera

% Options for timegate
globalOptions.useTimegate = false; %boolean, enables timegate
globalOptions.tgStart = 0; 
globalOptions.tgEnd = inf; 

% Options to enable plugins
globalOptions.enableCandidate  = true;
globalOptions.enableRefinement = true;
globalOptions.enableTracking   = false;
globalOptions.enablePostproc   = false;

end