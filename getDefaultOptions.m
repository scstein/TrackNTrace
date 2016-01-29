function [globalOptions, TNToptions] = getDefaultOptions()
% [globalOptions] = getDefaultGlobalOptions 
% 
% Set default options for TrackNTrace program and return them
% to main script.
% 
% OUTPUT:
%     globalOptions: struct containing general options (filename,
%     processing modes etc.)


% TNT internal options, these have nothing to do with the movie and are not
% saved along the other options
TNToptions.enableParallelProcessing = true;
TNToptions.closeMatlabpoolOnExit = false;
TNToptions.defaultCandidatePlugin = 'Cross correlation';
TNToptions.defaultFittingPlugin = 'TNT Fitter';
TNToptions.defaultTrackingPlugin = 'TNT NearestNeighbor';

% General Options of the GUI
globalOptions.filename_movies = ''; % This directory is chosen as startup when clicking to select movies
globalOptions.filename_dark_movie = ''; % This directory is chosen as startup when clicking to select a movie
globalOptions.previewMode = true; %enable/disable testMode by setting true/false.
globalOptions.firstFrame = 1; %integer, first movie frame to consider. 
globalOptions.lastFrame = inf; %integer, last movie frame to consider. Put inf to read to the end.
globalOptions.firstFrameTesting = 1; %integer, first movie frame to consider during testMode
globalOptions.lastFrameTesting = 50; %integer, last movie frame to consider during testMode. Put inf to read to the end.

% Options for photon conversion
globalOptions.usePhotonConversion = false; %boolean, enables photon conversion (absolutely necessary for MLE)
globalOptions.photonBias = 100; %integer, camera A/D count floor which is added to the image to remove negativ counts
globalOptions.photonSensitivity = 4; %double, electrons per image count, depends on pre-amp setting and readout mode. Check your camera manual and spec sheet!
globalOptions.photonGain = 100; %integer, EMCCD gain of camera

% Options for tracking particle positions
globalOptions.enableTracking = true;

end