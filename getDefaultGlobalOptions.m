function [globalOptions] = getDefaultGlobalOptions()
% [candidateOptions,fittingOptions,trackingOptions] = setOptionsTrackNTrace 
% 
% Set default options for TrackNTrace program and return them
% to main script. See relevant entries for details.
% 
% OUTPUT:
%     globalOptions: struct containing general options (filename,
%     processing modes etc.)

% General Options
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