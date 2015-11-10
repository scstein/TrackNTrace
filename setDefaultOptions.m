function [generalOptions] = setDefaultOptions()
% [candidateOptions,fittingOptions,trackingOptions] = setOptionsTrackNTrace 
% 
% Set default options for TrackNTrace program and return them
% to main script. See relevant entries for details.
% 
% OUTPUT:
%     generalOptions: struct containing general options (filename,
%     processing modes etc.)

% General Options
generalOptions.filename_movies = ''; % This directory is chosen as startup when clicking to select movies
generalOptions.filename_dark_movie = ''; % This directory is chosen as startup when clicking to select a movie
generalOptions.previewMode = true; %enable/disable testMode by setting true/false.
generalOptions.firstFrame = 1; %integer, first movie frame to consider. 
generalOptions.lastFrame = inf; %integer, last movie frame to consider. Put inf to read to the end.
generalOptions.firstFrameTesting = 1; %integer, first movie frame to consider during testMode
generalOptions.lastFrameTesting = 50; %integer, last movie frame to consider during testMode. Put inf to read to the end.

% Options for candidate detection behaviour
generalOptions.fitForward = true; %boolean, fit movie forwards or backwards
generalOptions.calculateCandidatesOnce = false; %boolean, if true, routine calculates candidates only once and tries refinement afterwards. DO NOT USE for diffusing particles
generalOptions.averagingWindowSize = 20; %integer, average of X frames is calculated and used to find first candidates if calculateCandidatesOnce = true

% Options for photon conversion
generalOptions.usePhotonConversion = false; %boolean, enables photon conversion (absolutely necessary for MLE)
generalOptions.photonBias = 100; %integer, camera A/D count floor which is added to the image to remove negativ counts
generalOptions.photonSensitivity = 4; %double, electrons per image count, depends on pre-amp setting and readout mode. Check your camera manual and spec sheet!
generalOptions.photonGain = 100; %integer, EMCCD gain of camera

% Options for tracking particle positions
generalOptions.enableTracking = true;



end