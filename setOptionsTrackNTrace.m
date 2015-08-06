function [candidateOptions,fittingOptions,trackingOptions] = setOptionsTrackNTrace()
% [candidateOptions,fittingOptions,trackingOptions] = setOptionsTrackNTrace 
% 
% Set options for TrackNTrace program and return them
% to main script. See relevant entries for details.
% 
% OUTPUT:
%     candidateOptions: struct containing options necessary for finding
%     particle localization candidates
%     
%     fittingOptions: struct containing options for Gaussian fitting
%     routine
%     
%     trackingOptions: struct containing options for tracking fitted
%     Positions


% Options for particle location candidates
%...general
candidateOptions.sigma = 1.3; %double, standard deviation in [pixel] of Gaussian pattern used for fitting. For non-moving particles in focus, the lower bound for sigma is lambda/(NA*4*sqrt(2*log(2)))
candidateOptions.fitForward = true; %boolean, fit movie forwards or backwards
candidateOptions.calculateCandidatesOnce = false; %boolean, if true, routine calculates candidates only once and tries refinement afterwards. DO NOT USE for diffusing particles
candidateOptions.averagingWindowSize = 10; %integer, average of X frames is calculated and used to find first candidates if calculateCandidatesOnce = true
candidateOptions.firstFrame = 1; %integer, first movie frame to consider
candidateOptions.lastFrame = 5000; %integer, last movie frame to consider

%...for normalized cross-correlation
candidateOptions.useCrossCorrelation = false; %boolean, either use cross correlation for finding particle candidates or intensity weighted filtering. The latter is recommended for diffusing particles
candidateOptions.corrThresh = 0.5; %double [0,1], threshold for cross correlation. higher means fewer candidates

%...for intensity filtering
candidateOptions.particleRadius = 3; %integer, expected radius of particles in [pixel] for intensity weighted filtering
candidateOptions.intensityThreshold = 5; %integer [0,100), only top X% of pixel intensity values are considered as candidates. Lower means fewer but better candiates
candidateOptions.intensityPtestVar = 0.02; %double [0,1], pixel is compared to Gaussian background at significance value of (1-X)*100%. Lower means fewer but better candidates
candidateOptions.backgroundCalculationInterval = 10; %integer, local background gets calculated every X frames to improve performance


% Options for fitting candidates
fittingOptions.fitSigma = true; %boolean, enables fitting of standard deviation. Literature recommends fixed sigma value for non-moving particles
fittingOptions.usePixelIntegratedFit = true; %boolean, enables pixel integrated Gaussian as fitting model (recommended). Disable for slighter better performance
fittingOptions.useMLErefine = true; %boolean, enables MLE fitting of particle positions. Performance intensive but more accurate in some cases


% Options for tracking particle positions
trackingOptions.method = 'utrack'; %string, choice of tracker (simpletracker, utrack, or track_cg)
trackingOptions.probDim = 2; %integer, problem dimension, 2 for xy-position. any tractable dimension (z, amplitude, pattern, background, ...) adds another dimension
trackingOptions.maxRadius = 6; %double, CAREFUL! maximum search radius of particle connections
trackingOptions.maxGap = 0; %intger, maximum allowed gap between two particle observations. Set to 0 if no gap closing desired
trackingOptions.minTrackLength = 2; %integer, minimum desired length of trajectories
trackingOptions.splitMovieParts = 5; %integer, splits positions array in X parts and tracks those individually. Might be necessary for utrack due to large memory consumption
trackingOptions.linkingMatrix = 'NearestNeighbor'; %string, choice of linking algorithm of simpletracker, 'Hungarian' or 'NearestNeighbor'
trackingOptions.verbose = false; %boolean, enables logging if tracker supports it

end