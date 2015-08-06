function [] = testCandidateSearch(movie, img_dark, candidateOptions,frame)
% testCandidateSearch(movie, img_dark, candidateOptions,frame)
% Helper function to test candidate search parameters, does not return
% anything. Refer to RunTrackNTrace or the respective candidate search
% functions for help


if ~isempty(img_dark)
    img = double(movie(:,:,frame))+img_dark;
else
    img = double(movie(:,:,frame));
end

[~] = findSpotCandidates(img, candidateOptions.sigma, candidateOptions.corrThresh,true);
[~] = findSpotCandidates_MOSAIC(img,candidateOptions.particleRadius,candidateOptions.intensityThreshold,candidateOptions.intensityPtestVar,0,true,true);

% as usual, remove global variables set by findSpotCandidates_MOSAIC
clear global img_bck_mean img_bck_std

