function [plugin] = plugin_crossCorrelation()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Cross correlation';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 1;

% The functions this plugin implements
mainFunc =  @findCandidates_crossCorrelation;

% Description of output parameters
outParamDescription = {'x';'y';'Corr'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Candidate detection based on matching a Gaussian PSF template to the image using normalized cross correlation.';

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
plugin.add_param('PSFsigma',...
    'float',...
    {1, 0, inf},...
    'Standard deviation of the PSF in pixels. sigma = FWHM/(2*sqrt(2*log(2))).');
plugin.add_param('CorrThreshold',...
    'float',...
    {0.4, 0, 1},...
    'Threshold between 0 and 1.');
plugin.add_param('distanceFactor',...
    'float',...
    {0.1, 2, inf},...
    'The search window size is 2*round(PSFsigma*distanceFactor)+1.');

end


%   -------------- User functions --------------

function candidatePos = findCandidates_crossCorrelation(img, options, currentFrame)
%Wrapper function for cross correlation candidate finding. Refer to
%matchSpot below or tooltips above to obtain information on input and
%output variables.
%
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidateOptions: Struct of input parameters provided by GUI.
%
% OUTPUT:
%     candidatePos - cell of Nx2 matrix of particle candidate positions
%     [column pixel, row pixel] without subpixel position. Middle of upper
%     left pixel would be [1,1].

sigma = options.PSFsigma;
CORR_THRESH = options.CorrThreshold;

% Candidate selection by normalized cross correlation
model = gaussian2d(round(10*sigma),sigma,true);
pattRows = size(model,1);
pattCols = size(model,2);

MIN_DIST = round(options.distanceFactor*sigma);
[ match_data, ~ ] = matchPattern( img, model, CORR_THRESH, MIN_DIST);

% convert to pixel position
match_data(:,1) = match_data(:,1)+ceil(pattRows/2);
match_data(:,2) = match_data(:,2)+ceil(pattCols/2);

% ATTENTION! Here we switch the coordinates
% candidatePos(:,1) is columns (x-coordinate)
% candidatePos(:,2) is rows    (y-coordinate)
candidatePos = match_data(:,[2,1]);

end


function f=gaussian2d(N,sigma, pixIntegrated)
%Usage   g_img = gaussian2d(N,sigma, pixIntegrated)
% N: grid size
% sigma: standard deviation
% pixIntegrated: Should the gaussian be integrated over the pixel area?
%
% The returned gaussian is normalized.

if nargin<3
    pixIntegrated = false;
end

[x, y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));

if(pixIntegrated)
    f = pi*sigma^2/2 * ( erfc(-(x-0+1/2)/(sqrt(2)*sigma)) - erfc(-(x-0-1/2)/(sqrt(2)*sigma))    )...
        .* ( erfc(-(y-0+1/2)/(sqrt(2)*sigma)) - erfc(-(y-0-1/2)/(sqrt(2)*sigma))   ) ;
else
    f= exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
end

f=f./sum(f(:));  %Normalize
end


function [ match_data, match_img ] = matchPattern( img, patt, CORR_THRESH, MIN_DIST)
% Detects a single pattern 'patt' withing the image 'img'. Only matches above the
% correlation threshold threshold 'CORR_THRESH' are taken as valid. Each
% returned match is the best candidate within a box window with edge
% length 2*MIN_DIST+1 centered around the candidate.
%
% USAGE: [ match_data, match_img ] = match_pattern(img, patt, CORR_THRESH, MIN_DIST)
%
% Input:
%    img         - The image to search in
%    patt        - The pattern to search for
%    CORR_THRESH - Minimum (absolute) normalized correlation value for a
%                  match to be valid. Valid values are between 0 and 1.
%    MIN_DIST    - Only the best candidates within a box window with edge
%                  length 2*MIN_DIST+1 centered around the candidate are taken.
%    FFTimg      - (optional) fft2 of img if known. see normxcorr2_preFFT for more info.
%    FFTpatt     - (optional) fft2 of pattern if known. see normxcorr2_preFFT for more info.
%
% Output:
%    match_data  - Mx3 matrix, where M is the number of valid matches. The
%                  columns correspond to (rowShift, colShift, correlation
%                  value). For rowShift = colShift = 0 the upper left pixel
%                  of the pattern lies ontop the upper left pixel of the
%                  image. The pattern centers can be calculated as
%                  (rowShift + pattRows/2, colShift + pattCols/2).
%    match_img   - Image of size (imgRows-pattRows)x(imgCols-pattCols)
%                  showing pattern shifts with pixels intensities
%                  corresponding to the correlation values.
%
% Note: Only patterns fully embedded in the image are detected.

%Note: CORR_THRESH depends on 'empty' space around the pattern.
% For larger zero padding and the presence of noise in the image
% this should be reduced, as the normalized correlation drops down.
% That said, zero padding seems to improve the performance of the detection
% (while sacrificing the usable border).

% Author: Simon Christoph Stein
% Date: 2014
% EMail: scstein@phys.uni-goettingen.de

imgRows = size(img,1);
imgCols = size(img,2);

pattRows = size(patt,1);
pattCols = size(patt,2);

MIN_DIST = min(max(1, MIN_DIST),min(imgRows,imgCols)); % Distance should be at least 1 and not greater than image size

% -- Pattern detection by cross-correlation --
% NOTE: Normalized cross correlation should to be used for feature matching,
% as otherwise the amplitude of the measured image will influence the
% result. For example: a spot much brighter than the pattern in the image
% will always lead to a local maximum.

% CC = xcorr2(img-mean(img(:)), patt);
CC = normxcorr2(patt,img); % normxcorr2 requires a smaller template than the image!
CC = CC(pattRows:imgRows, pattCols:imgCols); % cut out valid part of the cross-correlation (pattern fully embedded)

% Set all points below the thresold to the minimum
mask_overTOL = abs(CC)>CORR_THRESH;


% -- Detect local maxima within sliding window --
localmax_idx = localMax(CC, MIN_DIST);
mask_localmax = zeros(size(CC));
mask_localmax(localmax_idx) = 1;

mask = mask_localmax & mask_overTOL;


% Gather the positions of all local maxima
max_ind = find(mask);
[max_rPos, max_cPos] = ind2sub(size(CC), max_ind);
max_shifts = [max_rPos, max_cPos] -1; % (rowShift, colShift) of maxima

match_data = [max_shifts, CC(mask)]; % (rowShift, colShift, corrVal)


% image of matches with their correlation value
match_img = []; %not used here
% match_img = CC;
% match_img(~mask) = 0;
% match_img = padarray(match_img,[(pattRows-1)/2, (pattCols-1)/2]);
end



function max_indices = localMax(img, MIN_DIST)
% Returns indices of highest local maximum within a box window of edge
% length 2*MIN_DIST+1 around each pixel. Connected pixels with equal values
% are considered maxima as well.
neigh = ones(2*MIN_DIST+1); % Neighborhood matrix
img_dilated = imdilate(img, neigh);
max_indices = find(img == img_dilated);
end