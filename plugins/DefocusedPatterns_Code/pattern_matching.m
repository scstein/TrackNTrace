function [match_data, im_detect] = pattern_matching(img, patterns, CORR_THRESH, MIN_DIST, FFTimg, FFTpatterns, compute_im_detect, use_parallel_processing)
% Detects a set of patterns within the given image. Only matches above the
% correlation threshold threshold 'CORR_THRESH' are taken as valid. Each
% returned match is the best candidate within a box window with edge
% length 2*MIN_DIST+1 centered around the candidate.
%
% USAGE: [match_data, im_theo] = pattern_matching(img, patterns, CORR_THRESH, MIN_DIST)
%
% Input:
%    img         - The image to search in.
%    patterns    - (Length X Width X #Patterns) matrix holding the patterns to search for.
%    CORR_THRESH - Minimum (absolute) normalized correlation value for a
%                  match to be valid. Valid values are between 0 and 1.
%    MIN_DIST    - Only the best candidates within a box window with edge
%                  length 2*MIN_DIST+1 centered around the candidate are taken.
%    FFTimg      - (optional) FFT of img if known. see normxcorr2_preFFT for more info.
%    FFTpatterns - (optional) Rotated FFTs of patterns if known. see normxcorr2_preFFT for more info.
%    compute_theo_img - Wether or not to compute the composite image of detected patterns
%    use_parallel_processing - If true, patterns are fitted in parallel  using all available workers (make sure there are some available)
%
% Output:
%    match_data  - Cell array of length '#Patterns'. The k's cell holds a
%                  M(k)x5 matrix for the corresponding pattern in patterns(:,:,k),
%                  where M(k) is the number of valid matches. The columns
%                  correspond to (rowShift, colShift, integrated intensity, background, correlation value).
%                  The (per-pixel) background is estimated as the median along the pattern border.
%                  The intensity is the summed up area of the pattern match corrected with the estimated background.
%                  NOTE: This only makes sense for patterns with low intensity on the border.
%                  For rowShift = colShift = 0 the upper left pixel of the
%                  pattern lies ontop the upper left pixel of the image.
%                  The pattern centers can be calculated as
%                  (rowShift + pattRows/2, colShift + pattCols/2).
%    im_detect   - Composite image showing the detected patterns.
%
% Note: Only patterns fully embedded in the image are detected.

%Note: CORR_THRESH depends on 'empty' space around the patterns.
% For larger zero padding and the presence of noise in the image
% this should be reduced, as the normalized correlation drops down.
% That said, zero padding seems to improve the performance of the detection
% (while sacrificing the usable border).


% NOTE NOTE: May use the fact that different parameters for an emitter
% result in a different intensity? -> likelihood changes

% NOTE: (x) Improve speed of normxcorr2 by precomputing the ffts of the image and the patterns 
%       for freqxcorr so that they may be reused. (currently the fft of the
%       image is computed again and again for every patterns correlation!)
%       (x) For processing movies the ffts of the patterns will also be
%       recomputed every time!

% Author: Simon Christoph Stein
% Date: 2014
% EMail: scstein@phys.uni-goettingen.de

imgRows = size(img,1);
imgCols = size(img,2);

pattRows = size(patterns,1);
pattCols = size(patterns,2);

nrPatterns = size(patterns,3);


% -- Pattern matching --

% Check if FFTs are supplied by user
if nargin < 5; FFTimg = []; end;
if nargin < 6; FFTpatterns = []; end;
if (nargin < 7 || isempty(compute_im_detect)); compute_im_detect = false; end
if (nargin < 8 || isempty(use_parallel_processing)); use_parallel_processing = false; end

% Match every single pattern individually and store the results
match_data = cell(nrPatterns,1);
match_images = zeros(imgRows-pattRows+1, imgCols-pattCols+1, nrPatterns);

if(use_parallel_processing)
    parfor iPatt = 1:nrPatterns
        if isempty(FFTpatterns)
            [ m_data, m_img ] = match_pattern(img, patterns(:,:,iPatt), CORR_THRESH, MIN_DIST, FFTimg);
        else
            [ m_data, m_img ] = match_pattern(img, patterns(:,:,iPatt), CORR_THRESH, MIN_DIST, FFTimg, FFTpatterns(:,:,iPatt));
        end    
        match_data{iPatt} = m_data;
        match_images(:,:,iPatt) = m_img;
    end
else
    for iPatt = 1:nrPatterns
        if isempty(FFTpatterns)
            [ m_data, m_img ] = match_pattern(img, patterns(:,:,iPatt), CORR_THRESH, MIN_DIST, FFTimg);
        else
            [ m_data, m_img ] = match_pattern(img, patterns(:,:,iPatt), CORR_THRESH, MIN_DIST, FFTimg, FFTpatterns(:,:,iPatt));
        end    
        match_data{iPatt} = m_data;
        match_images(:,:,iPatt) = m_img;
    end
end


% > Filter results: <
% Take the best match (over all patterns) within a box window of edge length 2*MIN_DIST+1.

% Project best matches into one image.
%   best_img is the image of maximal correlation values
%   idx_img corresponds to the pattern ID of the max of every pixel
[best_img, idx_img] = max(match_images,[], 3);


% > locate highest maximum within search window <
zero_mask = (best_img==0);

% (x) Using Jans localmax function
% halfw_col = MIN_DIST; % search window size
% halfw_row = MIN_DIST; % search window size
%filtered_best = localmax2d(best_img,(1+halfw_col):(size(best_img,2)+halfw_col),(1+halfw_row):(size(best_img,1)+halfw_row),false);

% (x) Using imdilate
tmp_idx = localMax(best_img, MIN_DIST);
filtered_best = zeros(size(best_img));
filtered_best(tmp_idx) = 1;

filtered_best(zero_mask) = 0;

% The filtered_best image is 1 if a pixel is the best (maximum correlation
% value) candidate in the search box and 0 otherwise.
% Remove the deleted (non-optimal) candidates by setting their values to zero.
idx_img( ~filtered_best) = 0;
best_img( ~filtered_best) = 0;


% -- Rebuild match data --
% (shifts and correlation values) with the remaining matches
for iPatt = 1:nrPatterns
    patt_idx = (idx_img == iPatt);
    ind = find(patt_idx); % indices of all matches for iPatt
    [max_rPos, max_cPos] = ind2sub(size(patt_idx), ind);
    max_shifts = [max_rPos, max_cPos] -1; % (rowShift, colShift) of maxima
    
    match_data{iPatt} =   [max_shifts, best_img(ind)]; % (rowShift, colShift, corrValue)
end

% Estimate background and intensity
for iPatt = 1:nrPatterns
    matches = match_data{iPatt};
    matches_A_BG = zeros(size(matches,1),5);
    for iMatch = 1:size(matches,1)
        rowShift = matches(iMatch, 1);
        colShift = matches(iMatch, 2);
        
        % Estimate background as median of the pattern border
        toprow    = img(rowShift+1, colShift+1:colShift+pattCols);
        bottomrow = img(rowShift+pattRows, colShift+1:colShift+pattCols);
        leftrow   = img(rowShift+2:rowShift+pattRows-1, colShift+1);
        rightrow  = img(rowShift+2:rowShift+pattRows-1, colShift+pattCols-1);
        BG = median([toprow(:); bottomrow(:); leftrow(:); rightrow(:)]);
        
        BGsum = pattRows*pattCols*BG;
        
        % Intensity is summed up region of the pattern minus the background
        Intensity = sum(sum(img(rowShift+1:rowShift+pattRows, colShift+1:colShift+pattCols),1),2) - BGsum;
        
        % rowShift, colShift, integrated intensity, background, correlation value
        matches_A_BG(iMatch, :) = [rowShift, colShift, Intensity, BG, matches(iMatch,3)];
    end
    
    match_data{iPatt} = matches_A_BG;
end


% -- Generate the composite image with identified patterns --

% Add the matches up to generate the composite image. (if requested)
% The pixel values of a patterns binary representation are set to its ID.
if compute_im_detect
    im_detect = zeros(imgRows, imgCols);
    for iPatt = 1:nrPatterns
        matches = match_data{iPatt};
        for iMatch = 1:size(matches,1)
            rowShift = matches(iMatch, 1);
            colShift = matches(iMatch, 2);
            Intensity = matches(iMatch,3);
            im_detect(rowShift+1:rowShift+pattRows, colShift+1:colShift+pattCols) = im_detect(rowShift+1:rowShift+pattCols, colShift+1:colShift+pattCols) + Intensity*patterns(:,:,iPatt);
        end
    end
else
    im_detect = NaN;
end


end



function max_indices = localMax(img, MIN_DIST)
% Returns indices of highest local maximum within a box window of edge
% length 2*MIN_DIST+1 around each pixel. Connected pixels with equal values
% are considered maxima as well.
    neigh = ones(2*MIN_DIST+1); % Neighborhood matrix
    img_dilated = imdilate(img, neigh);
    max_indices = find(img == img_dilated);
end


