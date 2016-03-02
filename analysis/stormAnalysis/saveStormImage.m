function [storm_img] = saveStormImage( posData, filename, movieSize, minval, maxval, mag)
% [storm_img] = saveStormImage( posData, filename_movie, imsize, minval, maxval, mag)
% Generates a STORM like image from the input data.
%
% Input:
%    posData: Either fitData or trajectoryData.
%    filename: filename to save to with "_storm.png" appended.
%    movieSize: Size of the original movie [width,height].
% Optional input: (use '[]' to leave unspecified)
%    minval: Value shown as black in the image.
%    maxval: Value shown as white in the image.
%    mag: Magnification of original pixels.
%         With mag=2 the generated pixels are half the size of the original
%         ones. | default: 8 (like RapidStorm)
%
% Output:
%    storm_img: The generated storm image

if nargin < 6 || isempty(mag)
    mag = 8; % 8x magnification is used by Rapidstorm
end

%% Generate STORM image
if(iscell(posData)) % This is fitData not trajectoryData
    positions = vertcat(posData{:});
    positions = positions(:,1:2);
else
    positions = posData(:,3:4);
end

% Visualize the data
% Set the histogram centers. NOTE: Our coordinates are integer at pixel
% centers, i.e. the most upper left pixels center is at (1,1). If you
% magnify by factor 2 (split 1 pixel into 4) the new center coordinates are
% at 0.75,1.25,1.75,2,25 (first two pixels) and so on.. Think about it!
centers = {0.5+1/(2*mag):1/mag:movieSize(1)+0.5, 0.5+1/(2*mag):1/mag:movieSize(2)+0.5};
storm_img = hist3(positions, centers);
storm_img = storm_img.'; % Image needs to be transposed .. its a MATLAB coordinate thing

if nargin<4 || isempty(minval)
    minval = min(storm_img(:));
end

if nargin<5 || isempty(maxval)
    maxval = max(storm_img(:));
end

h = figure;
imagesc(storm_img, [minval, maxval]);
axis image
colormap hot
colorbar
axis off

% % Saving to disk
if ~isempty(filename)
    [~, name, ~] = fileparts(filename);
    outname = [name '_storm.png'];
    fprintf('Saving image %s ..\n', outname);
    save_tiff(outname, storm_img, 'uint',16);
    save_img_col(outname, storm_img, minval, maxval);
end

end

