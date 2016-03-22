% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
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
function [img_out] = filterImageGeneral(img,filterKernel,isSeparable,shape)
% [IMG_OUT] = filterImageGeneral(IMG,KERNEL,SEPARABLE,SHAPE)
% Filter a 2D image with a 1D or 2D filter using conv2.
% 
% INPUT:
%     IMG: 2D array of pixel intensities.
%     
%     KERNEL: If not a cell, it's assumed to be either a 1D row or a 2D
%     symmetric array which is parsed to conv2 together with the image. If
%     it is a cell, then the function will return the difference of the
%     images filtered with the first and the second cell array entry where
%     the second entry is assumed to be of larger size. In this case, the
%     outcoming filtered image will always have the same dimensions as
%     conv2(IMG,KERNEL_2,'valid').
%     
%     SEPARABLE: Boolean. If a filter is separable, that is, it can be
%     written as K = k'*k, filtering can be calculated as 2x1D convolutions
%     which is faster in most cases than a full 2D convolution. KERNEL has
%     to be a row vector in this case or a cell array thereof.
%     
%     SHAPE: String. See conv2, can be 'full', 'same', or 'valid'. For
%     subtraction filtering the shape is always 'valid'.
%     
% OUTPUT:
%     IMG_OUT: Filtered 2D array with dimensions according to SHAPE. See
%     conv2 for details.


if isSeparable
    if iscell(filterKernel)
        %in this case we apply 2 filters of different size, ignoring
        %shape
        padding_valid = size(filterKernel{2},1)/2;
        img_out = conv2(conv2(img,filterKernel{1},'same'),filterKernel{1}.','same') - ...
            conv2(conv2(img,filterKernel{2},'same'),filterKernel{2}.','same');
        img_out = img_out(ceil(padding_valid):end-floor(padding_valid),ceil(padding_valid):end-floor(padding_valid));
    else
        img_out = conv2(conv2(img,filterKernel,shape),filterKernel.',shape);
    end
else
    if iscell(filterKernel)
        %in this case we apply 2 filters of different size, ignoring
        %shape
        padding_valid = size(filterKernel{2},1)/2;
        img_out = conv2(img,filterKernel{1},'same')-conv2(img,filterKernel{2},'same');
        img_out = img_out(ceil(padding_valid):end-floor(padding_valid),ceil(padding_valid):end-floor(padding_valid));
    else
        img_out = conv2(img,filterKernel,shape);
    end
end

