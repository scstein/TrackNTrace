function [bgMean,bgStd] = spatialMovAveBG(imageLast5,imageSizeX,imageSizeY)

%the function in its current form assigns blocks of 11x11 pixels the
%same background values, for the sake of speed
%
% Copyright (C) 2014 LCCB 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%define pixel limits where moving average can be calculated
startPixelX = 16;
endPixelX = max(imageSizeX - 15,startPixelX);
startPixelY = 16;
endPixelY = max(imageSizeY - 15,startPixelY);

%allocate memory for output
bgMean = NaN(imageSizeX,imageSizeY);
bgStd = bgMean;

%go over all pixels within limits
for iPixelX = startPixelX : 11 : endPixelX
    for iPixelY = startPixelY : 11 : endPixelY
        
        %get local image
        imageLocal = imageLast5(iPixelX-15:min(iPixelX+15,imageSizeX),iPixelY-15:min(iPixelY+15,imageSizeY),:);
           
        %estimate robust mean and std
        %first remove NaNs representing cropped regions
        imageLocal = imageLocal(~isnan(imageLocal));
        if ~isempty(imageLocal)
            [bgMean1,bgStd1] = robustMean(imageLocal(:));
            bgStd1 = max(bgStd1,eps);
        else
            bgMean1 = NaN;
            bgStd1 = NaN;
        end
        
        %put values in matrix representing image
        bgMean(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5) = bgMean1;
        bgStd(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5) = bgStd1;
        
    end
end

%find limits of actual pixels filled up above
% firstFullX = find(~isnan(bgMean(:,startPixelY)),1,'first');
% lastFullX = find(~isnan(bgMean(:,startPixelY)),1,'last');
% firstFullY = find(~isnan(bgMean(startPixelX,:)),1,'first');
% lastFullY = find(~isnan(bgMean(startPixelX,:)),1,'last');
firstFullX = startPixelX - 5;
lastFullX = iPixelX + 5;
firstFullY = startPixelY - 5;
lastFullY = iPixelY + 5;

%patch the rest
for iPixelY = firstFullY : lastFullY
    bgMean(1:firstFullX-1,iPixelY) = bgMean(firstFullX,iPixelY);
    bgMean(lastFullX+1:end,iPixelY) = bgMean(lastFullX,iPixelY);
    bgStd(1:firstFullX-1,iPixelY) = bgStd(firstFullX,iPixelY);
    bgStd(lastFullX+1:end,iPixelY) = bgStd(lastFullX,iPixelY);
end
for iPixelX = 1 : imageSizeX
    bgMean(iPixelX,1:firstFullY-1) = bgMean(iPixelX,firstFullY);
    bgMean(iPixelX,lastFullY+1:end) = bgMean(iPixelX,lastFullY);
    bgStd(iPixelX,1:firstFullY-1) = bgStd(iPixelX,firstFullY);
    bgStd(iPixelX,lastFullY+1:end) = bgStd(iPixelX,lastFullY);
end

