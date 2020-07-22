% -------------------------------------------------------------------------
% Matlab version of direct cross-correlation drift correction method
% Input:    coords:             localization coordinates [x y frame], 
%           segpara:            segmentation parameters (time wimdow, frame)
%           imsize:             image size (pixel)
%           pixelsize:          camera pixel size (nm)
%           binsize:            spatial bin pixel size (nm)
% Output:   coordscorr:         localization coordinates after correction [xc yc] 
%           finaldrift:         drift curve 
% By Yina Wang @ Hust 2013.10.18
% -------------------------------------------------------------------------

function [coordscorr, finaldrift] = DCC(coords, segpara, imsize, pixelsize, binsize )

zoomfactor = pixelsize / binsize;
ntotalframe = max(coords(:,3));

coordscorr = zeros(size(coords,1),size(coords,2));
finaldrift = zeros(ntotalframe,2);

%% bin localizations, frame shifts calculation
nbinframe = floor(ntotalframe / segpara);
    
imshift = zeros(nbinframe-1,1);
fitresult = zeros(5,nbinframe-1);

i=1;
index = coords(:,3)>(i-1)*segpara+1 & coords(:,3)<=i*segpara;
imorigin = BinLocalizations(coords(index,1:2), imsize, zoomfactor);
autocorr = CrossCorrelation(imorigin,imorigin);
[yorigin, xorigin] = GaussianFit(autocorr);

for i = 2:nbinframe
    index = coords(:,3)>(i-1)*segpara+1 & coords(:,3)<=i*segpara;
    imbim = BinLocalizations(coords(index,1:2), imsize, zoomfactor);
        
    corr = CrossCorrelation(imorigin,imbim);
    [y,x] = GaussianFit(corr);
    imshift(i-1,1) = xorigin-x;
    imshift(i-1,2) = yorigin-y;
end

%% drift interpolation
indexinterp=zeros(nbinframe+2,1);
indexinterp(1)=1;
indexinterp(nbinframe+2)=ntotalframe;
indexinterp(2:nbinframe+1)=round(segpara/2):segpara:segpara*nbinframe-1;

finaldrift(:,1) = interp1(indexinterp,[0 0 imshift(:,1)' imshift(nbinframe-1,1)],1:ntotalframe,'spline');
finaldrift(:,2) = interp1(indexinterp,[0 0 imshift(:,2)' imshift(nbinframe-1,2)],1:ntotalframe,'spline');
finaldrift = finaldrift./zoomfactor;

for i = 1:ntotalframe
    index = find(coords(:,3)==i);
    coordscorr(index,1) = coords(index,1)-finaldrift(i,1);
    coordscorr(index,2) = coords(index,2)-finaldrift(i,2);
end
coordscorr(:,3) = coords(:,3);

end

