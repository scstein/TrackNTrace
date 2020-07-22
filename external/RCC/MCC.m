
% -------------------------------------------------------------------------
% Matlab version of Mean cross-correction method of OE2012 paper
% Input:    coords:             localization coordinates [x y frame], 
%           segpara:            segmentation parameters (time wimdow, frame)
%           imsize:             image size (pixel)
%           pixelsize:          camera pixel size (nm)
%           binsize:            spatial bin pixel size (nm)
% Output:   coordscorr:         localization coordinates after correction [xc yc] 
%           finaldrift:         drift curve (save b matrix for other analysis)
% By Yina Wang @ Hust 2013.10.21
% -------------------------------------------------------------------------

function [coordscorr, finaldrift, bsave] = MCC(coords, segpara, imsize, pixelsize, binsize )
    
zoomfactor = pixelsize / binsize;
ntotalframe = max(coords(:,3));

coordscorr = zeros(size(coords,1),size(coords,2));
finaldrift = zeros(ntotalframe,2);

%% bin localizations, frame to frame drift calculation
nbinframe = floor(ntotalframe / segpara);

imshift = zeros(nbinframe * (nbinframe - 1)/2,2); % shifts between any two bin steps in x and y direction
drift=zeros(nbinframe,2);

flag=1;
for i = 1:nbinframe-1
    index = coords(:,3)>(i-1)*segpara+1 & coords(:,3)<=i*segpara;
    imorigin = BinLocalizations(coords(index,1:2), imsize, zoomfactor);
    autocorr = CrossCorrelation(imorigin,imorigin);
    [yorigin, xorigin] = GaussianFit(autocorr);
    
    for j = i+1:nbinframe
        index = coords(:,3)>(j-1)*segpara+1 & coords(:,3)<=j*segpara;
        imbin = BinLocalizations(coords(index,1:2), imsize, zoomfactor);
        corr = CrossCorrelation(imorigin,imbin);
        [y,x] = GaussianFit(corr);
        
        imshift(flag,1) = xorigin-x;
        imshift(flag,2) = yorigin-y;
        drift(i,1)=drift(i,1)-imshift(flag,1);
        drift(i,2)=drift(i,2)-imshift(flag,2);
        flag=flag+1;
    end
    
    clear index
    for j = 1:i-1
        index=(nbinframe-1)*(j-1)-(j-2)*(j-1)/2-j+i;
        drift(i,1)=drift(i,1)+imshift(index,1);
        drift(i,2)=drift(i,2)+imshift(index,2);
    end
end
for j = 1:nbinframe-1
    index=(nbinframe-1)*(j-1)-(j-2)*(j-1)/2-j+nbinframe;
    drift(nbinframe,1)=drift(nbinframe,1)+imshift(index,1);
    drift(nbinframe,2)=drift(nbinframe,2)+imshift(index,2);
end
drift=drift./nbinframe;

%% save results
bsave=imshift;

%% drift interpolation
indexinterp=zeros(nbinframe+2,1);
indexinterp(1)=1;
indexinterp(nbinframe+2)=ntotalframe;
indexinterp(2:nbinframe+1)=round(segpara/2):segpara:segpara*nbinframe-1;

drift(:,1)=drift(:,1)-drift(1,1);
drift(:,2)=drift(:,2)-drift(1,2);
finaldrift(:,1) = interp1(indexinterp,[0 drift(:,1)' drift(nbinframe,1)],1:ntotalframe,'spline');
finaldrift(:,2) = interp1(indexinterp,[0 drift(:,2)' drift(nbinframe,2)],1:ntotalframe,'spline');
finaldrift = finaldrift./zoomfactor;

for i = 1:ntotalframe
    index = find(coords(:,3)==i);
    coordscorr(index,1) = coords(index,1)-finaldrift(i,1);
    coordscorr(index,2) = coords(index,2)-finaldrift(i,2);
end
coordscorr(:,3) = coords(:,3);

end
