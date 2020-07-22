
% -------------------------------------------------------------------------
% Matlab version of 2D redundant drift correction method
% Input:    coords:             localization coordinates [x y frame], 
%           segpara:            segmentation parameters (time wimdow, frame)
%           imsize:             image size (pixel)
%           pixelsize:          camera pixel size (nm)
%           binsize:            spatial bin pixel size (nm)
%           rmax:               error threshold for re-calculate the drift (pixel)
% Output:   coordscorr:         localization coordinates after correction [xc yc] 
%           finaldrift:         drift curve (save A and b matrix for other analysis)
% By Yina Wang @ Hust 2013.09.09
% -------------------------------------------------------------------------

function [coordscorr, finaldrift, A,b] = RCC(coords, segpara, imsize, pixelsize, binsize, rmax)
    
zoomfactor = pixelsize / binsize;
ntotalframe = max(coords(:,3));

coordscorr = zeros(size(coords,1),size(coords,2));
finaldrift = zeros(ntotalframe,2);

%% bin localizations, frame to frame drift calculation
nbinframe = floor(ntotalframe / segpara);

imshift = single(zeros(nbinframe * (nbinframe - 1)/2,2)); % shifts between any two bin steps in x and y direction
A = single(zeros( nbinframe*(nbinframe-1)/2, nbinframe-1 )); %set equations

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
        A(flag,i:j-1) = 1;
        flag=flag+1;
    end
end


%% solve overdetermined equations and refine the calculated drift
drift = pinv(A)*imshift;
err = A*drift-imshift;
b=imshift;

for i=1:nbinframe*(nbinframe-1)/2
    rowerr(i,1) = sqrt(err(i,1)^2+err(i,2)^2);
end
rowerr(:,2)=1:1:nbinframe*(nbinframe-1)/2;
rowerr=flipud(sortrows(rowerr,1));
clear index
index=rowerr(find(rowerr(:,1)>rmax*zoomfactor),2);
for i=1:size(index,1)
    flag = index(i);
    tmp=A;
    tmp(flag,:)=[];
    if rank(tmp,1)==(nbinframe-1)
        A(flag,:)=[];
        b(flag,:)=[];
        sindex=find(index>flag);
        index(sindex)=index(sindex)-1;
    else
        tmp=A;
    end
end
drift1 = pinv(A)*b;
for i =1:size(drift,1)
drift(i,1)=sum(drift1(1:i,1));
drift(i,2)=sum(drift1(1:i,2));
end

%% drift interpolation
indexinterp=zeros(nbinframe+2,1);
indexinterp(1)=1;
indexinterp(nbinframe+2)=ntotalframe;
indexinterp(2:nbinframe+1)=round(segpara/2):segpara:segpara*nbinframe-1;

finaldrift(:,1) = interp1(indexinterp,[0 0 drift(:,1)' drift(nbinframe-1,1)],1:ntotalframe,'spline');
finaldrift(:,2) = interp1(indexinterp,[0 0 drift(:,2)' drift(nbinframe-1,2)],1:ntotalframe,'spline');
finaldrift = finaldrift./zoomfactor;

for i = 1:ntotalframe
    index = find(coords(:,3)==i);
    coordscorr(index,1) = coords(index,1)-finaldrift(i,1);
    coordscorr(index,2) = coords(index,2)-finaldrift(i,2);
end
coordscorr(:,3) = coords(:,3);

end