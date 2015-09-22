function [xbins,y,hx] = Hist1D(X,h0,pic,bin,center,zero,hint)
% [xbins,y,hx] = HIST1D(X,h0,pic,bin,center,zero,hint)

% Constructs a 1D histogram out of the vector X with various binning rules.

% INPUT:
	% X: vector array of data to be binned
	% h0: double, bin size. Can be set to 0 for automatic binning (recommended)
	% pic: boolean, plot histogram
	% bin: string, sets binning rule. Available options are Freedman-Diaconis, Sturges rule, Scott's rule, and the standard Microsoft Excel rule
	% center: boolean, set to true for normaling binning and false for edge binning (using histc)
	% zero: center histogram symmetrically around 0
	% hint: boolean, set to true to obtain rounded bins. Useful for counting whole numbers (pixels, photons, ...)

% OUTPUT:
	% xbins: vector of bins
	% y: frequency of bins
	% hx: bin size
	

X = X(:);
Xmin = min(X);
Xmax = max(X);

if hint
    Xmin = round(Xmin);
    Xmax = round(Xmax);
end

switch bin 
    case 'freedman'
        hx = 2*iqr(timeseries(X))/length(X)^(1/3); %Freedman-Diaconis
        k = ceil((Xmax-Xmin)/hx);
    case 'sturges'
        k = ceil(log2(length(X))+1); %Sturges
        hx = (Xmax-Xmin)/k;
    case 'scott' %Scott
        hx = 3.5*std(X)/length(X)^(1/3);
        k = ceil((Xmax-Xmin)/hx);
    otherwise %Excel, Origin?, ...
        k = round(sqrt(length(X)));
        hx = (Xmax-Xmin)/k;
end
        
if h0>0 %set minimal bin distance
    hx = h0;
    k = ceil((Xmax-Xmin)/hx);
end

if hint
    hx = round(hx);
    k = ceil((Xmax-Xmin)/hx);
end
    

xbins = Xmin:hx:Xmin+k*hx; %in every case, this will be one bin too many
if rem(Xmax-max(xbins),2)>0
    xbins = xbins(1:end-1); %thus, cut off last bin ...
end
xbins = xbins+(Xmax-max(xbins))/2; % ... and shift the bins so that the center is the central bin (odd) or between both central bins (even)
if zero
    edge = max(abs(Xmin),Xmax);
    xbins = 0:hx:edge;
    xbins = [-xbins(end:-1:2),xbins];
end

%count events. Note that hist, histc counts out-of-boundary values to first
%or last bin respectively
if center
    y = hist(X(:),xbins(:));
else
    y = histc(X(:),xbins(:)); %histc is used for counting edges
end

if pic
    figure;
    bar(xbins,y);
end