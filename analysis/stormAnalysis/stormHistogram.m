function [stormMap,stormRaw] = stormHistogram(inputData,movieSize,pixelSize,histogramType,filename,mag,nHistogramSteps,MLEenabled,locPrecision,minIntens,maxIntens)
% [stormMap] = stormHistogram(inputData,movieSize,pixelSize,histogramType,filename,mag,nHistogramSteps,MLEenabled,locPrecision,minIntens,maxIntens)
% Create STORM histogram from TrackNTrace data using various histogramming
% methods.
% 
% INPUT:
%     inputData: Either fittingData (cell array) or trackingData(2D double
%     array) which come out of TrackNTrace, OR TNT file name. In the latter
%     case, trackingData is used if available.
%     
%     movieSize: 1D integer array, [WIDTH,HEIGHT] of STORM movie in
%     [pixels].
%     
%     pixelSize: Double, pixel size in [nm].
%     
%     histogramType: String, either 'weighted', 'jitter', 'gauss' or
%     something else in which case the routine defaults to simply binning
%     localizations.
%     'WEIGHTED': All localizations are binned with a weight of 1/precision
%     (the localization precision) instead of 1.
%     'JITTER': Draw noisy positions out of a Gaussian distribution where
%     the means are the positions and the standard deviations are the
%     localization precision. These positions are binned and an average of
%     <N> jittered histograms is calculated, default: N=4.
%     'GAUSS': Add Gaussian PSFs with an amplitude of 1 and a standard
%     deviation equal to the localization precision to the image instead if
%     binning positions.
%     
%     filename: String, name of a file where the results and a PNG image
%     are saved to. Can be empty.
%     
%     mag: Double, magnification or 'zoom' factor of STORM histogram.
%     
%     nHistogramSteps: Integer, number of steps used for jittered
%     histogram.
%     
%     MLEenabled: Boolean, set to true if MLE fitting was used.
%     
%     locPrecision: Double, overrides localization precision to a fixed
%     value. Must be in [nm].
%     
%     minIntens: Double in [0,1] range, sets minimum contrast of STORM
%     image.
%     
%     maxIntens: Double in [0,1] range, sets maximum contrast of STORM
%     image.
%     
%     
% OUTPUT:
%     stormMap: Normalized STORM histogram.
%     
%     stormRaw: Non-normalized STORM histogram.

%% Parse input

if nargin < 4
    filename = [];
end

if isempty(mag) || nargin < 5
    mag = 8;
end
mag = round(mag);

if isempty(nHistogramSteps) || nargin < 6
    nHistogramSteps = 4;
end

if isempty(MLEenabled) || nargin < 7
    MLEenabled = false;
end

if nargin < 8
    locPrecision = [];
end

guessLocPrecision = false;
if MLEenabled
    guessLocPrecision = true;
end

if isempty(minIntens) || isempty(maxIntens) || nargin < 9
    minIntens = 0;
    maxIntens = 0.8;
end
minIntens1 = min([minIntens,maxIntens]);
maxIntens = max([minIntens,maxIntens]);
minIntens = minIntens1; clear minIntens1;


%% Filter data

% Parse TNT data, either trackingData or fittingData
if ischar(inputData)
    try
        load(inputData,'trackingData');
        fittingData = inputData(:,3:end);
    catch
        load(inputData,'fittingData');
        fittingData = vertcat(inputData{:});
    end
else
    if iscell(inputData)
        %... then it's fittingData
        fittingData = vertcat(inputData{:});
    else
        fittingData = inputData(:,3:end);
    end
end

% Calculate localization precision if necessary
pos = fittingData(:,1:2)-0.5; %shift to: upper left pixel center = (0.5,0.5)
if guessLocPrecision
    % !!locPrecision is given in pixels!!
    if size(fittingData,2)<7
        sigma_sq = fittingData(:,6).^2;
    else
        sigma_sq = mean(fittingData(:,6)+fittingData(:,7),2).^2;
    end
    N = fittingData(:,4).*2*pi.*sigma_sq; %total number of photons
    tau = 2*pi*fittingData(:,5).*(sigma_sq+1/12)./N;
    locPrecision = sqrt((sigma_sq+1/12)./N.*(1+4*tau+sqrt(2*tau./(1+4*tau)))); %Rieger et al, DOI 10.1002/cphc.201300711
else
    if isempty(locPrecision)
        disp(sprintf('Localization precision estimation not possible, switchting to default STORM histogram.\n'));
        histogramType = 'default';
    else
        locPrecision = repmat(locPrecision/pixelSize,size(pos,1),1);
    end
end

%% Calculate STORM histogram

switch(histogramType)
    case 'weighted'
        % weighted: create histogram with weights 1/locPrecision instead of 1
        valid_pos = ~isinf(locPrecision)&~isnan(locPrecision);
        pos = pos(valid_pos,:);
        locPrecision = locPrecision(valid_pos,:);
        stormRaw = createHistogram(pos,1/mag,1./locPrecision);
        
    case 'jitter'
        % create jittered histogram where an average histogram of noisy positions is calculated
        valid_pos = ~isinf(locPrecision)&~isnan(locPrecision);
        pos = pos(valid_pos,:);
        locPrecision = locPrecision(valid_pos,:);
        
        pos_jitter = pos+randn(size(pos,1),2).*[locPrecision,locPrecision];
        stormRaw = createHistogram(pos_jitter,1/mag,ones(size(pos_jitter,1),1));
        
        for i=2:nHistogramSteps
            pos_jitter = pos+randn(size(pos,1),2).*[locPrecision,locPrecision];
            stormRaw = stormRaw+createHistogram(pos_jitter,1/mag,ones(size(pos_jitter,1),1));
        end
        stormRaw = stormRaw/nHistogramSteps;
        
    case 'gauss'
        % add gaussian functions with sigma = locPrecision to empty image
        valid_pos = ~isinf(locPrecision)&~isnan(locPrecision);
        pos = pos(valid_pos,:);
        locPrecision = locPrecision(valid_pos,:);
        
        pos_idx = floor(pos*mag)+1;
        stormRaw = zeros(movieSize(2)*mag,movieSize(1)*mag);
        for iPos = 1:size(pos,1)
            sigma = locPrecision(iPos)*mag;
            halfw = ceil(3*sigma);
            if pos_idx(iPos,1)-halfw<1 || pos_idx(iPos,1)+halfw>size(stormRaw,2) || pos_idx(iPos,2)-halfw<1 || pos_idx(iPos,2)+halfw>size(stormRaw,1)
                % out of bounds?
                continue
            end
            [img] = gaussianMask(rem(pos(iPos,1),1),rem(pos(iPos,2),1),sigma,sigma,halfw);
            stormRaw(pos_idx(iPos,2)-halfw:pos_idx(iPos,2)+halfw,pos_idx(iPos,1)-halfw:pos_idx(iPos,1)+halfw) = stormRaw(pos_idx(iPos,2)-halfw:pos_idx(iPos,2)+halfw,pos_idx(iPos,1)-halfw:pos_idx(iPos,1)+halfw)...
                +img;
        end
        
    otherwise
        % default to simple count histogram
        stormRaw = createHistogram(pos,1/mag,ones(size(pos,1),1));
        
end %SWITCH

%% Plot and save normalized image
stormMap = (stormRaw-min(stormRaw(:)))./(max(stormRaw(:))-min(stormRaw(:)));

h=figure;
x = ((1:size(stormMap,2))-0.5)/mag*pixelSize;
y = ((1:size(stormMap,1))-0.5)/mag*pixelSize;
imagesc(x,y,stormMap,[minIntens,maxIntens]); colormap hot; axis image; colorbar;
xlabel('x [nm]');
ylabel('y [nm]');

if ~isempty(filename)
    [folder, name, ~] = fileparts(filename);
    outname = [folder,filesep,name '_storm.png'];
    fprintf('Saving image %s ..\n', outname);
    print(h,outname,'-dpng','-r300');
    save([outname(1:end-4),'_raw.mat'],'stormRaw','stormMap','x','y');
end

%% Helper functions

    function [histogram] = createHistogram(XY,binWidth,XYweights)
        out_of_bounds = XY(:,1)<0 | XY(:,1)>movieSize(2) | XY(:,2)<0 | XY(:,2)>movieSize(1);
        XY = XY(~out_of_bounds,:);
        XYweights = XYweights(~out_of_bounds,:);
        idx = floor(XY/binWidth)+1;
        idx = sub2ind(movieSize*mag,idx(:,2),idx(:,1));
        
        histogram = zeros(movieSize(2)*mag,movieSize(1)*mag);
        for jPos = 1:numel(idx)
            histogram(idx(jPos)) = histogram(idx(jPos))+XYweights(jPos);
        end
        
    end


    function [img] = gaussianMask(x,y,sigma_x,sigma_y,halfw)
        [x_grid,y_grid] = meshgrid (-halfw:halfw,-halfw:halfw);
        
        img = 1/2*(erfc(-(x_grid-x+0.5)/(sqrt(2)*sigma_x))-erfc(-(x_grid-x-0.5)/(sqrt(2)*sigma_x))).*...
            (erfc(-(y_grid-y+0.5)/(sqrt(2)*sigma_y))-erfc(-(y_grid-y-0.5)/(sqrt(2)*sigma_y)));
    end


end %MAIN FUNCTION






