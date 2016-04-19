function [stormMap,stormRaw] = stormHistogram(inputData,movieSize,pixelSize,histogramType,filename,superResMag,nHistogramSteps,photonConversion,locPrecision,minIntens,maxIntens)
% [stormMap,stormRaw] = stormHistogram(inputData,movieSize,pixelSize,histogramType,filename,superResMag,nHistogramSteps,photonConversion,locPrecision,minIntens,maxIntens)
% Easier call: [stormMap,stormRaw]
% Create STORM histogram from TrackNTrace data using various histogramming
% methods.
% 
% INPUT:
%     inputData: Either refinementData (cell array) or trackingData(2D double
%     array) which come out of TrackNTrace, OR TNT file name. In the latter
%     case, trackingData is used if available.
%     
%     movieSize: 1D integer array, [WIDTH,HEIGHT] of STORM movie in
%     [pixels].
%     
%     pixelSize: Pixel size in [nm].
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
%     deviation equal to the localization precision to the isuperResMage instead if
%     binning positions.
%     
%     filename: String, name of a file where the results and a PNG isuperResMage
%     are saved to. Can be empty.
%     
%     superResMag: Magnification or 'zoom' factor of STORM histogram.
%     
%     nHistogramSteps: Integer, number of steps used for jittered
%     histogram.
%     
%     photonConversion: Boolean, set to true if ADU counts were converted
%     to photons.
%     
%     locPrecision: Overrides localization precision to a fixed
%     value. Must be in [nm].
%     
%     minIntens: Float in [0,1] range, sets minimum contrast of STORM
%     isuperResMage.
%     
%     maxIntens: Float in [0,1] range, sets maximum contrast of STORM
%     isuperResMage.
%     
%     
% OUTPUT:
%     stormMap: Normalized STORM histogram.
%     
%     stormRaw: Non-normalized STORM histogram.

%% Parse input

% No input data? Try to do things automatically as best as possible
if nargin < 4 || isempty(inputData)
    [inputData, path] = uigetfile({'*.mat*', 'MATLAB TrackNTrace file'},'Choose TrackNTrace file to analyze.');
    inputData = [path,inputData];
    load(inputData,'movieSize'); movieSize = movieSize(1:2);
    
    if ~exist('pixelSize','var') || isempty(pixelSize)
        warning off backtrace
        warning('Pixel size not given, setting to 100 nm.');
        warning on backtrace
        pixelSize = 100;
    end
    
    if ~exist('histogramType','var') || isempty(histogramType)
        warning off backtrace
        warning('Histogram type not given, setting to ''gauss''.');
        warning on backtrace
        histogramType = 'gauss';
    end
    
    load(inputData,'globalOptions');
    if exist('globalOptions','var')
        if isfield(globalOptions,'usePhotonConversion;')
            photonConversion = globalOptions.usePhotonConversion;
        end
    end
    
    load(inputData,'filename_movie');
    if exist('filename_movie','var')
        filename = filename_movie;
    end
        
end

    

if ~exist('filename','var') || isempty(filename)
    filename = [];
end

if ~exist('superResMag','var') || isempty(superResMag) 
    superResMag = 8;
end
superResMag = round(superResMag);

if ~exist('nHistogramSteps','var') || isempty(nHistogramSteps)
    nHistogramSteps = 4;
end

if ~exist('photonConversion','var') || isempty(photonConversion) 
    photonConversion = false;
end

if ~exist('locPrecision','var')
    locPrecision = [];
end

guessLocPrecision = false;
if photonConversion
    guessLocPrecision = true;
end

if ~exist('minIntens','var') || ~exist('maxIntens','var') || isempty(minIntens) || isempty(maxIntens)
    minIntens = 0;
    maxIntens = 0.6;
end
minIntens1 = min([minIntens,maxIntens]);
maxIntens = max([minIntens,maxIntens]);
minIntens = minIntens1; clear minIntens1;


%% Filter data

% Parse TNT data, either trackingData or refinementData
if ischar(inputData)
    load(inputData,'trackingData','refinementData');
    if exist('trackingData','var');
        refinementData = trackingData(:,3:end); %#ok<NODEF>
        clear trackingData;
    else
        refinementData = vertcat(refinementData{:}); %#ok<NODEF>
    end
else
    if iscell(inputData)
        %... then it's refinementData
        refinementData = vertcat(inputData{:});
    else
        refinementData = inputData(:,3:end);
    end
end

% Calculate localization precision if necessary
pos = refinementData(:,1:2)-0.5; %shift to: upper left pixel center = (0.5,0.5)
if guessLocPrecision
    % !!locPrecision is given in pixels!!
    if size(refinementData,2)<7
        sigma_sq = refinementData(:,6).^2;
    else
        sigma_sq = mean(refinementData(:,6)+refinementData(:,7),2).^2;
    end
    N = refinementData(:,4).*2*pi.*sigma_sq; %total number of photons
    tau = 2*pi*refinementData(:,5).*(sigma_sq+1/12)./N;
    locPrecision = sqrt((sigma_sq+1/12)./N.*(1+4*tau+sqrt(2*tau./(1+4*tau)))); %Rieger et al, DOI 10.1002/cphc.201300711
else
    if isempty(locPrecision)
        disp(sprintf('Localization precision estimation not possible, switchting to default STORM histogram.\n')); %#ok<DSPS>
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
        stormRaw = createHistogram(pos,1/superResMag,1./locPrecision);
        
    case 'jitter'
        % create jittered histogram where an average histogram of noisy positions is calculated
        valid_pos = ~isinf(locPrecision)&~isnan(locPrecision);
        pos = pos(valid_pos,:);
        locPrecision = locPrecision(valid_pos,:);
        
        pos_jitter = pos+randn(size(pos,1),2).*[locPrecision,locPrecision];
        stormRaw = createHistogram(pos_jitter,1/superResMag,ones(size(pos_jitter,1),1));
        
        for i=2:max(2,nHistogramSteps)
            pos_jitter = pos+randn(size(pos,1),2).*[locPrecision,locPrecision];
            stormRaw = stormRaw+createHistogram(pos_jitter,1/superResMag,ones(size(pos_jitter,1),1));
        end
        stormRaw = stormRaw/nHistogramSteps;
        
    case 'gauss'
        % add gaussian functions with sigma = locPrecision to empty isuperResMage
        valid_pos = ~isinf(locPrecision)&~isnan(locPrecision);
        pos = pos(valid_pos,:);
        locPrecision = locPrecision(valid_pos,:);
        
        pos_idx = floor(pos*superResMag)+1;
        stormRaw = zeros(movieSize(2)*superResMag,movieSize(1)*superResMag);
        for iPos = 1:size(pos,1)
            sigma = locPrecision(iPos)*superResMag;
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
        stormRaw = createHistogram(pos,1/superResMag,ones(size(pos,1),1));
        
end %SWITCH

%% Plot and save normalized isuperResMage
stormMap = (stormRaw-min(stormRaw(:)))./(max(stormRaw(:))-min(stormRaw(:)));

h=figure;
x = ((1:size(stormMap,2))-0.5)/superResMag*pixelSize;
y = ((1:size(stormMap,1))-0.5)/superResMag*pixelSize;
imagesc(x,y,stormMap,[minIntens,maxIntens]); colormap hot; axis image; colorbar;
xlabel('x [nm]');
ylabel('y [nm]');

if ~isempty(filename)
    [folder, name, ~] = fileparts(filename);
    if isempty(folder)
        outname = [name '_storm.png'];
    else
        outname = [folder,filesep,name '_storm.png'];
    end
    fprintf('Saving isuperResMage %s ..\n', outname);
    print(h,outname,'-dpng','-r300');
    save([outname(1:end-4),'_raw.mat'],'stormRaw','stormMap','x','y');
end

%% Helper functions

    function [histogram] = createHistogram(XY,binWidth,XYweights)
        out_of_bounds = XY(:,1)<0 | XY(:,1)>movieSize(2) | XY(:,2)<0 | XY(:,2)>movieSize(1);
        XY = XY(~out_of_bounds,:);
        XYweights = XYweights(~out_of_bounds,:);
        idx = floor(XY/binWidth)+1;
        idx = sub2ind(movieSize*superResMag,idx(:,2),idx(:,1));
        
        histogram = zeros(movieSize(2)*superResMag,movieSize(1)*superResMag);
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






