function [ fitData ] = locateParticles( movie, img_dark, candidateOptions, fittingOptions)
% [ fitData ] = locateParticles( movieStack, imgCorrection, candidateOptions, fittingOptions)
% Find locations of bright spots in an image, in this case a movie of
% fluorescent molecules, and fit a Gaussian distribution to these spots to
% obtain position, amplitude, background level and standard deviation.
%
% INPUT:
%     movieStack: 3D array of intensity values (y,x,N) where N is the
%     number of frames. All images are treated piecewise and only converted
%     to double when needed to avoid memory overlow
%
%     imgCorrection: 2D array of intensity correction values which is added
%     to each movie image to correct for non-isotropic camera readout.
%     Leave empty [] if no correction is needed
%
%     candidateOptions: struct of input options used to find localization
%     candidates. See RunTrackNTrace.m for details
%
%     fittingOptions: struct of input options used to fit localization
%     candidates. See RunTrackNTrace.m for details
%
%
% OUTPUT:
%     fitData: 1D cell array of of Gaussian distribution parameters for all
%     found particles with one cell per frame. The column order in one cell
%     is [x,y,A,B,sigma,flag], the line order is the particle index. The
%     positions x (img column) and y (img row) are not corrected by a middle pixel
%     shift, the center of the top left pixel is [1.0,1.0]. A is the
%     unnormalized amplitude of a Gaussian distribution A*exp(...)+B with
%     constant background B. flag is the exit flag of the fitting routine,
%     see psfFit_Image.m for details.

% Parse inputs
global img_bck_mean img_bck_std %global variable is needed for intensity weighted filtering to improve performance as background is only calculated every bck_interval frames
img_bck_mean = zeros(size(movie,1),size(movie,2));
img_bck_std = zeros(size(movie,1),size(movie,2));
bck_interval = candidateOptions.backgroundCalculationInterval;


fit_forward = candidateOptions.fitForward;
fit_sigma = fittingOptions.fitSigma;
usePixelIntegratedFit = fittingOptions.usePixelIntegratedFit;

useMLE = fittingOptions.useMLErefine;
if useMLE
    bias = fittingOptions.photonBias;
    photon_factor = fittingOptions.photonSensitivity/fittingOptions.photonGain;
end

calc_once = candidateOptions.calculateCandidatesOnce;
cand_corr = candidateOptions.useCrossCorrelation;
if ~isempty(img_dark) %if correction image is provided, do it
    correct_dark = true;
else
    correct_dark = false;
end
nrFrames = size(movie,3);
halfw = round(4*candidateOptions.sigma); %fitting routine window


% Trace first frame
if calc_once %if candidate search only takes place once, an average image is calculated to achieve higher SNR
    if fit_forward
        img = mean(double(movie(:,:,1:candidateOptions.averagingWindowSize)),3);
    else
        img = mean(double(movie(:,:,end:-1:end-candidateOptions.averagingWindowSize+1)),3);
    end
else
    if fit_forward
        img = double(movie(:,:,1));
    else
        img = double(movie(:,:,end));
    end
end

if correct_dark %use dark image stack to correct image for camera artifacts
    img = img+img_dark;
end

if useMLE
    if ~correct_dark
        img = img-bias;
    end
    img = img*photon_factor;
end


if cand_corr %find candidates either by cross correlation or intensity filtering
    candidatePos = findSpotCandidates(img, candidateOptions.sigma, candidateOptions.corrThresh,false);
else
    candidatePos = findSpotCandidates_MOSAIC(img,candidateOptions.particleRadius,candidateOptions.intensityThreshold,candidateOptions.intensityPtestVar,0,true,false); %0: disable neughbours
end

nrCandidates = size(candidatePos,1);
fitData = cell(nrFrames,1);

if nrCandidates>0
    fitData_temp = psfFit_Image( img, candidatePos.', [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLE, halfw, candidateOptions.sigma );
    if fit_forward
        fitData(1) = {fitData_temp(:,fitData_temp(end,:)==1).'}; %only keep fits with positive exit flag
    else
        fitData(nrFrames) = {fitData_temp(:,fitData_temp(end,:)==1).'};
    end
else
    if calc_once
        fprintf('No particles in first frame, switching to always calc. candidates.\n');
    end
    calc_once = false;
end



% Fit the other frames backwards
msgAccumulator = ''; % Needed for rewindable command line printing (rewPrintf subfunction)
startTime = tic;
elapsedTime = [];
lastElapsedTime = 0;

if(fit_forward)
    for iF = 2:nrFrames %first frame has already been dealt with
        elapsedTime = toc(startTime);
        
        % Output process every 0.5 seconds
        if( (elapsedTime-lastElapsedTime) > 0.5)
            rewindMessages();
            rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iF*(nrFrames-iF)/60),  floor(mod(elapsedTime/iF*(nrFrames-iF),60)))
            rewPrintf('Fitting frame %i/%i\n',iF,nrFrames)
            
            lastElapsedTime = elapsedTime;
        end
        
        if correct_dark
            if useMLE
                img = (double(movie(:,:,iF))+img_dark)*photon_factor;
            else
                img = double(movie(:,:,iF))+img_dark;
            end
        else
            if useMLE
                img = (double(movie(:,:,iF))-bias)*photon_factor;
            else
                img = double(movie(:,:,iF));
            end
        end
        
        if calc_once %in this case, fitted positions of the last frame serve as candidates for this frame
            fitData_temp = psfFit_Image( img, fitData{iF-1}.', [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLE, halfw, candidateOptions.sigma );
            fitData(iF) = {fitData_temp(:,fitData_temp(end,:)==1).'};
        else %otherwise, find new candidates
            if cand_corr
                candidatePos = findSpotCandidates(img, candidateOptions.sigma, candidateOptions.corrThresh,false);
            else
                candidatePos = findSpotCandidates_MOSAIC(img,candidateOptions.particleRadius,candidateOptions.intensityThreshold,candidateOptions.intensityPtestVar,0,~mod(iF-1,bck_interval),false);
            end
            
            nrCandidatesNew = size(candidatePos,1); %to remove empty entries, let's keep track of the largest amount of particles in one frame
            if nrCandidatesNew>0
                fitData_temp = psfFit_Image( img, candidatePos.', [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLE, halfw, candidateOptions.sigma );
                fitData(iF) = {fitData_temp(:,fitData_temp(end,:)==1).'};
            end
        end
    end
else %otherwise, we go backward in time
    for iF = nrFrames-1:-1:1
        elapsedTime = toc(startTime);
        
        % Output process every 0.5 seconds
        if( (elapsedTime-lastElapsedTime) > 0.5)
            rewindMessages();
            rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/(nrFrames-iF)*iF/60),  floor(mod(elapsedTime/(nrFrames-iF)*iF,60)))
            rewPrintf('Fitting frame %i/%i\n',nrFrames-iF+1,nrFrames)
            
            lastElapsedTime = elapsedTime;
        end
        
        if correct_dark
            if useMLE
                img = (double(movie(:,:,iF))+img_dark)*photon_factor;
            else
                img = double(movie(:,:,iF))+img_dark;
            end
        else
            if useMLE
                img = (double(movie(:,:,iF))-bias)*photon_factor;
            else
                img = double(movie(:,:,iF));
            end
        end
        
        if calc_once
            fitData_temp = psfFit_Image( img, fitData{iF+1}.', [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLE, halfw, candidateOptions.sigma );
            fitData(iF) = {fitData_temp(:,fitData_temp(end,:)==1).'};
        else
            if cand_corr
                candidatePos = findSpotCandidates(img, candidateOptions.sigma, candidateOptions.corrThresh,false);
            else
                candidatePos = findSpotCandidates_MOSAIC(img,candidateOptions.particleRadius,candidateOptions.intensityThreshold,candidateOptions.intensityPtestVar,0,~mod(nrFrames-iF,bck_interval),false);
            end
            
            nrCandidatesNew = size(candidatePos,1);
            if nrCandidatesNew>0
                fitData_temp = psfFit_Image( img, candidatePos.', [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLE, halfw, candidateOptions.sigma );
                fitData(iF) = {fitData_temp(:,fitData_temp(end,:)==1).'};
            end
        end
    end
end


rewindMessages();
rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iF*(nrFrames-iF)/60),  floor(mod(elapsedTime/iF*(nrFrames-iF),60)))
rewPrintf('Fitting done.\n');


    function rewPrintf(msg, varargin)
        % Rewindable message printing: Print msg and cache it.
        % Usage is analogous to sprintf.
        msg = sprintf(msg, varargin{:});
        msgAccumulator = [msgAccumulator, msg];
        fprintf(msg);
    end

    function rewindMessages()
        % Remove cached messages from command line, reset cache
        reverseStr = repmat(sprintf('\b'), 1, length(msgAccumulator));
        fprintf(reverseStr);
        
        msgAccumulator = '';
    end

clear global img_bck_mean img_bck_std %get rid of global variables again
end

