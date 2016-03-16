function [ fittingData, fittingOptions ] = fitParticles( movieStack, darkImage, globalOptions, fittingOptions, candidateData)
% [ fittingData ] = locateParticles( movieStack, darkImage, globalOptions, fittingOptions, candidateData)
% Find locations of bright spots in an image, in this case a movie of
% fluorescent molecules, and fit a Gaussian distribution to these spots to
% obtain position, amplitude, background level and standard deviation.
%
% INPUT:
%     movieStack: 3D array of intensity values (y,x,N) where N is the
%     number of frames. All images are treated piecewise and only converted
%     to double when needed to avoid memory overflow.
%
%     darkImage: 2D array of intensity correction values which is added
%     to each movie image to correct for non-isotropic camera readout.
%     Leave empty [] if no correction is needed.
%
%     globalOptions: struct of input options for this function, se
%     setDefaultOptions or TrackNTrace manual for details.
%
%     fittingOptions: struct of input options used to fit localization
%     candidates. See respective plugin function for details.
%
%     candidateData: Nx1 cell array of candidate positions where 
%        N is the number of analyzed frames. Each cell contains a 
%        KxP double array, where P is the number of model parameters
%        and K is the maximum amount of particles in the respective frame. 
%        Each row represents a unique candidate fit and the column order is 
%        [x, y]. These columns are mandatory for the fitting to work. 
%        Plugins can output extra data.
%
%
% OUTPUT:
%     fittingData: Nx1 cell array of fitted positions where N is
%        the number of analyzed frames. Each cell contains a KxP
%       double array, where P is the number of model parameters and K
%        is the maximum amount of particles in the respective frame. 
%        Each row represents a unique fit and the column order according to 
%        the model PSF is [\mu_x, \mu_y, \mu_z, Amplitude, Background].
%        These five columns are mandatory for most trackers to work correctly
%        in the next step. Plugins can output extra data.
%
%     fittingOptions: see above

global imgCorrection;
global parallelProcessingAvailable

% Parse inputs
if ~isempty(darkImage) %if correction image is provided, do it
    imgCorrection = darkImage;
    
    %if image is converted to photons, correction image has to be converted too
    %but the bias is already included. Add the bias again to avoid later
    %confusion
    if globalOptions.usePhotonConversion
        imgCorrection = (imgCorrection+globalOptions.photonBias)*(globalOptions.photonSensitivity/globalOptions.photonGain);
    end
end

nrFrames = size(movieStack,3);
fittingData = cell(nrFrames,1);

%Call plugins init function
if ~isempty(fittingOptions.initFunc)
    fittingOptions = fittingOptions.initFunc(fittingOptions);
end


% Try parallel processing of plugins main function
parallelProcessing_Failed = false;
if parallelProcessingAvailable && fittingOptions.useParallelProcessing
    try      
        imgCorrectionLocal = imgCorrection; % Need a copy for parallel processing
        
        fprintf('TNT: Fitting candidates using parallel processing (Frame by Frame).\n');
        startTime = tic;        
        parfor iFrame = 1:nrFrames
            img = correctMovie_Parallel(movieStack(:,:,iFrame), globalOptions, imgCorrectionLocal);
            nrCandidates = size(candidateData{iFrame},1);
            if nrCandidates>0
                fittingData(iFrame) = {fittingOptions.mainFunc(img,candidateData{iFrame},fittingOptions,iFrame)};
            end
        end
        totalTime = toc(startTime);
        fprintf('TNT: Time elapsed %im %is.\n',floor(totalTime/60), floor(mod(totalTime,60)));
    catch err
        warning off backtrace
        warning('Parallel execution failed. Switching to serial execution.\n Error: %s.',err.message);
        warning on backtrace
        parallelProcessing_Failed = true;
    end
end

% Standard serial processing of plugins main function
if not(parallelProcessingAvailable) || not(fittingOptions.useParallelProcessing) || parallelProcessing_Failed
    if parallelProcessingAvailable && not(fittingOptions.useParallelProcessing)
        fprintf('TNT: Fitting candidates (parallel processing disabled by plugin).\n');
    else
        fprintf('TNT: Fitting candidates.\n');
    end
    
    msgAccumulator = ''; % Needed for rewindable command line printing (rewPrintf subfunction)
    startTime = tic;
    elapsedTime = [];
    lastElapsedTime = 0;
    
    
    for iFrame = 1:nrFrames %first frame has already been dealt with
        elapsedTime = toc(startTime);
        
        % Output process every 0.5 seconds
        if( (elapsedTime-lastElapsedTime) > 0.5)
            rewindMessages();
            rewPrintf('TNT: Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iFrame*(nrFrames-iFrame)/60),  floor(mod(elapsedTime/iFrame*(nrFrames-iFrame),60)))
            rewPrintf('TNT: Fitting frame %i/%i\n',iFrame,nrFrames)
            
            lastElapsedTime = elapsedTime;
        end
        
        img = correctMovie(movieStack(:,:,iFrame));
        
        nrCandidates = size(candidateData{iFrame},1);
        if nrCandidates>0
            fittingData(iFrame) = {fittingOptions.mainFunc(img,candidateData{iFrame},fittingOptions,iFrame)};
        end
    end
    rewindMessages();
    rewPrintf('TNT: Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iFrame*(nrFrames-iFrame)/60),  floor(mod(elapsedTime/iFrame*(nrFrames-iFrame),60)))
end


% Call plugins post-processing function
if ~isempty(fittingOptions.postFunc)
    [fittingData,fittingOptions] = fittingOptions.postFunc(fittingData,fittingOptions);
end

fprintf('TNT: Fitting done.\n');


% Verify the outParamDescription, make it fit to the data if neccessary
fittingOptions = verifyOutParamDescription(fittingData, fittingOptions);

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

end

