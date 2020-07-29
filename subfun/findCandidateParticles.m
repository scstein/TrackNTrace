function [ candidateData, candidateOptions ] = findCandidateParticles( movieStack, darkImage, globalOptions, candidateOptions)
% [ candidateData, candidateOptions ] = findCandidateParticles( movieStack, darkImage, globalOptions, candidateOptions)
% Find rough estimate of locations of bright spots in an image, in this
% case a movie of fluorescent molecules, for later refinement.
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
%     candidateOptions: struct of input options used to find localization
%     candidates. See respective plugin function for details.
%
%
% OUTPUT:
%     candidateData: Nx1 cell array of fitted positions where 
%        N is the number of analyzed frames. Each cell contains a 
%        KxP double array, where P is the number of model parameters
%        and K is the maximum amount of particles in the respective frame. 
%        Each row represents a unique candidate fit and the column order is 
%        [x, y]. These columns are mandatory for the refinement to work. 
%        Plugins can output extra data.
%     used in later fitParticles routine.
%
%     candidateOptions: see above

global imgCorrection;
global parallelProcessingAvailable

% Parse inputs
if ~isempty(darkImage) %if correction image is provided, do it
    imgCorrection = darkImage;
    
    %if image is converted to photons, correction image has to be converted too
    %but the bias is already included. Add the bias again to avoid later
    %confusion
    if globalOptions.usePhotonConversion
        imgCorrection = (imgCorrection+double(globalOptions.binFrame)*globalOptions.photonBias)*(globalOptions.photonSensitivity/globalOptions.photonGain);
    end
end

nrFrames = size(movieStack,3);
candidateData = cell(nrFrames,1);

totalTime_start = tic;

tic;
%Call plugins init function
if ~isempty(candidateOptions.initFunc)
    candidateOptions = candidateOptions.initFunc(candidateOptions);
end
initTime = toc;

% Try parallel processing of plugins main function
parallelProcessing_Failed = false;
if parallelProcessingAvailable && candidateOptions.useParallelProcessing
    try       
        imgCorrectionLocal = imgCorrection; % Need a copy for parallel processing
        
        fprintf('TNT: Locating candidates using ''%s''(parallel processing).\n', candidateOptions.plugin_name);
        mainTime_start = tic;
        parfor iFrame = 1:nrFrames
            img = correctMovie_Parallel(movieStack(:,:,iFrame), globalOptions, imgCorrectionLocal);
            candidateData(iFrame) = {candidateOptions.mainFunc(img,candidateOptions,iFrame)};
        end
        mainTime = toc(mainTime_start);
    catch err
        warning off backtrace
        warning('Parallel execution failed. Switching to serial execution.\n Error: %s.',err.message);
        warning on backtrace
        parallelProcessing_Failed = true;
    end
end

% Standard serial processing of plugins main function
if not(parallelProcessingAvailable) || not(candidateOptions.useParallelProcessing) || parallelProcessing_Failed
    if parallelProcessingAvailable && not(candidateOptions.useParallelProcessing)
        fprintf('TNT: Locating candidates using ''%s'' (parallel processing disabled by plugin).\n', candidateOptions.plugin_name);
    else
        fprintf('TNT: Locating candidates.\n');
    end    
    
    msgAccumulator = ''; % Needed for rewindable command line printing (rewPrintf subfunction)
    mainTime_start = tic;
    mainTime = [];
    lastElapsedTime = 0;
    
    for iFrame = 1:nrFrames %first frame has already been dealt with
        mainTime = toc(mainTime_start);
        
        % Output process every 0.5 seconds
        if( (mainTime-lastElapsedTime) > 0.5)
            rewindMessages();
            rewPrintf('TNT: Time elapsed %im %is - to go: %im %is\n', floor(mainTime/60), floor(mod(mainTime,60)),  floor(mainTime/iFrame*(nrFrames-iFrame)/60),  floor(mod(mainTime/iFrame*(nrFrames-iFrame),60)))
            rewPrintf('TNT: Locating candidates in frame %i/%i\n',iFrame,nrFrames)
            
            lastElapsedTime = mainTime;
        end
        
        img = correctMovie(movieStack(:,:,iFrame));
        
        candidateData(iFrame) = {candidateOptions.mainFunc(img,candidateOptions,iFrame)};
    end
    rewindMessages();
    rewPrintf('TNT: Time elapsed %im %is - to go: %im %is\n', floor(mainTime/60), floor(mod(mainTime,60)),  floor(mainTime/iFrame*(nrFrames-iFrame)/60),  floor(mod(mainTime/iFrame*(nrFrames-iFrame),60)))
end

% Call plugins post-processing function
tic;
if ~isempty(candidateOptions.postFunc)
    [candidateData,candidateOptions] = candidateOptions.postFunc(candidateData,candidateOptions);
end
postTime = toc;


totalTime = toc(totalTime_start);
fprintf('TNT: Candidate search took %im %is (init: %im %is, main: %im %is, post: %im %is).\n', floor(totalTime/60), floor(mod(totalTime,60)),floor(initTime/60), floor(mod(initTime,60)) ,floor(mainTime/60), floor(mod(mainTime,60)) ,floor(postTime/60), floor(mod(postTime,60)));

% Verify the outParamDescription, make it fit to the data if neccessary
candidateOptions = verifyOutParamDescription(candidateData, candidateOptions);


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

