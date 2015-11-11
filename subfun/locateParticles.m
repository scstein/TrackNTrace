function [ fitData ] = locateParticles( movieStack, imgCorrection, generalOptions, candidateOptions, fittingOptions)
% [ fitData ] = locateParticles( movieStack, imgCorrection, candidateOptions, fittingOptions)
% Find locations of bright spots in an image, in this case a movie of
% fluorescent molecules, and fit a Gaussian distribution to these spots to
% obtain position, amplitude, background level and standard deviation.
%
% INPUT:
%     movieStack: 3D array of intensity values (y,x,N) where N is the
%     number of frames. All images are treated piecewise and only converted
%     to double when needed to avoid memory overflow.
%
%     imgCorrection: 2D array of intensity correction values which is added
%     to each movie image to correct for non-isotropic camera readout.
%     Leave empty [] if no correction is needed.
% 
%     generalOptions: struct of input options for this function, se
%     setDefaultOptions or TrackNTrace manual for details.
%
%     candidateOptions: struct of input options used to find localization
%     candidates. See respective plugin function for details.
%
%     fittingOptions: struct of input options used to fit localization
%     candidates. See respective plugin function for details.
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
fitForward = generalOptions.fitForward;
usePhoton = generalOptions.usePhotonConversion;
if usePhoton
    photon_bias = generalOptions.photonBias;
    photon_factor = generalOptions.photonSensitivity/generalOptions.photonGain;
end

calcOnce = generalOptions.calculateCandidatesOnce;
if ~isempty(imgCorrection) %if correction image is provided, do it
    correctDark = true;
else
    correctDark = false;
end
nrFrames = size(movieStack,3);

candidateFun = candidateOptions.functionHandle;
fittingFun = fittingOptions.functionHandle;


% Trace first frame
if calcOnce %if candidate search only takes place once, an average image is calculated to achieve higher SNR
    if fitForward
        img = mean(double(movieStack(:,:,1:generalOptions.averagingWindowSize)),3);
    else
        img = mean(double(movieStack(:,:,end:-1:end-generalOptions.averagingWindowSize+1)),3);
    end
else
    if fitForward
        img = double(movieStack(:,:,1));
    else
        img = double(movieStack(:,:,end));
    end
end

%convert counts to photons
if usePhoton
    img = (img-photon_bias)*photon_factor;
end

%use dark image stack to correct image for camera artifacts
if correctDark
    if usePhoton
        imgCorrection = (imgCorrection+photon_bias)*photon_factor;
    end
    img = img+imgCorrection;
end

% call candidate function
candidatePos = candidateFun(img,candidateOptions);
nrCandidates = size(candidatePos,1);


fitData = cell(nrFrames,1);
if nrCandidates>0
    fitData_temp = fittingFun(img,candidatePos,fittingOptions);
    if fitForward
        fitData(1) = fitData_temp;
    else
        fitData(nrFrames) = fitData_temp;
    end
else
    if calcOnce
        fprintf('No particles in first frame, switching to always calc. candidates.\n');
    end
    calcOnce = false;
end



% Fit the other frames backwards
msgAccumulator = ''; % Needed for rewindable command line printing (rewPrintf subfunction)
startTime = tic;
elapsedTime = [];
lastElapsedTime = 0;

global iLocF %loop variable can be known by any function

if(fitForward)
    for iLocF = 2:nrFrames %first frame has already been dealt with
        elapsedTime = toc(startTime);
        
        % Output process every 0.5 seconds
        if( (elapsedTime-lastElapsedTime) > 0.5)
            rewindMessages();
            rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iLocF*(nrFrames-iLocF)/60),  floor(mod(elapsedTime/iLocF*(nrFrames-iLocF),60)))
            rewPrintf('Fitting frame %i/%i\n',iLocF,nrFrames)
            
            lastElapsedTime = elapsedTime;
        end
        
        if correctDark
            if usePhoton
                img = (double(movieStack(:,:,iLocF))-photon_bias)*photon_factor+imgCorrection;
            else
                img = double(movieStack(:,:,iLocF))+imgCorrection;
            end
        else
            if usePhoton
                img = (double(movieStack(:,:,iLocF))-photon_bias)*photon_factor;
            else
                img = double(movieStack(:,:,iLocF));
            end
        end
        
        if calcOnce %in this case, fitted positions of the last frame serve as candidates for this frame
            fitData(iLocF) = fittingFun(img,fitData{iLocF-1},fittingOptions);
        else %otherwise, find new candidates
            candidatePos = candidateFun(img,candidateOptions);
            nrCandidatesNew = size(candidatePos,1); %to remove empty entries, let's keep track of the largest amount of particles in one frame
            if nrCandidatesNew>0
                fitData(iLocF) = fittingFun(img,candidatePos,fittingOptions);
            end
        end
    end
else %otherwise, we go backward in time
    for iLocF = nrFrames-1:-1:1
        elapsedTime = toc(startTime);
        
        % Output process every 0.5 seconds
        if( (elapsedTime-lastElapsedTime) > 0.5)
            rewindMessages();
            rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/(nrFrames-iLocF)*iLocF/60),  floor(mod(elapsedTime/(nrFrames-iLocF)*iLocF,60)))
            rewPrintf('Fitting frame %i/%i\n',nrFrames-iLocF+1,nrFrames)
            
            lastElapsedTime = elapsedTime;
        end
        
        if correctDark
            if usePhoton
                img = (double(movieStack(:,:,iLocF))-photon_bias)*photon_factor+imgCorrection;
            else
                img = double(movieStack(:,:,iLocF))+imgCorrection;
            end
        else
            if usePhoton
                img = (double(movieStack(:,:,iLocF))-photon_bias)*photon_factor;
            else
                img = double(movieStack(:,:,iLocF));
            end
        end
        
        if calcOnce
            fitData(iLocF) = fittingFun(img,fitData{iLocF+1},fittingOptions);
        else
            candidatePos = candidateFun(img,candidateOptions);
            nrCandidatesNew = size(candidatePos,1);
            if nrCandidatesNew>0
                fitData(iLocF) = fittingFun(img,candidatePos,fittingOptions);
            end
        end
    end
end


rewindMessages();
rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iLocF*(nrFrames-iLocF)/60),  floor(mod(elapsedTime/iLocF*(nrFrames-iLocF),60)))
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

clear global iLocF
end

