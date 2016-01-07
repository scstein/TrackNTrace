function [ fitData, fittingOptions ] = fitParticles( movieStack, darkImage, globalOptions, fittingOptions, candidateData)
% [ fitData ] = locateParticles( movieStack, darkImage, globalOptions, fittingOptions, candidateData)
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
%     candidateData: Cell array of position estimates used in fitting routine
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
% 
%     fittingOptions: see above

global imgCorrection;

% Parse inputs
if ~isempty(darkImage) %if correction image is provided, do it
    imgCorrection = darkImage;
    
    %if image is converted to photons, correction image has to be converted too
    %but the bias is already included. Add the bias again to avoid later
    %confusion
    if globalOptions.usePhotonConversion
        imgCorrection = (imgCorrection+globalOptions.photonBias)*globalOptions.photonFactor;
    end
end

nrFrames = size(movieStack,3);
fitData = cell(nrFrames,1);

%Call fit init function
if ~isempty(fittingOptions.initFunc)
    fittingOptions = fittingOptions.initFunc(fittingOptions);
end


msgAccumulator = ''; % Needed for rewindable command line printing (rewPrintf subfunction)
startTime = tic;
elapsedTime = [];
lastElapsedTime = 0;

for iLocF = 1:nrFrames %first frame has already been dealt with
    elapsedTime = toc(startTime);
    
    % Output process every 0.5 seconds
    if( (elapsedTime-lastElapsedTime) > 0.5)
        rewindMessages();
        rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iLocF*(nrFrames-iLocF)/60),  floor(mod(elapsedTime/iLocF*(nrFrames-iLocF),60)))
        rewPrintf('Fitting frame %i/%i\n',iLocF,nrFrames)
        
        lastElapsedTime = elapsedTime;
    end
    
    img = correctMovie(movieStack(:,:,iLocF));
    
    nrCandidates = size(candidateData{iLocF},1);
    if nrCandidates>0
        fitData(iLocF) = {fittingOptions.mainFunc(img,candidateData{iLocF},fittingOptions,iLocF)};
    end
end

if ~isempty(fittingOptions.postFunc)
    [fitData,fittingOptions] = fittingOptions.postFunc(fitData,fittingOptions);
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

end

