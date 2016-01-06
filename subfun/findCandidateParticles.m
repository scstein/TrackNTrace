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
%     candidateData: 1D cell array of xy position estimates (2D row array)
%     used in later fitParticles routine.
% 
%     candidateOptions: see above

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
candidateData = cell(nrFrames,1);

%Call candidate init function
if ~isempty(candidateOptions.initFunc)
    candidateOptions = candidateOptions.initFunc(candidateOptions);
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
    
    candidateData(iLocF) = candidateOptions.perFrameFunc(img,candidateOptions,iLocF);
end

if ~isempty(candidateOptions.postFunc)
    [candidateData,candidateOptions] = candidateOptions.postFunc(candidateData,candidateOptions);
end

rewindMessages();
rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iLocF*(nrFrames-iLocF)/60),  floor(mod(elapsedTime/iLocF*(nrFrames-iLocF),60)))
rewPrintf('Candidate search done.\n');


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

