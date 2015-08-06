function [ fitData ] = fitMovie( movie, candidateOptions, fittingOptions)
% - Fit a series of candidates in a given movie frame by frame.
% - Candidates are extracted from an average of "windowSize" frames.
% - Movie can be fitted forwards or backwards. In the second case the
% candidates are determined from the last "windowsSize" frames and the fit
% propagates initial conditions backwards through the movie.
% 
% Input:
%     candidateOptions.sigma
%     candidateOptions.corrThresh
%     candidateOptions.windowSize
%     candidateOptions.fitForward
% 
%     fittingOptions.fitSigma
%     fittingOptions.usePixelIntegratedFit
%     fittingOptions.useMLErefine
%
% Output:
%     fitData: Matrix 6xNxF. Fitted values for the 6 parameters (xpos, ypos,
%              A, BG, sigma) for N candidates in F frames.
%
% HINT: Run this function once or twice on a very small movie to let MATLAB
% perform internal optimizations. The run on the big data might be
% significiantly smaller after that..

% Input parsing
fit_forward = candidateOptions.fitForward; % Fit in forward direction
fit_sigma = fittingOptions.fitSigma;
usePixelIntegratedFit = fittingOptions.usePixelIntegratedFit;
                        fittingOptions.usePixelIntegratedFit
useMLErefine = fittingOptions.useMLErefine;

% Take average of windowSize frames to select candidates
if(fit_forward)
    avg = mean(movie(:,:,1:candidateOptions.windowSize),3);
else
    avg = mean(movie(:,:,end-candidateOptions.windowSize:end),3);
end

% Compute candidates = [col,row,Amp]
candidatePos = findSpotCandidates(avg, candidateOptions.sigma, candidateOptions.corrThresh,false);
nrCandidates = size(candidatePos,1);
fprintf('Identified %i candidates for fitting.\n', nrCandidates);

% Compute drift by fitting
nrFrames = size(movie,3);
fitData = zeros(6,nrCandidates,nrFrames); %to save [col,row,Amp,BG] per Candidate per frame

% Start with first fit or last frame
halfw = round(4*candidateOptions.sigma);
if(fit_forward)
    fitData(:,:,1) = psfFit_Image( movie(:,:,1), candidatePos', [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLErefine, halfw, candidateOptions.sigma );
else
    fitData(:,:,nrFrames) = psfFit_Image( movie(:,:,end), candidatePos', [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLErefine, halfw, candidateOptions.sigma );
end


% Fit the other frames backwards
msgAccumulator = ''; % Needed for rewindable command line printing (rewPrintf subfunction)
startTime = tic;
elapsedTime = [];
lastElapsedTime = 0;

if(fit_forward)
    for iF = 2:nrFrames
        elapsedTime = toc(startTime);
        
        % Output process every 0.5 seconds
        if( (elapsedTime-lastElapsedTime) > 0.5)
            rewindMessages();
            rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iF*(nrFrames-iF)/60),  floor(mod(elapsedTime/iF*(nrFrames-iF),60)))
            rewPrintf('Fitting frame %i/%i\n',iF,nrFrames)
            
            lastElapsedTime = elapsedTime;
        end
        
        fitData(:,:,iF) = psfFit_Image( movie(:,:,iF), fitData(:,:,iF-1), [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLErefine, halfw, candidateOptions.sigma );
    end
else
    for iF = nrFrames-1:-1:1
        elapsedTime = toc(startTime);
        
        % Output process every 0.5 seconds
        if( (elapsedTime-lastElapsedTime) > 0.5)
            rewindMessages();
            rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/(nrFrames-iF)*iF/60),  floor(mod(elapsedTime/(nrFrames-iF)*iF,60)))
            rewPrintf('Fitting frame %i/%i\n',nrFrames-iF+1,nrFrames)
            
            lastElapsedTime = elapsedTime;
        end
        
        fitData(:,:,iF) = psfFit_Image( movie(:,:,iF), fitData(:,:,iF+1), [1,1,1,1,fit_sigma], usePixelIntegratedFit, useMLErefine, halfw, candidateOptions.sigma );
    end
end


rewindMessages();
rewPrintf('Time elapsed %im %is - to go: %im %is\n', floor(elapsedTime/60), floor(mod(elapsedTime,60)),  floor(elapsedTime/iF*(nrFrames-iF)/60),  floor(mod(elapsedTime/iF*(nrFrames-iF),60)))
rewPrintf('Fitting frame %i/%i',iF,nrFrames)
rewPrintf(' .. done\n')


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
        fprintf([reverseStr]);
        
        msgAccumulator = '';
    end

end

