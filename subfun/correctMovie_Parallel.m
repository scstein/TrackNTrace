function [correctedStack] = correctMovie_Parallel(movieStack, globalOptions, imgCorrection)
% Correct movie frame with dark image and/or convert image counts to
% photons.
%
% Photon conversion is performed using
%    I_phot =(I_ADC-bias)*sensitivity/gain.
% where these parameters can be found in the camera specification sheet.
%
% The _Parallel version of this function can be called in parfor/spmd statements 
% and can thus not use global variables.
% 
% INPUT:
%     movieStack: 3D array of intensity values (y,x,N) where N is the
%     number of frames. All images are treated piecewise and only converted
%     to double when needed to avoid memory overflow.
%
%     globalOptions: struct housing the following parameters
%        - photonBias: Image count bias of camera.% 
%        - photonSensitivity: Sensitivity of camera, see camera specification for details.
%        - photonGain: Gain of camera during recording.
%        - binFrame: number of raw frames summed to one frame in movieStack
%
%     imgCorrection: 2D correction image created from a movie of dark
%     images taken with closed shutter. See helper/CalculateDark.m
%     
% OUTPUT:
%     correctedStack: 3D array of corrected intensity values.

correctDark = true;
if nargin < 2 || isempty(imgCorrection)
    correctDark = false;
end

usePhoton = globalOptions.usePhotonConversion;

nImages = size(movieStack,3);

if usePhoton && isinf(globalOptions.binFrame)
    if globalOptions.photonBias == 0
        % prevents errors (nan output) when Bias is 0 and Binning is inf.
        globalOptions.binFrame = 1;
    else 
        % For Bias>0, we would need to know how many frames are acutually binned
        % -> Acces to importPlugin/write to metadata. Warn for now
        globalOptions.binFrame = 1;
        warning('Error correcting the bias of the movie. Infinite frame binning not compatible with bias correction. Assuming frame binning 1.');
    end
end

if correctDark
    if usePhoton
        correctedStack = (double(movieStack)-double(globalOptions.binFrame)*globalOptions.photonBias)*(globalOptions.photonSensitivity/globalOptions.photonGain)+repmat(imgCorrection,[1,1,nImages]);
    else
        correctedStack = double(movieStack)+repmat(imgCorrection,[1,1,nImages]);
    end
else
    if usePhoton
        correctedStack = (double(movieStack)-double(globalOptions.binFrame)*globalOptions.photonBias)*(globalOptions.photonSensitivity/globalOptions.photonGain);
    else
        correctedStack = double(movieStack);
    end
end

end