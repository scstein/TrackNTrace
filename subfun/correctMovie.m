function [correctedStack] = correctMovie(movieStack)
% Correct movie frame with dark image and/or convert image counts to
% photons.
% 
% INPUT:
%     movieStack: 3D array of intensity values (y,x,N) where N is the
%     number of frames. All images are treated piecewise and only converted
%     to double when needed to avoid memory overflow.    
%     
% OUTPUT:
%     correctedStack: 3D array of corrected intensity values.


global globalOptions
global imgCorrection

% from globalOptions, imgCorrection:
% imgCorrection: 2D correction image created from a movie of dark
% images taken with closed shutter. See calculateDark.m or manual for
% details. The variable referred to here is initialized in
% locateParticles.m and changed form the original if photon conversion is
% enabled. See locateParticles.m for details.
% 
% photonBias: Image count bias of camera.
% 
% photonFactor: Sensitivity divided by gain, see camera specification
% for details.



correctDark = true;
if nargin < 2 || isempty(imgCorrection)
    correctDark = false;
end

usePhoton = globalOptions.usePhotonConversion;

nImages = size(movieStack,3);


if correctDark
    if usePhoton
        correctedStack = (double(movieStack)-globalOptions.photonBias)*globalOptions.photonFactor+repmat(imgCorrection,[1,1,nImages]);
    else
        correctedStack = double(movieStack)+repmat(imgCorrection,[1,1,nImages]);
    end
else
    if usePhoton
        correctedStack = (double(movieStack)-globalOptions.photonBias)*globalOptions.photonFactor;
    else
        correctedStack = double(movieStack);
    end
end

end