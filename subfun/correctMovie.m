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

% Here we simply call the parallel version (which is not allowed to use
% global variables) to avoid code duplication. Overhead of one additional
% function call is negligible.
correctedStack = correctMovie_Parallel(movieStack,globalOptions,imgCorrection);

end