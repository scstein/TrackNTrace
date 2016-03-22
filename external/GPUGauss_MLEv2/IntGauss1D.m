% \file IntGauss1D.m
% \author Peter Relich
% \date June 02, 2011
% \brief 1-D integrated Gaussian function invoked in the gaussMLEv2 files
% 

%% IntGauss1D function, returns gauss function of input parameters
function AvInt = IntGauss1D( ii, x, sigma )
    % \brief computes average of off-centered error function integrals
    % 
    % Inputs:
    % \ii counting parameter, sets off-set in error function
    % \x position parameter in error function
    % \sigma sigma value of the PSF
    % 
    % Output: avInt, the average error function integral value
    
    norm = 0.5 ./ (sigma.^2);
    AvInt = 0.5 * ( erf( (ii - x + 0.5).*sqrt(norm) ) - erf( (ii - x -0.5).*sqrt(norm) ) );
end