% \file DerivativeIntGauss1D.m
% \author Peter Relich
% \date June 02, 2011
% \brief Derivative Gaussian function invoked in the gaussMLEv2 files
% 

%% Derivative of IntGauss1D
function [dudt, d2udt2] = DerivativeIntGauss1D( ii, x, sigma, N, PSFy )
    % \brief computes the 1st and 2nd derivatives of the 1D gaussian
    % 
    % Inputs:
    % \ii index location (used as an off-set) from the external loop
    % \x x-coordinate
    % \sigma sigma
    % \N number of Photons
    % \PSFy Point Spread Function coefficient that remains constant in this derivative
    % 
    % Outputs:
    % \Idudt first derivative
    % \Id2udt2 second derivative
    
    a = exp(-0.5*( (ii + 0.5 - x)./sigma ).^2);
    b = exp(-0.5*( (ii - 0.5 - x)./sigma ).^2);
    
    dudt = -N.*PSFy.*(a - b) ./ (sqrt(2.0*pi).*sigma);
    d2udt2 = -N.*( (ii + 0.5 - x).*a - (ii - 0.5 - x).*b ).*PSFy ./ ( sqrt(2.0*pi).*sigma.^3 );
end