% \file DerivativeIntGauss1DSigma.m
% \author Peter Relich
% \date June 02, 2011
% \brief Sigma Derivative function invoked in the gaussMLEv2 files
% 

%% Derivative of IntGauss with 1DSigma
function [dudt, d2udt2] = DerivativeIntGauss1DSigma( ii, x, Sx, N, PSFy )
    % \brief Compute derivative of the 1D gaussian
    %
    % Inputs:
    % \ii index location (used as an off-set) from the external loop
    % \x x-coordinate
    % \Sx Sigma in the x direction
    % \N Number of Photons
    % \PSFy Point Spread Function coefficient that remains constant in this derivative
    % 
    % Outputs:
    % \dudt first derivative
    % \d2udt2 second derivative
    
    ax = exp(-0.5*( (ii + 0.5 - x)./Sx ).^2);
    bx = exp(-0.5*( (ii - 0.5 - x)./Sx ).^2);
    
    dudt = -N.*(ax.*(ii + 0.5 - x) - bx.*(ii - 0.5 - x)).*PSFy ./ (sqrt(2.0*pi).*Sx.^2);
    
    d2udt2 = -2*dudt./Sx - N.*(ax.*(ii + 0.5 - x).^3 - bx.*(ii - 0.5 - x).^3)...
        .*PSFy ./ (sqrt(2.0*pi).*Sx.^5);
end