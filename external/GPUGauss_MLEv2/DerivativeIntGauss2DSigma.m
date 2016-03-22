% \file DerivativeIntGauss2DSigma.m
% \author Peter Relich
% \date June 02, 2011
% \brief 2D sigma derivative function invoked in the gaussMLEv2 files
% 

%% Derivative of IntGauss with 2DSigma
function [dudt, d2udt2] = DerivativeIntGauss2DSigma( ii, jj, x, y, S, N, PSFx, PSFy )
    % \brief compute the derivate of the 2D gaussian
    %
    % Inputs:
    % \ii index location (used as an off-set) from the external loop
    % \jj index location (used as an off-set) from the external loop
    % \x x-coordinate
    % \y y-coordinate
    % \S Sigma value
    % \N Number of Photons
    % \PSFx Point Spread Function coefficient that remains constant in a specific derivative
    % \PSFy Point Spread Function coefficient that remains constant in a specific derivative
    %
    % Outputs:
    % \dudt   First Derivative
    % \d2udt2 Second Derivative
    
    [dSx, ddSx] = DerivativeIntGauss1DSigma( ii, x, S, N, PSFy );
    [dSy, ddSy] = DerivativeIntGauss1DSigma( jj, y, S, N, PSFx );
    
    dudt = dSx+dSy;
    d2udt2 = ddSx+ddSy;
end