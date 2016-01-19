% \file d2alphadz2.m
% \author Peter Relich
% \date May 24, 2011
% \brief second alpha derivative function invoked in the DerivativeIntGauss2Dz file
% 

%% Second derivative of Alpha function
function d2_alph = d2alphadz2( z, Ax, Bx, d )
    % \brief computes second derivative of alpha in relation to z
    %
    % Inputs:
    % \z z-coordinate
    % \Ax Constant input for astygmatic imaging
    % \Bx Constant input for astygmatic imaging
    % \d Constant input for astygmatic imaging
    %
    % Output: d2_alpha, second derivative of alpha
    
    d2_alph = 2.0/d^2 + 6.0*Ax*z/d^3 + 12.0*Bx*z.^2/d^4;
end