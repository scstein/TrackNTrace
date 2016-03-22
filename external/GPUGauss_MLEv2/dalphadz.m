% \file dalphadz.m
% \author Peter Relich
% \date May 24, 2011
% \brief alpha derivative function invoked in the DerivativeIntGauss2Dz file
% 

%% First derivative of Alpha Function
function d_alph = dalphadz( z, Ax, Bx, d )
    % \brief computes first derivative of alpha in relation to z
    % 
    % Inputs:
    % \z z-coordinate
    % \Ax Constant input for astygmatic imaging
    % \Bx Constant input for astygmatic imaging
    % \d Constant input for astygmatic imaging
    %
    % Output: d_alph, derivative of alpha
    
    d_alph = 2.0*z/d^2 + 3.0*Ax*z.^2/d^3 + 4.0*Bx*z.^3/d^4;
end