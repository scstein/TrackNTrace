% \file alpha.m
% \author Peter Relich
% \date May 24, 2011
% \brief alpha function invoked in the DerivativeIntGauss2Dz file
% 

%% Alpha function
function alph = alpha( z, Ax, Bx, d )
    % \brief computes coefficient for alpha, is this coefficient used?
    %
    % Inputs:
    % \z z-coordinate
    % \Ax Constant input for astygmatic imaging
    % \Bx Constant input for astygmatic imaging
    % \d Constant input for astygmatic imaging
    %
    % Ouput: alpha
    
    alph = 1 + (z/d).^2 + Ax*(z/d).^3 + Bx*(z/d).^4;
end
