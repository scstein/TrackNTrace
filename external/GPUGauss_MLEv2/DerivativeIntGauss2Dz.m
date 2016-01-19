% \file DerivativeIntGauss2Dz.m
% \author Peter Relich
% \date May 24, 2011
% \brief z derivative of Gaussian function invoked in the gaussMLEv2 files
% 

%% Derivative of IntGauss with wrt. Z in the
function [pPSFx, pPSFy, dudt, d2udt2] = DerivativeIntGauss2Dz( ii, jj, theta, PSFSigma_x, PSFSigma_y, Ax, Ay, Bx, By, gamma, d )
    % \brief Compute derivatives of the 2D gaussian w.r.t z
    %
    % Inputs:
    % \ii index location (used as an off-set) from the external loop
    % \jj index location (used as an off-set) from the external loop
    % \theta array of coordinates utilized [x,y,N,bg,z]
    % \PSFSigma_x Point Spread Function Sigma in the x direction
    % \PSFSigma_y Point Spread Function Sigma in the y direction
    % \Ax Constant input for astygmatic imaging, x bias
    % \Ay Constant input for astygmatic imaging, y bias
    % \Bx Constant input for astygmatic imaging, x bias
    % \By Constant input for astygmatic imaging, y bias
    % \gamma Constant input for astygmatic imaging
    % \d Constant input for astygmatic imaging
    % 
    % Outputs:
    % \pPSFx Point Spread Function in the x direction
    % \pPSFy Point Spread Function in the y direction
    % \dudt first derivative in z
    % \d2udt2 second derivative in z
    
    z = theta(5,:);
    
    alphax = alpha( z-gamma, Ax, Bx, d );
    alphay = alpha( z+gamma, Ay, By, d );
    
    Sx = PSFSigma_x*sqrt(alphax);
    Sy = PSFSigma_y*sqrt(alphay);
    
    PSFx = IntGauss1D( ii, theta(1,:), Sx );
    PSFy = IntGauss1D( jj, theta(2,:), Sy );
    pPSFx = PSFx;
    pPSFy = PSFy;
    
    [dudt(1,:), ddx] = DerivativeIntGauss1D( ii, theta(1,:), Sx, theta(3,:), PSFy );
    [dudt(2,:), ddy] = DerivativeIntGauss1D( jj, theta(2,:), Sy, theta(3,:), PSFx );
    [dSx, ddSx] = DerivativeIntGauss1DSigma( ii, theta(1,:), Sx, theta(3,:), PSFy );
    [dSy, ddSy] = DerivativeIntGauss1DSigma( jj, theta(2,:), Sy, theta(3,:), PSFx );
    
    dSdalpha_x = PSFSigma_x ./ (2.0*sqrt(alphax));
    dSdalpha_y = PSFSigma_y ./ (2.0*sqrt(alphay));
    dSdzx = dSdalpha_x.*dalphadz( z-gamma, Ax, Bx, d );
    dSdzy = dSdalpha_y.*dalphadz( z+gamma, Ay, By, d );
    dudt(5,:) = dSx.*dSdzx + dSy.*dSdzy;
    
    d2udt2(1,:) = ddx;
    d2udt2(2,:) = ddy;
    
    d2Sdalpha2_x = -PSFSigma_x ./ (4.0*(alphax).^1.5);
    d2Sdalpha2_y = -PSFSigma_y ./ (4.0*(alphay).^1.5);
    
    ddSddzx = d2Sdalpha2_x.*(dalphadz( z-gamma, Ax, Bx, d )).^2 + dSdalpha_x.*d2alphadz2( z-gamma, Ax, Bx, d );
    ddSddzy = d2Sdalpha2_y.*(dalphadz( z+gamma, Ay, By, d )).^2 + dSdalpha_y.*d2alphadz2( z+gamma, Ay, By, d );
    
    d2udt2(5,:) = ddSx.*(dSdzx.^2) + dSx.*ddSddzx + ddSy.*(dSdzy.^2) + dSy.*ddSddzy;

end