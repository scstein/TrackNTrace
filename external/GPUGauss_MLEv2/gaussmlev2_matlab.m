%gaussmlev2_matlab  matlab implementation of gaussmlev2
%
%   see gaussmlev2.m for help and calling convention

% \file GaussMLEv2.m
% \author Peter Relich
% \date May 31, 2011

function [P, CRLB, LL] = gaussmlev2_matlab(data, PSFSigma, iterations, fittype, Ax, Ay, Bx, By, gamma, d)

    switch fittype
        case 1
            [P, CRLB, LL]=GaussMLEv2(data,PSFSigma,iterations);
        case 2
            [P, CRLB, LL]=GaussMLEv2_sigma(data,PSFSigma,iterations);
            return;
        case 3
            [P, CRLB, LL]=GaussMLEv2_z(data,PSFSigma,iterations,Ax,Ay,Bx,By,gamma,d);
            return;
        case 4
            [P, CRLB, LL]=GaussMLEv2_sigmaxy(data,PSFSigma,iterations);
            return;
    end
    
end


