function m = gaussian_irf(n, mu, fwhm)
% GAUSSIAN_IRF simulate gaussian IRF
%   m = p.gaussian_irf(m, mu, fwhm)
%
%   N number of output channels
%   MU means
%   FWHM widths of the gaussian at half maximum
m = photonscore_mex(photonscore.Const.FN_GAUSSIAN_IRF, n, mu, fwhm);
