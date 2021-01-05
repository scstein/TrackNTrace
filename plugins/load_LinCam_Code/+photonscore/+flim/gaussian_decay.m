function m = gaussian_decay(n, mu, fwhm, tau)
% GAUSSIAN_DECAY simulate gaussian-exp convolution
%   m = p.gaussian_decay(m, mu, fwhm, tau)
%
%   N number of output channels
%   MU mean
%   FWHM width of the gaussian at half maximum
%   TAU lifetimes to simulate

if numel(mu) == 1 && numel(fwhm) == 1 && numel(tau) > 1
    n = repmat(n, size(tau));
    mu = repmat(mu, size(tau));
    fwhm = repmat(fwhm, size(tau));
end

m = photonscore_mex(photonscore.Const.FN_GAUSSIAN_DECAY, ...
    n, mu, fwhm, tau);

