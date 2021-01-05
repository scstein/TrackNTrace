function m = convolve(irf, irf_shif, channels, tau, tau_ref)
% EXP_CONVOLVE Compute exponential decay convolution integral with IRF.
%   m = p.exp_convolve(irf, irf_shift, n, tau) convole
%
%   IRF vector of instrument respone function
%   IRF_SHIFT offset of IRF
%   CHANNELS number of output channels
%   TAU vector or scalar of lifetime
%   TAU_REF reference lifetime
switch nargin
case 4
  m = photonscore_mex(photonscore.Const.FN_FLIM_CONV, ...
      irf, irf_shif, channels, tau, -1);
case 5
  m = photonscore_mex(photonscore.Const.FN_FLIM_CONV, ...
      irf, irf_shif, channels, tau, tau_ref);
otherwise
  error('Wrong number of arguments');
end
