function [medi, mea] = medimean(fl, range)
if nargin == 1
  range = [0 65535];
end

[medi, mea] = photonscore_mex(photonscore.Const.FN_FLIM_MEDIMEAN, ...
    fl.image, fl.time, uint16(range));
end
