function h = hist_1d(x, from, to, bins, varargin)
h_type = 'double';
if ~isempty(varargin)
  h_type = varargin{1};
end
h = zeros(bins, 1, h_type);

photonscore_mex(photonscore.Const.FN_HIST, h, x, from, to)
