function h = edges_hist_1d(x, edges, w)
if nargin == 2
  w = ones(size(x));
end

h = photonscore_mex(photonscore.Const.FN_EDGES_HIST_1D, x, edges, w);
end
