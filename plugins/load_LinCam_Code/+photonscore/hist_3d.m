function h = hist_3d(x, x_from, x_to, x_bins, y, varargin)
y_from = x_from;
y_to = x_to;
y_bins = x_bins;
z_from = x_from;
z_to = x_to;
z_bins = x_bins;
h_type = 'double';
if length(varargin) == 1 % z
  z = varargin{1};
elseif length(varargin) == 2 % z, h_type
  z = varargin{1};
  h_type = varargin{2};
elseif length(varargin) == 4 % z, z_from, z_to, z_bins
  z = varargin{1};
  z_from = varargin{2};
  z_to = varargin{3};
  z_bins = varargin{4};
elseif length(varargin) == 5 % z, z_from, z_to, z_bins, h_type
  z = varargin{1};
  z_from = varargin{2};
  z_to = varargin{3};
  z_bins = varargin{4};
  h_type = varargin{5};
elseif length(varargin) == 7 % y_from, y_to, y_bins, z, z_from, z_to, z_bins
  y_from = varargin{1};
  y_to = varargin{2};
  y_bins = varargin{3};
  z = varargin{4};
  z_from = varargin{5};
  z_to = varargin{6};
  z_bins = varargin{7};
elseif length(varargin) == 8 % ...same.. + h_type
  y_from = varargin{1};
  y_to = varargin{2};
  y_bins = varargin{3};
  z = varargin{4};
  z_from = varargin{5};
  z_to = varargin{6};
  z_bins = varargin{7};
  h_type = varargin{8};
else
  error('Inconsistent number of inputs');
end
h = zeros(x_bins, y_bins, z_bins, h_type);
photonscore_mex(photonscore.Const.FN_HIST, ...
    h, x, x_from, x_to, y, y_from, y_to, z, z_from, z_to);

