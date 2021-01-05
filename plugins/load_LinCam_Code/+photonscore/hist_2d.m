function h = hist_2d(x, x_from, x_to, x_bins, y, varargin)
if length(varargin) <= 4
  y_from = x_from;
  y_to = x_to;
  y_bins = x_bins;
  h_type = 'double';
  if length(varargin) == 1
    h_type = varargin{1};
  elseif length(varargin) == 2
    y_from = varargin{1};
    y_to = varargin{2};
    h_type = 'double';
  elseif length(varargin) == 3
    y_from = varargin{1};
    y_to = varargin{2};
    if ischar(varargin{3})
      h_type = varargin{3};
    else
      y_bins = varargin{3};
    end
  elseif length(varargin) == 4
    y_from = varargin{1};
    y_to = varargin{2};
    y_bins = varargin{3};
    h_type = varargin{4};
  end
  h = zeros(x_bins, y_bins, h_type);
  photonscore_mex(photonscore.Const.FN_HIST, ...
      h, x, x_from, x_to, y, y_from, y_to);
else
  error('Inconsistent number of inputs');
end

