function o = sort(x, x_min, x_max, x_bins, y, y_min, y_max, y_bins, dt)
if nargin == 6
  % sort(x, from, to, bins, y, dt)
  ax = x;
  ax_min = x_min;
  ax_max = x_max;
  ax_bins = x_bins;
  ay = y;
  ay_min = x_min;
  ay_max = x_max;
  ay_bins = x_bins;
  adt = y_min;
elseif nargin == 9
  % full form
  ax = x;
  ax_min = x_min;
  ax_max = x_max;
  ax_bins = x_bins;
  ay = y;
  ay_min = y_min;
  ay_max = y_max;
  ay_bins = y_bins;
  adt = dt;
else
  error 'wrong number of arguments';
end

[h, t] = photonscore_mex(photonscore.Const.FN_FLIM_SORT, ...
    ax, ax_min, ax_max, ax_bins, ...
    ay, ay_min, ay_max, ay_bins, ...
    adt);

o.image = h;
o.time = t;
