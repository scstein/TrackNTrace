function [x1 model] = fit_decay(irf, decay, x, from, to)
% FIT
%   IRF vector of IRF
%   DECAY decay counts
%   X initial guess
%   FROM, TO lower and upper bounds constraints
%
%   X, FROM and TO fields:
%     irf_shift
%     background
%     a
%     tau
%     tau_ref
irf = double(irf);
decay = int32(decay);

[components, fits] = size(x.tau);

if fits ~= 1 && fits ~= size(decay, 2)
  error('Number of fits is inconsistent');
end

x0 = zeros(3+2*components, fits);
lo = zeros(3+2*components, fits);
up = zeros(3+2*components, fits);

scalars = {
    {1, 'irf_shift',  0},...
    {2, 'tau_ref',    0},...
    {3, 'background', 0},...
    {4, 'a',          1},...
    {5, 'tau',        10}...
};

a_t_range = (0:(components-1))*2;

for i=1:length(scalars)
  s = scalars{i};
  ix = s{1};
  name = s{2};
  x_value = GetField(x, name, s{3});
  lo_value = GetField(from, name, x_value);
  up_value = GetField(to, name, x_value);
  if strcmp(name, 'a') || strcmp(name, 'tau')
    ix = ix + a_t_range;
  end
  x0(ix, :) = x_value;
  lo(ix, :) = lo_value;
  up(ix, :) = up_value;
end

fit = photonscore_mex(photonscore.Const.FN_FLIM_FIT_DECAY, ...
    irf, decay, lo, x0, up, components);

x1.irf_shift = fit(1, :);
x1.tau_ref = fit(2, :);
x1.background = fit(3, :);
x1.a = fit(4 + a_t_range, :);
x1.tau = fit(5 + a_t_range, :);
x1.likelihood = fit(end, :);

if nargout == 2
  model = photonscore.flim.model_of_fit_decay(irf, decay, x1);
end

end

function value = GetField(struct, field_name, default_value)
try
    value = struct.(field_name);
catch
    value = default_value;
end
end
