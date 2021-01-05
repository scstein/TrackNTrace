function [x1 model] = admixture_em(model, counts, x0, varargin)
% ADMIXTURE_EM resolve COUNTS mixture with a MODEL starting from X0
%
%   MODEL matrix of size CHANNELS x COMPONENTS
%   COUNTS integer matrix of size CHANNELS x RHS
%   X0 initial guess of size COMPONENTS x RHS
%
%   Optional arguments:
%     'iterations', default 500
%     'method', 0, 1, 2 or default 3
%     'inf', infinite pixel compensation, default is (1-sum(model, 1))
ip = inputParser;
addRequired(ip,'model');
addRequired(ip,'counts');
addRequired(ip,'x0');
addParameter(ip, 'iterations', 500);
addParameter(ip, 'method', 3);
addParameter(ip, 'inf', 1 - double(sum(model, 1)'));
parse(ip, model, counts, x0, varargin{:});
p = ip.Results;

if ~isa(p.model, 'double')
    p.model = double(model);
end

if ~isa(p.counts, 'int32')
    p.counts = int32(p.counts);
end

if ~isa(p.x0, 'double')
    p.x0 = double(p.x0(:));
end

channels = size(p.model, 1);
components = size(p.x0, 1);
rhs = size(p.counts, 2);

if components ~= size(p.x0, 1) || rhs ~= size(p.x0, 2)
  error('Incorrect X0 size, expected (%d, %d), got (%d, %d)',...
        components, rhs, size(p.x0, 1), size(p.x0, 2));
end

if channels ~= size(p.counts, 1) || rhs ~= size(p.counts, 2)
  error('Incorrect COUNTS size, expected (%d, %d), got (%d, %d)',...
        channels, rhs, size(p.counts, 1), size(p.counts, 2));
end

% The number of columns of the MODEL can be COMPONENTS or RHS * COMPONENTS
if ~(components == size(p.model, 2) || rhs * components == size(p.model, 2))
  error('Incorrect MODEL size, expected (%d, %d) or (%d, %d), got (%d, %d)',...
        components, rhs, ...
        components, components * rhs, ...
        size(p.inf, 1), size(p.inf, 2));
end

if size(p.model, 2) ~= size(p.inf, 1) || 1 ~= size(p.inf, 2)
  error('Incorrect INF size (%d, 1), got (%d, %d)',...
        size(p.model, 2),...
        size(p.inf, 1), size(p.inf, 2));
end

model_size = numel(model) / size(x0, 2);
if model_size < 10000
  % Go with parallelized for a small number model sizes
  x1 = photonscore_mex(photonscore.Const.FN_FLIM_AM_EM, ...
        p.model, p.inf, p.counts, p.x0, ...
        p.iterations, p.method);
else
  % FIXME: call MKL backed solver for larger-scale models
  x1 = photonscore_mex_xl(0, ...
        p.model, p.inf, p.counts, p.x0, ...
        p.iterations, p.method);
end

if nargout == 2
  model = photonscore.flim.model_of_admixture_em(model, counts, x1);
end

end
