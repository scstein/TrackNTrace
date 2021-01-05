function x1 = tail_fit(decays, x0, iterations, accelerate)
if nargin == 4
  arg_iterations = iterations;
  arg_accelerate = accelerate;
end

if nargin == 3
  arg_iterations = iterations;
  arg_accelerate = 1;
end

if nargin == 2
  arg_iterations = 100;
  arg_accelerate = 1;
end

if size(decays, 1) == 1
  decays = decays(:);
end
fits = size(decays, 2);

x0 = x0(:);
if mod(length(x0), fits) ~=0
  error 'Inconsistent X0 size'
end

x1 = photonscore_mex(photonscore.Const.FN_FLIM_TAIL_FIT, ...
    double(decays), double(x0), arg_iterations, arg_accelerate);
end
