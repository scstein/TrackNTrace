function res = model_of_fit_decay(irf, decay, x1)
model = zeros(size(decay));
residuals = zeros(size(decay));
for i=1:size(model, 2)
  if size(irf, 2) > 1
    irf_i = irf(:,i);
  else
    irf_i = irf;
  end
  m = photonscore.flim.convolve(irf_i/sum(irf_i(:)),...
          x1.irf_shift(i),...
          size(model, 1), x1.tau(:, i), x1.tau_ref(i));
  model(:, i) = m * x1.a(:, i) + x1.background(i);
  decay_i = double(decay(:, i));
  residuals(:, i) = model(:, i) - decay_i;
  s = decay_i > 0;
  residuals(s, i) = residuals(s, i) ./ sqrt(decay_i(s));
end
res.model = model;
res.residuals = residuals;
end

