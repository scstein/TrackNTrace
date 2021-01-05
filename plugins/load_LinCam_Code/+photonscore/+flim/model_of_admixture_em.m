function m = model_of_admixture_em(model, decay, x1)
m.model = model * x1;
m.residuals = (m.model - decay);
s = decay > 0;
m.residuals(s) = m.residuals(s) ./ sqrt(decay(s));
m.decay = decay;
end

