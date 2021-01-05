function plot_fit_decay_model(model, irf, decay, ch)
if nargin == 3
  ch = 1000
end
x = 1:length(decay);
x = x*ch/1000;
subplot(4,1,[1 2 3]);
semilogy(x, irf);
hold on;
semilogy(x, decay);
semilogy(x, model.model);
legend 'IRF' 'Decay' 'Model'
hold off;
grid on
xlim([min(x) max(x)]);

subplot(4,1,4);
plot(x, model.residuals);
legend 'Residuals'
grid on
ylim([-6 6]);
xlim([min(x) max(x)]);
end

