function [x, m] = fit_decay_auto(r, d, exps)
batch = 100;
win = 4;
trials = 5;
rep = batch / win;
bg = mean(d(1:10));

% First trial
xl.background = ones(1, batch) * bg - sqrt(bg);
x0.background = ones(1, batch) * bg;
xu.background = ones(1, batch) * bg + sqrt(bg);

xl.tau_ref = ones(1, batch) * 0;
x0.tau_ref = ones(1, batch) * 2;
xu.tau_ref = ones(1, batch) * length(d) / 10;

xl.irf_shift = zeros(1, batch) - 20;
x0.irf_shift = zeros(1, batch) .* randn(1, batch);
xu.irf_shift = zeros(1, batch) + 20;

% 3-component fit
xl.tau = ones(exps, batch) * 3;
x0.tau = rand(exps, batch) * length(d) / 2 + 10;
xu.tau = ones(exps, batch) * length(d);

x0.a = ones(exps, batch);

x1 = photonscore.flim.fit_decay(...
    repmat(r(:), 1, batch), repmat(d(:), 1, batch),...
    x0, xl, xu);

for trial=1:trials
    s = (x1.likelihood == -1);
    x1.likelihood(s) = 10^22;
    [~, il] = sort(x1.likelihood);
    il = il(1:win);

    x0 = copy(x1, il, rep);
    if exist('x')
        if x0.likelihood(1) < x.likelihood
            x = copy(x0, 1, 1);
        end
    else
        x = copy(x0, 1, 1);
    end
    x0.irf_shift = x0.irf_shift .* (1 + randn(1, batch)/10);
    x0.irf_shift = min(xu.irf_shift, x0.irf_shift);
    x0.irf_shift = max(xl.irf_shift, x0.irf_shift);

    x0.tau_ref = x0.tau_ref .* (1 + randn(1, batch)/10);
    x0.tau_ref = min(xu.tau_ref, x0.tau_ref);
    x0.tau_ref = max(xl.tau_ref, x0.tau_ref);

    x0.background = ones(1, batch) * bg;
    x0.a = ones(exps, batch);

    x0.tau = x0.tau .* (1 + randn(exps, batch)/exps);
    x0.tau = min(xu.tau, x0.tau);
    x0.tau = max(xl.tau, x0.tau);

    x1 = photonscore.flim.fit_decay(...
        repmat(r(:), 1, batch), repmat(d(:), 1, batch),...
        x0, xl, xu);
end

[~, il] = min(x1.likelihood);
x = copy(x1, il, 1);

m = photonscore.flim.model_of_fit_decay(r, d, x);
end

function r = copy(x, ix, rep)
r.irf_shift = repmat(x.irf_shift(ix), 1, rep);
r.tau_ref = repmat(x.tau_ref(ix), 1, rep);
r.background = repmat(x.background(ix), 1, rep);
r.a = repmat(x.a(:, ix), 1, rep);
r.tau = repmat(x.tau(:, ix), 1, rep);
r.likelihood = repmat(x.likelihood(ix), 1, rep);
end
