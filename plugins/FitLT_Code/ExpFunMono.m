function [err, z] = ExpFunMono(p, t, y, mle)
%[err, z] = ExpFunMono(p, t, y, mle), minimised, faster version of ExpFunFull
% p   - parameters (bg tau1 amp1 tau2 amp2 ...)
% t   - time (x) axis, does not have to start at zero but has to be unifrom
% y   - function values (y-axis)
% mle - if true the MLE (negative log likelihood) is returned, otherwise chi^2
%
% err - error: negative log likelihood
% z   - total function value
%
% Christoph Thiele, 2020

t = t(:);
p = p(:);
y = y(:);

n = numel(t);
dt = (t(end)-t(1))/n; % assumes uniform time axis

z = exp(-t/p(2));
% % z = p(3)*z/sum(z) + p(1)/n;     % Norm numerically
z = p(3)/p(2)*dt*z + p(1)/n;         % Norm analytically

if mle
    % MLE
    ind = y>0;
    if all(ind)
        err = sum(-y.*log(z+1e-30)+z);
    else
        err = sum(-y(ind).*log(z(ind)+1e-30))+sum(z);% Uses the normalisation: sum(z)==p(1)+p(3)
    end
else
    % chi2
    err = sum((y-z).^2./abs(z))/numel(t);
end
% Plot for debugging
% semilogy(t, y, 'o', t, z); drawnow
end

