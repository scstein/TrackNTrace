function [err, c, zz, z] = ExpFunLSQ(p, t, y, mle)
%[err, z] = ExpFunLSQ(p, t, y, mle), minimised, faster version of ExpFun
% Does a nonnegative least-square fit of y with the lifetimes p.
% p   - parameters (tau1 tau2 ...)
% t   - time (x) axis, does not have to start at zero but has to be unifrom
% y   - function values (y-axis)
% mle - if true the MLE (negative log likelihood) is returned, otherwise chi^2
%
% err - error: negative log likelihood OR chi^2
% c   - compontens amplitudes (bg amp1 amp2 ...)
% zz  - compontens function value
% z   - total function value
%
% Christoph Thiele, 2020


t = t(:);
p = p(:)';
y = y(:);

ind = y>0;
n = numel(t);
dt = (t(end)-t(1))/n; % assumes uniform time axis

zz = [ones(size(y))/n, exp(-t./p)./p*dt];
c = lsqnonneg(zz,y);
z = zz*c;

if mle
    % MLE
    if all(ind)
        err = sum(-y.*log(z+1e-30)+z);
    else
        err = sum(-y(ind).*log(z(ind)+1e-30))+sum(z);% For t(1)==0 -> sum(z)==p(1)+p(3)
    end
else
    % chi2
    err = sum((y-z).^2./abs(z))/numel(t);
end
    
end

