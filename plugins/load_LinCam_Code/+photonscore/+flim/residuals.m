function [r, lsr, chi] = residuals(n, m)
% RESIDUALS returs weighted residuals vector, L*r and chi^2
%   r = p.residuals(n, m)
%
%   N counts
%   M model
s = n > 0;
r = (n - m);
r(s) = r(s) ./ sqrt(n(s));
if nargout > 1
  s = (n > 0) & (m > 0);
  lsr = 2 * sum(n(s) .* log(n(s) ./ m(s))) / length(n);
end
if nargout > 2
  chi = sum(r.^2) / length(n);
end
end

