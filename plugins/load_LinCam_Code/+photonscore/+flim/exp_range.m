function r = exp_range(from, to, count)
% EXP_RANGE generates exponentially spaced range
r = 0:(count-1);
r = exp(r / r(end)) - 1;
r = r / r(end) * (to - from) + from;
end

