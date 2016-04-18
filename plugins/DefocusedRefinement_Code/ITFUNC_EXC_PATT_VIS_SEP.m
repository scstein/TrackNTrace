function [ x ] = ITFUNC_EXC_PATT_VIS_SEP( x, ~, SEPdata, field )
% Visualization of the pattern with variables x and parameters pars.

mim( x(5)*RadialPatternPos_SEP(SEPdata, x(1), x(2), x(4), x(3), field) + x(6));
end

