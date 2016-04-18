function [ x ] = ITFUNC_PATT_VIS( x, pars )
% Visualization of the pattern with variables x and parameters pars.
mim( x(5)*PatternGenerateSimple( x(1), x(2), x(3), x(4), pars.NA,  pars.lamem, x(7), x(8), pars.pixelsize, pars.nn, pars.field) + x(6) );
x'
end

