function [ x ] = ITFUNC_PATT_VIS_SEP( x, pars )
% Visualization of the pattern with variables x and parameters pars.

mim( x(5)*PatternGenerateSimple_SEP( x(1), x(2), x(3), x(4), pars.field, pars.SetupParams, pars.SEPdata) + x(6));
% compare_images(pars.img, x(5)*PatternGenerateSimple_SEP( x(1), x(2), x(3), x(4), pars.field, pars.SetupParams, pars.SEPdata) + x(6));
end

