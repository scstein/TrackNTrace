function f = LSQ_PATT_FUNC_SEP(x, pars)

% img = pars.img;
% field = pars.field;
% SetupParams = pars.SetupParams;
% SEPdata = pars.SEPdata;

% xcent = x(1);
% ycent = x(2);
% be = x(3);
% al = x(4);
% Amp = x(5);
% BG = x(6);

im_theo = x(5)*PatternGenerateSimple_SEP( x(1), x(2), x(3), x(4), pars.field, pars.SetupParams, pars.SEPdata) + x(6);
mim(im_theo)

f = (pars.img-im_theo).^2;
f = sum(f(:));
% f = f(:);
end

