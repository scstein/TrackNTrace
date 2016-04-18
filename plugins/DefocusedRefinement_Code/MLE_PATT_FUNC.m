function f = MLE_PATT_FUNC(x, pars)

% img = pars.img;
% field = pars.field;
% SetupParams = pars.SetupParams;

% xcent = x(1);
% ycent = x(2);
% be = x(3);
% al = x(4);
% Amp = x(5);
% BG = x(6);
% NA = x(7)
% mag = x(8)
% focus = x(9)
% pixelsize = x(10)

% Compute MLE function
% im_theo = x(5)*PatternGenerateSimple_SEP( x(1), x(2), x(3), x(4), pars.field, pars.SetupParams, pars.SEPdata) + x(6);
[ im_theo ] = x(5)*PatternGenerateSimple( x(1), x(2), x(3), x(4), pars.NA,  pars.lamem, x(7), x(8), pars.pixelsize, pars.nn, pars.field) + x(6);
% mim(im_theo);

TINY = 1e-30;
f = pars.img.*log(im_theo + TINY) - im_theo;
f = -sum(f(:));
end

