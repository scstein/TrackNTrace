function f = MLE_EXC_PATT_FUNC_SEP(x, img, SEPdata, field)

% xcent = x(1);
% ycent = x(2);
% be = x(3);
% al = x(4);
% Amp = x(5);
% BG = x(6);

% Compute MLE function
im_theo = x(5)*RadialPatternPos_SEP(SEPdata, x(1), x(2), x(4), x(3), field) + x(6);
% mim(im_theo);
% pause

TINY = 1e-30;
f = img.*log(im_theo + TINY) - im_theo;
f = -sum(f(:));
end

