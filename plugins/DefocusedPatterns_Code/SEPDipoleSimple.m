function [ SEPdata, SetupParams ] = SEPDipoleSimple(nn, pixelsize, NA, lamem, mag, focus)

if nargin<1 || isempty(nn)
    nn = 30;
end

if nargin<2 || isempty(pixelsize)
    pixelsize = 24;
end

if nargin<3 || isempty(NA)
    NA = 1.49;
end

if nargin<4 || isempty(lamem)
    lamem = 0.690;
end

if nargin<5 || isempty(mag)
    mag = 400;
end

if nargin<6 || isempty(focus)
    focus = 0.9;
end

if nargin<7
    atf = [];
end
if nargin<8
    ring = [];
end

% These paremeters we do not touch
d1 = [];
d = 0.001;
d0 = [];
n1 = 1;
n = 1;
n0 = 1.52;
z = 0;

[~, ~, ~, rho, ~, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 1.5*max(nn)*pixelsize/mag], z, NA, n0, n, n1, d0, d, d1, lamem, mag, focus, atf, ring);
SEPdata.rho = rho;
SEPdata.fxx0 = fxx0;
SEPdata.fxx2 = fxx2;
SEPdata.fxz = fxz;
SEPdata.byx0 = byx0;
SEPdata.byx2 = byx2;
SEPdata.byz = byz;

SetupParams.nn = nn;
SetupParams.pixelsize = pixelsize;
SetupParams.NA = NA;
SetupParams.lamem = lamem;
SetupParams.mag = mag;
SetupParams.focus = focus;
end

