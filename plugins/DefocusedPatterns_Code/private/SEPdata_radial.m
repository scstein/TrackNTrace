function SEPdata = SEPdata_radial(pixel, NA, n0, n, n1, d0, d, d1, lamex,  focpos, atf, ring)

if nargin<12 || isempty(ring)
    ring = 'cos(psi).*rad';
end
    
if nargin<11 || isempty(atf)
    atf = [];
end

if nargin<10 || isempty(focpos)
    focpos = 0;
end

if nargin<9 || isempty(lamex)
    lamex = 0.64; % 640 nm laser excitation
end


if nargin<8 || isempty(d1)
    d1 = [];
end

if nargin<7 || isempty(d)
    d = 0.001;
end

if nargin<6 || isempty(d0)
    d0 = [];
end

if nargin<5 ||isempty(n1)
    n1 = 1;
end

if nargin<4 || isempty(n)
    n = 1;
end

if nargin<3 || isempty(n0)
    n0 = 1.52;
end

if nargin<2 || isempty(NA)
    NA = 1.49;
end

if nargin<1 || isempty(pixel)
    pixel = 0.04;
end

nn = round(0.6./pixel);
resolution = [lamex/0.02 lamex/0.001];
rhofield = [-lamex/resolution(1)/2 nn(1)*pixel(1)*1.1];
zfield = [0 0.01];
fd = 3e3;

over = inf;
maxm = 3;

[fxc, fxs, fyc, fys, fzc, fzs, rho, ~] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);            %#ok<*AGROW>
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
SEPdata.fxc   = fxc + fxc2;
SEPdata.fxs   = fxs + fxs2;
SEPdata.fyc   = fyc + fyc2;
SEPdata.fys   = fys + fys2;
SEPdata.fzc   = fzc + fzc2;
SEPdata.fzs   = fzs + fzs2;
SEPdata.rho   = rho;
SEPdata.nn    = nn;
SEPdata.pixel = pixel;

