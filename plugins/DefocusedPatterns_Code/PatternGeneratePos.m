function [int,Sx,Sy] = PatternGeneratePos(xcent,ycent,z, NA, n0, n, n1, d0, d, d1, lamem, mag, focus, atf, ring, pixel, nn, field, be, al, pic)

if nargin<16 || isempty(pixel)
    pixel = 24;
end
if nargin<14
    atf = [];
end
if nargin<15
    ring = [];
end
if nargin<21 || isempty(pic)
    pic = 0;
end
if nargin<18 || isempty(field)
    field = 100;
end
if nargin<19 || isempty(be)
    be= 0; % in-plane angle
end
if nargin<20 || isempty(al)
    al = pi/2; %  out-of-plane angle
end
if nargin<17 || isempty(nn)
    nn = 30;
end
if isempty(xcent)
    xcent = 0;
end
if isempty(ycent)
    ycent = 0;
end

if nargin<13 || isempty(focus)
    focus = 0.9;
end

if nargin<12 || isempty(mag)
    mag = 400;
end

if nargin<11 || isempty(lamem)
    lamem = 0.690;
end

if nargin<10 || isempty(d1)
    d1 = [];
end

if nargin<9 || isempty(d)
    d = 0.001;
end

if nargin<8 || isempty(d0)
    d0 = [];
end

if nargin<7 ||isempty(n1)
    n1 = 1;
end

if nargin<6 || isempty(n)
    n = 1;
end

if nargin<5 || isempty(n0)
    n0 = 1.52;
end

if nargin<4 || isempty(NA)
    NA = 1.49;
end
if numel(field)==1
    field = [field,field];
end


bck = Disk(nn);

fprintf('NA %f, mag %f, focus %f\n', NA, mag, focus)
[~, ~, ~, rho, ~, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 1.5*max(nn)*pixel/mag], z, NA, n0, n, n1, d0, d, d1, lamem, mag, focus, atf, ring);


[Sx,Sy] = meshgrid(-field(2)-xcent:field(2)-xcent,-field(1)+ycent:field(1)+ycent);
% [x y] = meshgrid(max(1,floor(Sx(1)+nn+1)):min(floor(Sx(end)-nn-1),2*nn+1),max(1,floor(Sy(1)+nn+1)):min(floor(Sy(end)-nn-1),2*nn+1));
% x = x-nn-1;
% y = y-nn-1;
% [x, y] = meshgrid((-nn:nn) + rem(xcent,1),(-nn:nn) + rem(ycent,1)); % subpixel shift
[x y] = meshgrid((-nn:nn),(-nn:nn));
p = angle(x+1i*y);
r = sqrt(x.^2+y.^2);
rho = rho/pixel;
ex = sin(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,fxx2,r,'cubic'))+cos(al)*cos(p-be).*interp1(rho,fxz,r,'cubic');
by = sin(al)*(interp1(rho,byx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,byx2,r,'cubic'))+cos(al)*cos(p-be).*interp1(rho,byz,r,'cubic');
ey = sin(al)*sin(2*(p-be)).*interp1(rho,fxx2,r,'cubic')+cos(al)*sin(p-be).*interp1(rho,fxz,r,'cubic');
bx = -(sin(al)*sin(2*(p-be)).*interp1(rho,byx2,r,'cubic')+cos(al)*sin(p-be).*interp1(rho,byz,r,'cubic'));
tmp = ex*cos(be)-ey*sin(be);
ey = ey*cos(be)+ex*sin(be);
ex = tmp;
tmp = by*cos(be)+bx*sin(be);
bx = bx*cos(be)-by*sin(be);
by = tmp;
int = real(ex.*conj(by) - ey.*conj(bx));
% int = bck.*int;

% Shift to final position
tmpint = zeros(2*field+1);
tmpint(1:size(int,1), 1:size(int,2)) = int;

% % Whole pixel shift using circshift
% xshift = field(2)-floor(size(int,2)/2)+floor(xcent);
% yshift = field(1)-floor(size(int,1)/2)-floor(ycent);
% int = circshift(tmpint,[yshift, xshift]);

% % Supixel shifts by interpolation usint imtransform
% xshift = field(2)-floor(size(int,2)/2)+xcent;
% yshift = field(1)-floor(size(int,1)/2)-ycent;
% xform = [ 1  0  0
%          0  1  0
%          xshift yshift  1 ];
% tform_translate = maketform('affine',xform);
% int = imtransform(tmpint, tform_translate,  'XData', [1 size(tmpint,2)], 'YData', [1 size(tmpint,1)]);

% % Supixel shifts by interpolation usint interp2
xshift = field(2)-floor(size(int,2)/2)+xcent;
yshift = field(1)-floor(size(int,1)/2)-ycent;
[xi, yi] = meshgrid([1:size(tmpint,2)]-xshift, [1:size(tmpint,1)]-yshift);
int = interp2(tmpint, xi, yi, 'linear', 0'); % Note cubic interpolation does not work at all (wtf?)

if pic == 1;
    mim(int)
end