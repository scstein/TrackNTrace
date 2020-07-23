function [int, Sx, Sy] = RadialPatternPos_SEP(SEPdata, xcent, ycent, al, be, field, pic)

if nargin<7 || isempty(pic)
    pic = 0;
end
if nargin<6 || isempty(field)
    field = 100;
elseif numel(field) == 1
    field = [field field];
end

if nargin<5 || isempty(be)
    be = 0; % in plane angle
end

if nargin<4 || isempty(al)
    al = pi/2; % out of plane angle
end

if nargin<3 || isempty(ycent)
    ycent = 0;
end

if nargin<2 || isempty(xcent)
    xcent = 0;
end


fxc     = SEPdata.fxc;
fxs     = SEPdata.fxs;
fyc     = SEPdata.fyc;
fys     = SEPdata.fys;
fzc     = SEPdata.fzc;
fzs     = SEPdata.fzs;
rho     = SEPdata.rho;
nn      = SEPdata.nn;
pixel   = SEPdata.pixel;

maxm = 3;

bck = Disk(nn);
[Sx,Sy] = meshgrid(-field(2)-xcent:field(2)-xcent,-field(1)+ycent:field(1)+ycent);
[x, y] = meshgrid((-nn:nn),(-nn:nn));
% rr = pixel*sqrt((x-rem(xcent,1)).^2+(y-rem(ycent,1)).^2);
% psi = angle((x-rem(xcent,1)) + 1i*(y-rem(xcent,1)));
rr = pixel*sqrt(x.^2+y.^2);
psi = angle(x + 1i*y);
mask(max(nn)-nn+1:max(nn)-nn+2*nn+1,max(nn)-nn+1:max(nn)-nn+2*nn+1) = interp1(rho(:,1),(fxc(:,1,1)*sin(al)*cos(be)+fyc(:,1,1)*sin(al)*sin(be)+fzc(:,1,1)*cos(al)) ,rr,'cubic',0);
for j=1:maxm
    mask(max(nn)-nn+1:max(nn)-nn+2*nn+1,max(nn)-nn+1:max(nn)-nn+2*nn+1) = mask(max(nn)-nn+1:max(nn)-nn+2*nn+1,max(nn)-nn+1:max(nn)-nn+2*nn+1) + ...
            interp1(rho(:,1),(fxc(:,1,j+1)*sin(al)*cos(be)+fyc(:,1,j+1)*sin(al)*sin(be)+fzc(:,1,j+1)*cos(al)) ,rr,'cubic',0).*cos(j*psi) + ...
            interp1(rho(:,1),(fxs(:,1,j)*sin(al)*cos(be)+fys(:,1,j)*sin(al)*sin(be)+fzs(:,1,j)*cos(al)) ,rr,'cubic',0).*sin(j*psi);
end

int = abs(mask).^2;
int = int.*bck;
int = int/sum(int(:));


% % Whole pixel shift to final position
tmpint = zeros(2*field+1);
tmpint(1:size(int,1), 1:size(int,2)) = int;

% % Supixel shifts by interpolation usint interp2
xshift = field(2)-floor(size(int,2)/2)+xcent;
yshift = field(1)-floor(size(int,1)/2)-ycent;
[xi, yi] = meshgrid([1:size(tmpint,2)]-xshift, [1:size(tmpint,1)]-yshift);
int = interp2(tmpint, xi, yi, 'linear', 0');
int = int/sum(int(:)); % Normalize Integral to 1


if pic
    mim(int)
end
end