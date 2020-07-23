function [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm)
                                                           
% function [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield,
% zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf,
% resolution, ring, maxm) calculates the electric field distribution of a 
% focused laser beam 
%
% input arguments:
% rhofield - two-element vector containing the minimum and maximum value of the rho coordinate
% zfield   - two-element vector containing the minimum and maximum value of the z coordinate
% NA       - numerical aperture
% fd       - focal distance
% n0       - refractive index vector of layers below sample
% n        - refractive index of sample medium
% n1       - refractive index vector of layers above sample
% d0       - thickness vector of layers below sample
% d        - thickness of sample medium
% d1       - thickness vector of layers above sample
% lamex    - excitation wavelength
% over     - laser characteristics: [beam waist position, beam waist radius, astigmatism value]
% focpos   - focus position [z, x displacement, y displacment, x inclination, y inclination]
% atf      - cover slide correction
%
% output arguments:
% exc      - excitation intensity distribution
% rho, z   - rho and z coordinate grid
% f**      - electric field component distributions
%
% field is given by
% ex = fxc(:,:,1) + fxc(:,:,2)*cos(phi) + ... + fxs(:,:,1)*sin(phi) + ...
% etc.

maxnum = 1e3;
if nargin < 17 || isempty(maxm)
    maxm = 2;
end
if nargin<15 || isempty(resolution)
    resolution = [20 20];
end
if length(resolution)==1
    resolution = resolution*[1 1];
end

if isempty(d) || d<max(zfield(:))
    d = max(zfield(:));
end

n0 = n0(:).'; n1 = n1(:).'; d0 = d0(:).'; d1 = d1(:).';

k0 = 2*pi/lamex*n0; k = 2*pi/lamex*n; 
d0 = 2*pi*d0/lamex; d1 = 2*pi*d1/lamex;

if length(rhofield)==1
    dz = lamex/resolution(2); 
    rhov = rhofield; 
    if length(zfield)>1
        zv = zfield(1) + (0.5:(zfield(2)-zfield(1))/dz)*dz;
    else
        zv = zfield;
    end
    if isempty(rhov)
        rhov = rhofield;
    end
    if isempty(zv)
        zv = zfield;
    end
    rho = rhov;
    z = zv;
elseif length(rhofield)==2
    drho = lamex/resolution(1);
    dz = lamex/resolution(2);
    rhov = rhofield(1) + (0.5:(rhofield(2)-rhofield(1))/drho)*drho; 
    if length(zfield)>1
        zv = zfield(1) + (0.5:(zfield(2)-zfield(1))/dz)*dz;
    else
        zv = zfield;
    end
    if isempty(rhov)
        rhov = rhofield;
    end
    if isempty(zv)
        zv = zfield;
    end
    [rho, z] = ndgrid(rhov, zv);
else
    rhov = rhofield(:,1)'; zv = zfield(1,:);
    rho = rhofield; z = zfield;
end

chimax = asin(NA/n0(1)); 
chimin = 0; 
dchi = (chimax-chimin)/maxnum;
chi = (chimin:dchi:chimax)';
c0 = cos(chi);
s0 = sin(chi);
ss = n0(1)/n*s0;
cc = sqrt(1-ss.^2);

if length(over)==1 % over = ww;
    ww = over;
    fac = s0.*sqrt(c0);
    om = [1 1]/ww^2;
elseif length(over)==2 % over = [z0, ww]
    if over(2)==0 % wide field illumination - not an elegant solution but it works
        fac = 0*chi;    
        fac(1) = 1;
        om = [1 1];
    else
        if over(1)==inf
            zeta = 0;
            ww = over(2);
        else
            zeta = over(1)*lamex/pi/over(2)^2;
            ww = over(2)*sqrt(1+zeta^2);
        end
        fac = s0.*sqrt(c0);
        om = [1 1]*(1-1i*zeta)/ww^2;
    end
elseif length(over)==3 % over = [ast, w0, 1]
    zeta = 2*over(1)*pi*over(2)^2/(NA*fd)^2;
    w0 = over(2)/sqrt(1+zeta^2);
    zeta = [-zeta zeta];
    om = (1-1i*zeta)./(1+zeta.^2)./w0^2/2;
    fac = s0.*sqrt(c0);
else % over = [z1, z2, w1, w2]
    zeta = over(1:2)*lamex/pi./over(3:4).^2;
    om = (1-1i*zeta)./(1+zeta.^2)./over(3:4).^2/2;
    fac = s0.*sqrt(c0);
end

[rp0,rs0,tp0,ts0] = Fresnel(n0(1)*c0,[n0 n],d0);
[rp0,rs0] = Fresnel(n*cc,[n n0(end:-1:1)],d0(end:-1:1));
[rp1,rs1] = Fresnel(n*cc,[n n1],d1);
if nargin>13 && ~isempty(atf)
    if length(atf)==1 % accounting for reflection losses when using water immersion; atf = ref index of cover slide
        [tmp, tms, tmp, tms] = Fresnel(n0(1)*c0,n0(1),atf);
        [mp, ms, mp, ms] = Fresnel(sqrt(atf^2-n0(1)^2*s0.^2),atf,n0(1));
    else % + aberration for water immersion; atf(2) = thickness mismatch
        [tmp, tms, tmp, tms] = Fresnel(n0(1)*c0,[n0(1) n0(1) atf(1)], -2*pi*atf(2)/lamex);        
        [mp, ms, mp, ms] = Fresnel(sqrt(atf(1)^2-n0(1)^2*s0.^2),[atf(1) atf(1) n], 2*pi*atf(2)/lamex);        
    end
    tp0 = tmp.*mp.*tp0;
    ts0 = tms.*ms.*ts0;
else
    atf = [];
end
tp = tp0./(1-rp1.*rp0.*exp(2*1i*k*cc.'*d));
ts = ts0./(1-rs1.*rs0.*exp(2*1i*k*cc.'*d));
tp = tp.';
ts = ts.';
rp = tp.*rp1.';
rs = ts.*rs1.';

phase1 = k*cc*zv - k0(1)*c0*focpos(1)*ones(size(zv));
phase2 = k*cc*(2*d-zv) - k0(1)*c0*focpos(1)*ones(size(zv));
% phase1 = k*cc*(zv - focpos(1));
% phase2 = k*cc*(2*d-zv-focpos(1));
ez1 = exp(1i*phase1);
ez2 = exp(1i*phase2);

barg = k*rhov'*ss';
row = ones(size(zv));

if length(focpos)<5
    focpos = [focpos zeros(1,5-length(focpos))];
end
rho0 = sqrt(sum(focpos(2:3).^2));
psi0 = angle(focpos(2)+1i*focpos(3));
rad = s0/max(s0);
if nargin>15 && ~isempty(ring) % ring has to be a function of rad and psi
    fun = ['fac.*exp(-om(1)*(fd*n*ss*cp-focpos(4)).^2 - om(2)*(fd*n*ss*sp-focpos(5)).^2 + 1i*k*ss*rho0*cos(psi-psi0)).*(' ring ');'];
else
    fun = 'fac.*exp(-om(1)*(fd*n*ss*cp-focpos(4)).^2 - om(2)*(fd*n*ss*sp-focpos(5)).^2 + 1i*k*ss*rho0*cos(psi-psi0));';
    ring = [];
end
for j=0:2*maxm
    psi = j/(2*maxm+1)*2*pi; cp = cos(psi); sp = sin(psi);
    eval(['ef = ' fun]);
    tmpxt(:,j+1) = (cp^2*tp.*cc+sp^2*ts).*ef;
    tmpxr(:,j+1) = (-cp^2*rp.*cc+sp^2*rs).*ef;
    tmpyt(:,j+1) = cp*sp*(tp.*cc-ts).*ef;
    tmpyr(:,j+1) = cp*sp*(-rp.*cc-rs).*ef;
    tmpzt(:,j+1) = -cp*ss.*tp.*ef;
    tmpzr(:,j+1) = -cp*ss.*rp.*ef;
end
tmpxt = fft(tmpxt,[],2)/(maxm+0.5);
tmpxr = fft(tmpxr,[],2)/(maxm+0.5);
tmpyt = fft(tmpyt,[],2)/(maxm+0.5);
tmpyr = fft(tmpyr,[],2)/(maxm+0.5);
tmpzt = fft(tmpzt,[],2)/(maxm+0.5);
tmpzr = fft(tmpzr,[],2)/(maxm+0.5);

for j=0:maxm
    jj = besselj(j,barg);
    if j==0
        fxc(:,:,1) = jj*((tmpxt(:,1)*row).*ez1 + (tmpxr(:,1)*row).*ez2);
        fyc(:,:,1) = jj*((tmpyt(:,1)*row).*ez1 + (tmpyr(:,1)*row).*ez2);
        fzc(:,:,1) = jj*((tmpzt(:,1)*row).*ez1 + (tmpzr(:,1)*row).*ez2);
    else
        fxc(:,:,j+1) = 1i^(-j)*jj*(((tmpxt(:,j+1)+tmpxt(:,end-j+1))*row).*ez1 + ((tmpxr(:,j+1)+tmpxr(:,end-j+1))*row).*ez2);
        fyc(:,:,j+1) = 1i^(-j)*jj*(((tmpyt(:,j+1)+tmpyt(:,end-j+1))*row).*ez1 + ((tmpyr(:,j+1)+tmpyr(:,end-j+1))*row).*ez2);
        fzc(:,:,j+1) = 1i^(-j)*jj*(((tmpzt(:,j+1)+tmpzt(:,end-j+1))*row).*ez1 + ((tmpzr(:,j+1)+tmpzr(:,end-j+1))*row).*ez2);
        fxs(:,:,j) = 1i^(-j+1)*jj*(((tmpxt(:,j+1)-tmpxt(:,end-j+1))*row).*ez1 + ((tmpxr(:,j+1)-tmpxr(:,end-j+1))*row).*ez2);
        fys(:,:,j) = 1i^(-j+1)*jj*(((tmpyt(:,j+1)-tmpyt(:,end-j+1))*row).*ez1 + ((tmpyr(:,j+1)-tmpyr(:,end-j+1))*row).*ez2);
        fzs(:,:,j) = 1i^(-j+1)*jj*(((tmpzt(:,j+1)-tmpzt(:,end-j+1))*row).*ez1 + ((tmpzr(:,j+1)-tmpzr(:,end-j+1))*row).*ez2);
    end
end

if nargout==1
    exc.maxnum = maxnum;
    exc.NA = NA;
    exc.fd = fd;
    exc.n0 = n0;
    exc.n = n;
    exc.n1 = n1;
    exc.d0 = d0/2/pi*lamex;
    exc.d = d;
    exc.d1 = d1/2/pi*lamex;
    exc.lambda = lamex;
    exc.over = over;
    exc.focpos = focpos;
    exc.atf = atf;
    exc.resolution = resolution;
    exc.ring = ring;
    exc.maxm = maxm;
    exc.rho = rho;
    exc.z = z;
    exc.fxc = fxc;
    exc.fyc = fyc;
    exc.fzc = fzc;
    exc.fxs = fxs;
    exc.fys = fys;
    exc.fzs = fzs;
    fxc = exc;
end

% FocusImage2D(exc.rho,exc.z,cat(4,cat(3,exc.fxc,exc.fxs),cat(3,exc.fyc,exc.fys),cat(3,exc.fzc,exc.fzs)))
% [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc([0 8], 0, 1.2, 3e3, 1.33, 1.33, 1.33, [], 0, [], 0.635, 700/5, 0); plot(rho,abs(fxc(:,1)).^2+abs(fyc(:,1)).^2+abs(fzc(:,1)).^2)