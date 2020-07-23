function [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, orient)

% SEPDipole calculates the emission pattern of dipoles while imaged with a camera
%
% [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, orient)
%
% input arguments:
% r        - two-element vector containing the minimum and maximum value of the rho coordinate
% z        - distance of molecule from surface
% NA       - numerical aperture
% n1       - refractive index vector (bottom to top) below molecule
% n        - refractive index of molecule's layer
% n2       - refractive index vector (bottom to top) above molecule
% d1       - thickness vector for layers below molecule
% d        - thickness of embedding layer
% d2       - thickness vector for layers above molecule
% lambda   - wavelength
% mag      - magnification
% focpos   - focus position along optical axis
% atf      - cover slide correction
% ring     - phase modulation as a function of sin(theta)/sin(theta_max)
% orient   - 3D orientation as given by [theta, phi] (only intx returned)
%
% output arguments:
% int*     - intensity distributions
% f**      - electric field component distributions
% b**      - magnetic field component distributions
% rho, phi - 2d coordinate fields

maxnum = 2e3;
ni = 1.0; % it is assumed that imaging medium is air

ki = 2*pi/lambda*ni; 

if length(rho)==2
    drho = lambda/50;
    rho = rho(1) + (0:(rho(2)-rho(1))/drho)*drho;
end
rho = mag*rho; rho = rho(:)';
    
chimin = 0;
chimax = asin(NA/mag/ni);
dchi = (chimax-chimin)/maxnum;
chi = chimin + (0.5:maxnum)*dchi;

ci = cos(chi);
si = sin(chi);
s0 = ni/n1(1)*mag*si;
c0 = sqrt(1-s0.^2);
psi = asin(s0);
[v,pc,ps] = DipoleL(psi,2*pi*z/lambda,n1,n,n2,2*pi*d1/lambda,2*pi*d/lambda,2*pi*d2/lambda);
v = v.'; ps = ps.'; pc = pc.';

if nargin>12 && ~isempty(atf)
    if length(atf)==1 % accounting for reflection losses when using water immersion; atf = ref index of cover slide
        [~, ~, tmp, tms] = Fresnel(n1(1)*c0,n1(1),atf);
        [~, ~, mp, ms] = Fresnel(sqrt(atf^2-n1(1)^2*s0.^2),atf,n1(1));
    else % + aberration for water immersion; atf(2) = thickness mismatch
        if 1
            [~, ~, tmp, tms] = Fresnel(n1(1)*c0,[n1(1) atf(1) atf(1)], 2*pi*atf(2)/lambda);
            [~, ~, mp, ms] = Fresnel(sqrt(atf(1)^2-n1(1)^2*s0.^2),[atf(1) n1(1) n1(1)], -2*pi*atf(2)/lambda);
        else
            w0 = n1(1)*c0;
            w = sqrt(atf(1)^2-n1(1)^2 + w0.^2);
            mp = exp(1i*2*pi*atf(2)/lambda*(w-w0));
            ms = mp;
            tmp = 1;
            tms = 1;
        end
    end
        v = tmp.*mp.*v;
        pc = tmp.*mp.*pc;
        ps = tms.*ms.*ps;
end

phase = -n1(1).*c0*focpos;
if nargin>13 && ~isempty(ring)
    rad = s0/max(s0);
    eval(['phase = phase + ' ring ';']);
end

%alternative formulation but same result:
%phase = ni*mag*focpos*s0./c0.*si-n1(1)./c0*focpos; 
ez = exp(1i*2*pi/lambda*phase);

fac = dchi*si.*sqrt(ci./c0); % aplanatic objective

barg = ki*si'*rho;
j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);

ezi = fac.*ci.*ez;
ezr = fac.*ez;

fxx0 = (ezi.*pc+ezr.*ps)*j0; % cos(0*phi)-component
fxx2 = -(ezi.*pc-ezr.*ps)*j2; % cos(2*phi)-component
fxz =  -2*1i*(ezi.*v)*j1; % cos(1*phi)-component

byx0 = (ezr.*pc+ezi.*ps)*j0; % cos(0*phi)-component
byx2 = -(ezr.*pc-ezi.*ps)*j2; % cos(2*phi)-component
byz = -2*1i*(ezr.*v)*j1; % cos(1*phi)-component

phi = 0:pi/100:2*pi; phi = phi'; col = ones(size(phi));
if nargin>14 && ~isempty(orient) % image of a molecule with orientation orient
    phi = phi - orient(2);
    inty = real((sin(orient(1))*(col*fxx0+cos(2*phi)*fxx2)+cos(orient(1))*cos(phi)*fxz).*...
        conj((sin(orient(1))*(col*byx0+cos(2*phi)*byx2))+cos(orient(1))*cos(phi)*byz)); 
    intz = real((sin(orient(1))*sin(2*phi)*fxx2+cos(orient(1))*sin(phi)*fxz).*...
        conj(sin(orient(1))*sin(2*phi)*byx2+cos(orient(1))*sin(phi)*byz)); 
    phi = phi + orient(2);
    intx = inty + intz;    
else
    intx = real((col*fxx0+cos(2*phi)*fxx2).*conj(col*byx0+cos(2*phi)*byx2)); % int for x-dipole along x-pol
    inty = real((sin(2*phi)*fxx2).*conj(sin(2*phi)*byx2)); % int for x-dipole along y-pol
    intz = col*real(fxz.*conj(byz)); % int for z-dipole along x-pol
end
