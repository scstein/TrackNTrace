function [v,pc,ps,tp,ts,tp1,ts1,fac] = DipoleL(theta,z,n0,n,n1,d0,d,d1)

% [v,pc,ps] = DipoleL(theta,z,n,d) calculates the electric field amplitudes of dipole radiation along emission angle theta 
% of a dipole at distance z from an interface within a layer 
% theta - direction of radiation downwards
% z  - molecule's distance from the bottom of its layer
% n0 - vector of refractive indices of the stack below the molecule's layer
% n  - refracive index of the molecule's layer
% n1 - vector of refractive indices of the stack above the the molecule's layer
% d0 - vector of layer thickness values of the stack below the molecule's layer ( length(d0)=length(n0)-1 )
% d  - thickness of molecule's layer
% d1 - vector of layer thickness values of the stack above the molecule's layer ( length(d1)=length(n1)-1 )

z = z(:)';
col = ones(size(z));
theta = abs(theta(:));
ind = theta<=pi/2;

v = zeros(length(theta),length(z));
pc = v; ps = v;

if sum(ind)>0
    tmp = theta(ind);
    w = sqrt(n^2-n0(1)^2*sin(tmp).^2);
    [rpu, rsu, tpu, tsu] = Fresnel(w,[n n1],d1);
    [rpd, rsd, tpd, tsd] = Fresnel(w,[n n0(end:-1:1)],d0(end:-1:1));    
    tp = (tpd./(1-rpu.*rpd.*exp(2*1i*w.'*d))).';
    ts = (tsd./(1-rsu.*rsd.*exp(2*1i*w.'*d))).';
    tp1 = tp.*(rpu.*exp(2*1i*w.'*d)).';
    ts1 = ts.*(rsu.*exp(2*1i*w.'*d)).';
    fac = sqrt(n0(1))*n0(1)*cos(tmp)./w;
    ez = exp(1i*w*z);
    % all components in a occordinate system with z-axis pointing away from surface
    v(ind,:) = ((fac.*n0(1)/n.*sin(tmp).*tp)*col).*ez + ((fac.*n0(1)/n.*sin(tmp).*tp1)*col)./ez;
    pc(ind,:) = ((fac.*w/n.*tp)*col).*ez - ((fac.*w/n.*tp1)*col)./ez;
    ps(ind,:) = ((fac.*ts)*col).*ez + ((fac.*ts1)*col)./ez;
end
