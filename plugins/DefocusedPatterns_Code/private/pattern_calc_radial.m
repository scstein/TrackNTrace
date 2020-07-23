pattern = 'radial';
pixel = 0.05;
be_res = 5; % minimum resolution of in-plane angle
al_res = 5; % minimum resolution of out-of-plane angle
cnt = 1;
for k=90:-al_res:0
    al = k/180*pi;
    if k==90
        jj = round(180/be_res);
        dbe = pi/jj;
    elseif k==0
        jj = 1;
        dbe = 0;
    else
        jj = round(sin(al)*360/be_res);
        dbe = 2*pi/jj;
    end
    for j=1:jj
        theta(cnt) = al;
        phi(cnt) = dbe*(j-1);
        cnt = cnt+1;
    end
end


NA = 1.49;
n0 = 1.52;
n = 1;
n1 = 1;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;
focpos = 0;
pic = 1;

[err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm, imm_c, mask, molim] = radialpatterns_focus(tag, pattern, pixel, theta, phi, NA, n0, n, n1, d0, d, d1, lamex, focpos, pic);