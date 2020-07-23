function [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm, imm_c] = FindPatternSep(im,mask,sup,bck,sze,tsh,fun,flag)
% FindPatternSep takes the intensity image and the computed model patterns
% to identify single molecules. 
%---
% Input variables:
% im – measured image mask – 3D stack of model patterns 
% sup – support of pattern matching, see sjk above; usually a square or circular disk
% bck – background bjk, usually the same as sup 
% sze – a threshold parameter which can be used to make the pattern recognition 
% more discriminating for suppressing false posi- tives (default 1)
% tsh – threshold value ?, see above (default 1)
% fun – string with function for deciding a match, default is 'cim>tsh*sqrt(err)'
% flag – if flag is 1, then pattern matching is suppressed ate the image borders 
% within a stripe of half the width of the model images; this is important to 
% prevent false recognitions in the border region (default 1).
%---
% Output variables:
% err – error image 
% xc, yc – center coordinates of found patterns
% sc – corresponding pattern index 
% len – number of found patterns 
% imm - reconstructed image with found model patterns with each pattern in
% a different matrix plane.
% imm_c – reconstructed image with found model patterns


if nargin<7 || isempty(fun)
    
    
    fun = 'cim>tsh*sqrt(err)';
end
n1 = size(mask,1);
n2 = size(mask,2);
n3 = size(mask,3);
if nargin<4 || isempty(bck)
    bck = ones(n1,n2,n3);
elseif size(bck,3)<size(mask,3)
    bck = repmat(bck,[1 1 n3]);
end
if nargin<3 || isempty(sup)
    sup = ones(n1,n2,n3);
elseif size(sup,3)<size(mask,3)
    sup = repmat(sup,[1 1 n3]);
end
if nargin<8 || isempty(flag)
    flag = true;
end
for j=1:n3
    bck(:,:,j) = bck(:,:,j)/sqrt(sum(sum(sup(:,:,j).*bck(:,:,j).^2)));
    mask(:,:,j) = mask(:,:,j).*(sup(:,:,j)>0);
    mask(:,:,j) = mask(:,:,j)./sqrt(sum(sum(sup(:,:,j).*mask(:,:,j).^2)));
end
err = inf*im;
cim = im;
bim = im;
sim = 1+0*im;
for s = 1:n3
    im0 = mConv2(im,squeeze(sup(:,:,s).*bck(:,:,s)));
    im02 = mConv2(im.^2,squeeze(sup(:,:,s)));
    crs = sum(sum(sup(:,:,s).*bck(:,:,s).*mask(:,:,s)));
    crs = inv([1 crs; crs 1]);
    im1 = mConv2(im,squeeze(sup(:,:,s).*mask(:,:,s)));
    tmperr = im02 - crs(1,1)*im0.*im0 - 2*crs(1,2)*im0.*im1 - crs(2,2)*im1.*im1;
    tmpbim = crs(1,1)*im0 + crs(1,2)*im1;
    tmpcim = crs(1,2)*im0 + crs(2,2)*im1;
    cim(tmperr<err) = tmpcim(tmperr<err);
    bim(tmperr<err) = tmpbim(tmperr<err);
    sim(tmperr<err) = s;
    err(tmperr<err) = tmperr(tmperr<err);
end

if nargout>4
    if nargin<5 || isempty(sze)
        sze = 1;
    end
    if nargin<6 || isempty(tsh)
        tsh = 1;
    end
    eval(['[cl,len]=mCluster(' fun ',sze);'])
    [a,b] = meshgrid(1:size(im,2),1:size(im,1));
    cnt = 1; xc = []; yc = []; bc = []; cc = []; sc = [];
    for j=1:length(len)
        ind = cl==j;
        tmpxc = a(ind);
        tmpyc = b(ind);
        tmperr = err(ind);
        tmpb = bim(ind);
        tmpc = cim(ind);
        tmps = sim(ind);
        ind = tmperr==min(tmperr);
        sc(cnt) = mean(tmps(ind));
        bc(cnt) = mean(tmpb(ind));
        cc(cnt) = mean(tmpc(ind));
        xc(cnt) = mean(tmpxc(ind));
        yc(cnt) = mean(tmpyc(ind));    
        cnt = cnt+1;
    end
    if ~isempty(xc)
        xc = round(xc);
        yc = round(yc);
        if flag
            ind = xc<(n2+1)/2 | xc>size(im,2)-(n2-1)/2 | yc<(n1+1)/2 | yc>size(im,1)-(n1-1)/2; 
            xc(ind)=[]; yc(ind)=[]; bc(ind) = []; cc(ind) = []; sc(ind) = []; len(ind) = [];
        end
    end
end

if nargout>10 
    if ~isempty(xc)
        imm = zeros(size(im,1), size(im,2), numel(xc));
        [a,b] = size(im);
        nn = (n1-1)/2;
        for k=1:length(xc)
%             imm(round(max(yc(k)-nn,1):min(yc(k)+nn,a)),round(max(xc(k)-nn,1):min(xc(k)+nn,b))) = imm(round(max(yc(k)-nn,1):min(yc(k)+nn,a)),round(max(xc(k)-nn,1):min(xc(k)+nn,b))) + ...
%                 cim(yc(k),xc(k))*mask(round(max(nn-yc(k)+2,1):min(2*nn+1-yc(k)-nn+a,2*nn+1)),round(max(nn-xc(k)+2,1):min(2*nn+1-xc(k)-nn+b,2*nn+1)),sim(yc(k),xc(k)));
%         end
             imm(round(max(yc(k)-nn,1):min(yc(k)+nn,a)),round(max(xc(k)-nn,1):min(xc(k)+nn,b)),k)=cim(yc(k),xc(k))*mask(round(max(nn-yc(k)+2,1):min(2*nn+1-yc(k)-nn+a,2*nn+1)),round(max(nn-xc(k)+2,1):min(2*nn+1-xc(k)-nn+b,2*nn+1)),sim(yc(k),xc(k)));
%         
        end
    else
        imm = [];
    end
end
imm_c=sum(imm,3);
end