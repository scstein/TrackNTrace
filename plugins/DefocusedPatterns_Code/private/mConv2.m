function y = mConv2(a, b)

% y = mConv2(a, b) convolves matrix a with b where it is assumed that b has smaller size than a.

[ma,na,bla] = size(a);
[mb,nb,bla] = size(b);
if rem(mb,2)==0
   mb = mb/2;
   mbb = 0;
else
   mb = (mb-1)/2;
   mbb = 1;
end
if rem(nb,2)==0
   nb = nb/2;
   nbb = 0;
else
   nb = (nb-1)/2;
   nbb = 1;
end
% fa = fft2([a(mb:-1:1,nb:-1:1) a(mb:-1:1,:) a(mb:-1:1,end:-1:end-nb+1); ...
%       a(:,nb:-1:1) a a(:,end:-1:end-nb+1); ...
%       a(end:-1:end-mb+1,nb:-1:1) a(end:-1:end-mb+1,:) a(end:-1:end-mb+1,end:-1:end-nb+1)]);
fa = fft2(cat(1,cat(2,a(mb:-1:1,nb:-1:1,:),a(mb:-1:1,:,:),a(mb:-1:1,end:-1:end-nb+1,:)), ...
    cat(2,a(:,nb:-1:1,:),a,a(:,end:-1:end-nb+1,:)), ...
    cat(2,a(end:-1:end-mb+1,nb:-1:1,:),a(end:-1:end-mb+1,:,:),a(end:-1:end-mb+1,end:-1:end-nb+1,:))));
mask = zeros(size(fa));
mask(1:mb+mbb,1:nb+nbb,:) = b(mb+1:end,nb+1:end,:);
mask(end-mb+1:end,1:nb+nbb,:) = b(1:mb,nb+1:end,:);
mask(1:mb+mbb,end-nb+1:end,:) = b(mb+1:end,1:nb,:);
mask(end-mb+1:end,end-nb+1:end,:) = b(1:mb,1:nb,:);
mask = fft2(mask);
tmp = real(ifft2(fa.*conj(mask)));
y = tmp(mb+1:end-mb,nb+1:end-nb,:);

