function [ y,x ] = ZeroPaddingFFT( corr, subboxsize, paddingfactor )

%crop out the correlation peak
boxsize = size(corr,1);
sy = boxsize/2;
sx = sy;
offsetx = floor(sx - subboxsize/2);
offsety = floor(sy - subboxsize/2);
if offsetx<0
    offsetx=0;
elseif offsetx>boxsize
    offsetx=boxsize-subboxsize;
end
if offsety<0
    offsety=0;
elseif offsety>boxsize
    offsety=boxsize-subboxsize;
end
subbox=corr(offsety+1:offsety+subboxsize,offsetx+1:offsetx+subboxsize);

%Fourier interpolation, zero-padding
data = fftshift(fft2(subbox));

padsize = subboxsize * paddingfactor;
padoffset = padsize/2-subboxsize/2;
padbox = zeros(padsize,padsize);
padbox(padoffset+1:padoffset+subboxsize,padoffset+1:padoffset+subboxsize)=data;
data=abs(ifft2(padbox));

[y1,x1]=find(data == max(max(data)));
y = mean(y1)/paddingfactor+offsety;
x = mean(x1)/paddingfactor+offsetx;

end

