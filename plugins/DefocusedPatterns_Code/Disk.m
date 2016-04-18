function y = Disk(m,opt)

if length(m)==1
	m = [m m];
end

[jx,jy] = meshgrid(-m(2):m(2),-m(1):m(1));

if nargin<2
    y = ceil(sqrt(jx.^2/m(1)^2+jy.^2/m(1)^2))<2;
else
    rad = sqrt(jx.^2 + jy.^2);
    eval(['y = Disk([' num2str(m(1)) ',' num2str(m(2)) ']).*' opt ';']);
end
y = y/sum(sum(y));

