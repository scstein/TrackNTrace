function f=gaussian2d(N,sigma, pixIntegrated)
%Usage   g_img = gaussian2d(N,sigma, pixIntegrated)
% N: grid size
% sigma: standard deviation
% pixIntegrated: Should the gaussian be integrated over the pixel area?
%
% The returned gaussian is normalized.

if nargin<3
    pixIntegrated = false;
end

[x, y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));

if(pixIntegrated)
    f = pi*sigma^2/2 * ( erfc(-(x-0+1/2)/(sqrt(2)*sigma)) - erfc(-(x-0-1/2)/(sqrt(2)*sigma))    )...
                     .* ( erfc(-(y-0+1/2)/(sqrt(2)*sigma)) - erfc(-(y-0-1/2)/(sqrt(2)*sigma))   ) ;
else
    f= exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
end

f=f./sum(f(:));