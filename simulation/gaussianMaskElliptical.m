function [img] = gaussianMaskElliptical(x,y,amp,sigma_x,sigma_y,angle,halfw,px_integrated)

[x_grid,y_grid] = meshgrid (-halfw:halfw,-halfw:halfw);

a = cos(angle)^2/(2*sigma_x^2)+sin(angle)^2/(2*sigma_y^2);
b = -sin(2*angle)/(4*sigma_x^2)+sin(2*angle)/(4*sigma_y^2);
c = sin(angle)^2/(2*sigma_x^2)+cos(angle)^2/(2*sigma_y^2);

if ~px_integrated || abs(angle)>1e-6
    img = amp*exp(-(x_grid-x).^2*a+2*b*(x_grid-x).*(y_grid-y)-c*(y_grid-x).^2);
else
    img = amp*pi*sigma_x*sigma_y/2*(erfc(-(x_grid-x+0.5)/(sqrt(2)*sigma_x))-erfc(-(x_grid-x-0.5)/(sqrt(2)*sigma_x))).*...
        (erfc(-(y_grid-y+0.5)/(sqrt(2)*sigma_y))-erfc(-(y_grid-y-0.5)/(sqrt(2)*sigma_y)));
end