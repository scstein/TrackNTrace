function [F,J] = fitMultiRGaussian(p,x,y,~)

if nargin ~=4
    error('prog:input','\nWrong number of input arguments, aborting.');
end

order = (numel(p)+1)/2;

F = zeros(size(x,1),1);
J = [];


if order == 1
    F = 1/p(1)*x.*exp(-x.^2/(2*p(1)));
else
    for i = 1:order-1
        F = F+abs(p(i))/p(order-1+i)*x.*exp(-x.^2/(2*p(order-1+i)));
    end
    F = F+(1-sum(abs(p(1:order-1))))/p(end)*x.*exp(-x.^2/(2*p(end)));
end


F = y-F;