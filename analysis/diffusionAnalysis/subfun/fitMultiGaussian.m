function [F,J] = fitMultiGaussian(p,x,y,use_v)

if nargin ~=4
    error('prog:input','Wrong number of input arguments, aborting. \n');
end

if use_v
    order = numel(p)/2;
else
    order = (numel(p)+1)/2;
end

F = zeros(size(x,1),1);
J = [];

if use_v
    if order == 1
        F = F+1/sqrt(2*pi*p(1))*exp(-(x-p(2)).^2/(2*p(1)));
    else
        for i = 1:order-1
            F = F+abs(p(i))/sqrt(2*pi*p(order-1+i))*exp(-(x-p(end)).^2/(2*p(order-1+i)));
        end
        F = F+abs(1-sum(abs(p(1:order-1))))/sqrt(2*pi*p(2*order-1))*exp(-(x-p(end)).^2/(2*p(2*order-1)));
    end
else
    if order == 1
        F = F+1/sqrt(2*pi*p(1))*exp(-x.^2/(2*p(1)));
    else
        for i = 1:order-1
            F = F+abs(p(i))/sqrt(2*pi*p(order-1+i))*exp(-x.^2/(2*p(order-1+i)));
        end
        F = F+abs(1-sum(abs(p(1:order-1))))/sqrt(2*pi*p(2*order-1))*exp(-x.^2/(2*p(2*order-1)));
    end
end

F = y-F;