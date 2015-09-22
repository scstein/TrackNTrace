function [p,S,mu] = polyfit3(x,y,n,nul,w)
%POLYFIT3 Fit polynomial to data, with new features.
%
% The new usage:
%   POLYFIT3(X,Y,N,NUL,W) finds the coefficients of a polynomial P(X) of
%   degree N that fits the data, P(X(I))~=Y(I), in a least-squares sense.
%   Any of the coefficients can be forced to be zero, and data can be weighted.
%
%   NUL is a vector with those coefficients forced to be zero in usual order.
%   W is a vector containing the weights of data.
%
%  POLYFIT3(X,Y,N,NUL) assumes that all weights of data are equal to 1.
%  POLYFIT3(X,Y,N,[],W) allows to weight data without forcing any coefficient to be zero.
%
% The original usage is still available:
%   POLYFIT3(X,Y,N) finds the coefficients of a polynomial P(X) of
%   degree N that fits the data, P(X(I))~=Y(I), in a least-squares sense.
%
%   The structure S contains the Cholesky factor of the Vandermonde
%   matrix (R), the degrees of freedom (df), and the norm of the
%   residuals (normr) as fields. Note that df also depends on NUL vector.
%   Residuals are weighted by sqrt(w/mean(w)).
%
%   [P,S,MU] = POLYFIT(X,Y,N,W) finds the coefficients of a polynomial
%   in XHAT = (X-MU(1))/MU(2) where MU(1) = mean(X) and MU(2) = std(X).
%   This centering and scaling transformation improves the numerical
%   properties of both the polynomial and the fitting algorithm.
%
%   Warning messages result if N is >= length(X), if X has repeated, or
%   nearly repeated, points, or if X might need centering and scaling.
%
%   See also POLYFIT, POLY, POLYVAL, ROOTS.
%   New features added by Antoni J. Canos, 12-15-03. Tested under MATLAB 6.0.0. R12.

% The regression problem is formulated in matrix format as:
%
%    W*y = W*V*p    or
%
%              3  2
%    W*y = W*[x  x  x  1] [p3
%                          p2
%                          p1
%                          p0]
%
% where the vector p contains the coefficients to be found.  For a
% 7th order polynomial, matrix V would be:
%
% V = [x.^7 x.^6 x.^5 x.^4 x.^3 x.^2 x ones(size(x))];
%

% DEMO (Not weighted data).
%x=linspace(0,11,12);
%noise=(rand(size(x))-0.5);
%y=3-x.^3+2*x.^5+noise;
%[p1,S] = polyfit3(x,y,5)
%w=ones(size(x));
%[p2,S] = polyfit3(x,y,5,[2 4 5],w)
%close all;
%figure;
%plot(x,y,'.r'); hold on;
%x=linspace(0,11);
%plot(x,polyval(p1,x),'b');
%plot(x,polyval(p2,x),'k');

% DEMO (Weighted data).
%x=linspace(-5,5,11);
%noise=(rand(size(x))-0.5);
%y=x.^4+noise; y(11)=y(11)+100;
%[p1,S] = polyfit3(x,y,4)
%w=ones(size(x)); w(11)=0.1;
%[p2,S] = polyfit3(x,y,4,[2 3 4 5],w)
%close all;
%figure;
%plot(x,y,'.r'); hold on;
%x=linspace(-5,5,11);
%plot(x,polyval(p1,x),'b');
%plot(x,polyval(p2,x),'k');

if ~isequal(size(x),size(y))
    error('X and Y vectors must be the same size.')
end

x = x(:);
y = y(:);

if nargin == 3
    nul=[];
    w=ones(size(x));
elseif nargin == 4
    w=ones(size(x));
end

%Normalized weights
w=w(:)./mean(w);

if nargout > 2
   mu = [mean(x); std(x)];
   x = (x - mu(1))/mu(2);
end

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1);
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

% To be used in polyval
Vfull=V;

V(:,nul)=[];

% Construct Weighted Vandermonde matrix.
W=diag(w);
WV=W*V;
Wy=W*y;

% To be used in polyval
WVfull=W*Vfull;

% Solve least squares problem, and save the Cholesky factor.
[Q,R] = qr(WV,0);
ws = warning('off'); 
p = R\(Q'*Wy);    % Same as p = WV\Wy;
warning(ws);
if size(R,2) > size(R,1)
   warning('Polynomial is not unique; degree >= number of data points.')
elseif condest(R) > 1.0e10
    if nargout > 2
        warning(sprintf( ...
            ['Polynomial is badly conditioned. Remove repeated data points.']))
    else
        warning(sprintf( ...
            ['Polynomial is badly conditioned. Remove repeated data points\n' ...
                '         or try centering and scaling as described in HELP POLYFIT.']))
    end
end
r = sqrt(w).*(y - V*p);
p = p.';         % Polynomial coefficients are row vectors by convention.
 
% To be used in polyval
p(setdiff([1:n+1],nul))=p;
p(nul)=0;
[Q,R] = qr(WVfull,0);

% S is a structure containing three elements: the Cholesky factor of the
% Vandermonde matrix, the degrees of freedom and the norm of the residuals.

S.R = R;
S.df = length(y) - (n+1)+length(nul);
S.normr = norm(r);
