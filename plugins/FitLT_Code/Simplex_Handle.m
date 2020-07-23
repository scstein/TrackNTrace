function [x, dx, steps, f] = Simplex_Handle(fhandle, x, xmin, xmax, tol, steps, varargin)

%	[x, dx, steps, f] = Simplex('F', X0, XMIN, XMAX, TOL, STEPS, VARARGIN) 
%	attempts to return a vector x and its error dx, so that x minimzes the 
%	function f = F(x) near the starting vector X0 under the conditions that 
% 	xmin <= x <= xmax.
%	TOL is the relative termination tolerance dF/F; (default = 1e-10)
%	STEPS is the maximum number of steps; (default = 200*number of parameters).
%	The returned value of STEPS is the actual number of performed steps.
%   Pass a negative value to silence the output. In any case abs(steps) is used.
%	Simplex allows for up to 10 additional arguments for the function F.
%	Simplex uses a Nelder-Mead simplex search method.
%
%   Modified by CT, 2017:
%   Parameters and constrains now can be NaN, to indicate that they should be fixed and ignored. 
%   Modified by CT, 2018:
%   The fitting is now skipped if the error of the initial values is 0.
%   (Happens for empty TCSPC).
%   Modified by CT, 2019:
%   Enforce constrains on reflection and expansion.
%   Modified by CT, 2020:
%   A negtive number of steps silences the output.
%
%#function TCSPC_Fun IRF_Fun 

x = x(:);
if nargin<5
    tol = 1e-10;
    if nargin<4
        xmax = Inf*ones(length(x),1);
        if nargin<3
            xmin = -Inf*ones(length(x),1);
        end
    end
elseif isempty(tol)
    tol = 1e-5;
end
if nargin<6
    steps = [];
end
if isempty(xmin)
    xmin = -Inf*ones(size(x));
end
if isempty(xmax)
    xmax = Inf*ones(size(x));
end
xmin = xmin(:);
xmax = xmax(:);
xmax(xmax<xmin) = xmin(xmax<xmin);
x(x<xmin) = xmin(x<xmin);
x(x>xmax) = xmax(x>xmax);
xfix = zeros(size(x));
tmp = xmin==xmax | isnan(xmin) & isnan(xmax);
xfix(tmp) = xmin(tmp);
mask = diag(~tmp);
mask(:, tmp) = [];
x(tmp) = [];
xmin(tmp) = [];
xmax(tmp) = [];

% if isa(fhandle,'function_handle')
%     fun = fhandle;
%     evalstr = 'fun';
% else
%     evalstr = fhandle;
% end
%fv = feval(fhandle,mask*x+xfix,extraArg{:});
% evalstr = [evalstr, '(mask*x+xfix'];
if nargin>6
%     evalstr = [evalstr, ',varargin{:}'];
    extraArg = varargin;
else
    extraArg = {};
end
% evalstr = [evalstr, ')'];

n = length(x);
if n==0 
	x = xfix;
	dx = zeros(size(xfix));
	steps = 0;
	return
end
if isempty(steps)
	steps = 200*n;
end

xin = x(:);
%v = 0.9*xin;
v = xin;
x(:) = min(xmax,max(xmin,v)); 
%fv = eval(evalstr); 
fv = feval(fhandle,mask*x+xfix,extraArg{:});
% Exit if error is 0. Otherwise the fitting loop becomes infinte.
if fv==0
    steps = 0;
    f = fv;
    return
end
for j = 1:n
	y = xin;
    if y(j) ~= 0
        y(j) = (1 +.2*rand)*y(j);
    else
        y(j) = 0.2;
    end
    if y(j)>=xmax(j)
        y(j) = xmax(j);
    end
    if y(j)<=xmin(j)
        y(j) = xmin(j);
    end
    v = [v y];
	x(:) = y; 
    %f = eval(evalstr);
    f = feval(fhandle,mask*x+xfix,extraArg{:});
	fv = [fv f];
end
[fv, j] = sort(fv);
v = v(:,j);
count = n+1;

% Parameter settings for Nelder-Meade
alpha = 1; beta = 1/2; gamma = 2;

% Begin of Nelder-Meade simplex algorithm
% 
while count < abs(steps)
	if (2*abs(fv(n+1)-fv(1))/(abs(fv(1))+abs(fv(n+1))) <= tol)
		break
	end

	% Reflection:
	vmean = sum(v(:, 1:n),2)/n; % Same as mean but 4 times faster
	vr = (1 + alpha)*vmean - alpha*v(:, n+1);
	x(:) = min(xmax,max(xmin,vr));                  % It might be better to constrain vr, but this version is consitent with the previous estimation of the parameter contribution.
    fr = feval(fhandle,mask*x+xfix,extraArg{:});
	count = count + 1; 
	vk = vr; fk = fr;

	if fr < fv(1) && all(xmin<=vr) && all(vr<=xmax)
		% Expansion:
		ve = gamma*vr + (1-gamma)*vmean;
		x(:) = min(xmax,max(xmin,ve));
        fe = feval(fhandle,mask*x+xfix,extraArg{:});
		count = count + 1;
		if fe < fv(1) && all(xmin<=ve) && all(ve<=xmax)
			vk = ve; fk = fe;
		end
	else
		vtmp = v(:,n+1); ftmp = fv(n+1);
		if fr < ftmp && all(xmin<=vr) && all(vr<=xmax)
			vtmp = vr; ftmp = fr;
		end
		% Contraction:
		vc = beta*vtmp + (1-beta)*vmean;
		x(:) = min(xmax,max(xmin,vc));
        fc = feval(fhandle,mask*x+xfix,extraArg{:});
		count = count + 1;
		if fc < fv(n) && all(xmin<=vc) && all(vc<=xmax)
			vk = vc; fk = fc;
		else
			% Shrinkage:
			for j = 2:n
				v(:, j) = (v(:, 1) + v(:, j))/2;
				x(:) = v(:, j);
                fv(j) = feval(fhandle,mask*x+xfix,extraArg{:});
			end
			count = count + n-1;
			vk = (v(:, 1) + v(:, n+1))/2;
			x(:) = vk;
            fk = feval(fhandle,mask*x+xfix,extraArg{:});
			count = count + 1;
		end
	end
	v(:, n+1) = vk;
	fv(n+1) = fk;
	[fv, j] = sort(fv);
	v = v(:,j);
end

x = v(:,1);
dx = abs(v(:,n+1)-v(:,1));
x = mask*x + xfix;
dx = mask*dx;
f = fv(1);
if steps>0 && count>=steps
	disp(['Warning: Maximum number of iterations (', int2str(steps),') has been exceeded']);
else
	steps = count;
end
