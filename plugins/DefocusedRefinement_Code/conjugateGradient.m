function [ p, fret ] = conjugateGradient( p0, fname, gname, FTOL, GSTEP, ITMAX, ITFUNC, verbose, varargin)
% Usage:  [p, fret] = conjugateGradient( p0, fname, gname, FTOL, ITMAX, ITFUNC,  a,b,[..])
% Conjugate gradient to minimize the function f(p, a,b,[...]) given by
% the function handle or name 'fname' with gradient 'gname' starting at the
% initial guess p0. Returns the position of the minimum 'p' and the
% corresponding value of the function at the minimum 'fret'.
%
% The function is minimized in the variables specified in
% the parameter vector p, while a,b, etc. are fixed paramter values. If the
% gradient is not supplied(gname==[]) a numerical one will be used.
%
%   Minimal Example - Find Minimum of 2d Parabola:
%   f  = @(x,a,b) a*(x(1)-1).^2 + (x(2)-b).^2; % function to optimize
%   a = 100;  b = 20;  % Set values of fixed paramters
%   x0 = [2,10];       % initial guess for min. position
%   ITMAX = 100;       % Maximum number of iterations
%   [x_min, fret] = conjugateGradient(x0, f, [], [], ITMAX, [], a,b);
%
% Input:
%    p0     -  Variable vector of initial values.
%    fname  -  Name or handle of function to optimize. The function must be
%              callable as fval = f(p,  a,b,[...]), where p is the vector
%              of variables you want to optimize in and a,b, etc. are fixed parameters.
%    gname  -  Gradient of function 'fname'. Must be callable as grad = g(p, a,b,[...]),
%              where grad is the ROW(!) vector matching the variables in p.
%              If gname is empty ('[]') a numerical gradient will be used.
%    FTOL   -  Function tolerance for minimum search. (default 1e-8). Giving
%              negative values forces search to continue up to ITMAX iterations.
%    GSTEP  -  Gradient step for bracketing (finding minima along the gradient line)
%              For too low values the minimization takes longer or might
%              not work well. For too high values the parameters might
%              explode.              
%    ITMAX  -  Maximum number of iterations to run. (default 100);
%    ITFUNC -  (OPTIONAL) Function handle or name callable as pNew = ITFUNC(p, a,b,[...])
%              that is called after every iteration. Can be used to visualize or alter
%              the result after every step (e.g. projection on some set of solutions).
%    verbose - boolean. If true (default) printf progress to the MATLAB console.
%    varargin: Comma seperated parameter list (a,b,[...]) neccessary for the
%              evaluation of function 'fname' and gradient 'gname'.
% Output:
%    p      -  Position (variable values) of the minimum of the function.
%    fret   -  Function value at the minimum

% Adapted and extended C Code of numerical recipes 3rd edition
% Author: Simon Christoph Stein
% Date: 2015

if nargin < 3 || isempty(gname)
    % If the gradient is not supplied, we use a numerical gradient of the
    % function for evaluation.
    gname = @(x,varargin) NUMERIC_GRAD(fname, x, varargin{:});
end

if nargin < 5 || isempty(FTOL)
    % Decreasing the tolarance increases speed
    % Minimum can be about 1e-8 for double, 1e-4 for single precision
    % (square root of machine precision)
    FTOL = 1e-8; % fractional tolerance to stop optimizing
end

if nargin < 6 || isempty(GSTEP)
    GSTEP = 1e-12; % gradient step for bracketing    
end

if nargin < 7 || isempty(ITMAX)
    ITMAX = 100; % max iterations
end

if nargin < 8 || isempty(ITFUNC)
    ITFUNC = @(x, varargin) x; % Function called after each iteration
end

if nargin < 9 || isempty(verbose)
    verbose = true;
end

EPSILON = 1e-18;
GTOL = 1e-18; % gradient tolerance


% Maximum number of restarts with smaller GSTEP
% Optimization is restartet in case of errors during function evaluation.
MAX_RETRY = 20; 
Retry_cnt = 0;

while true
    try
        p = p0(:);
        % compute function value
        fp = feval(fname, p, varargin{:});
        
        % compute gradient value
        g = -feval(gname, p, varargin{:});
        h = g;
        xi = g;
        
        reverseStr = ''; %% for one-line printing
        for its = 1:ITMAX
            if verbose
            % Command line message
                msg = sprintf('CJG: iteration %i/%i - fval: %f',its, ITMAX, fp);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            
            % Map function to 1D  for line search along search direction xi
            func1dim = @(alpha) feval(fname, p + alpha*xi, varargin{:});
            
            % Bracket the minimum
            % Inital points  [x] & [x + GSTEP*g]
            [a,b,c] = bracket(0, GSTEP, func1dim);
            
            
            % perform minimum search
            %   [alphaMin, fmin] = goldenSectionSearch(a,b,c, func1dim, 1e-3);
            [alphaMin, fret] = brentsMinimumSearch(a,b,c, func1dim, 1e-3);
            
            p = p + alphaMin*xi; % set p to found minimum
            
            % Invoke iteration function
            p = feval(ITFUNC, p, varargin{:});
            
            
            % One possible return
            if(2.0*abs(fret-fp) <= FTOL * (abs(fret) + abs(fp) + EPSILON) )
                if verbose
                    fprintf('\nCJG succeeded: Found minimum with respect to tolerance.\n')
                end
                return
            end
            
            fp = fret;
            xi = feval(gname, p,  varargin{:});
            den = max(fp,1.0);
            temp = abs(xi) .* max( [abs(p), ones(size(p,1),1)] , [], 2)./den ;
            test = max(temp);
            % The other possible return
            if (test < GTOL)
                if verbose
                    fprintf('\nCJG succeeded: Gradient small.\n')
                end
                return
            end
            gg = sum(g.^2);
            %     dgg = sum(xi.^2)   % Fletcher Reeves
            dgg = sum((xi + g).*xi);  % Polak Ribiere
            
            if gg == 0.0 %unlikely gradient exactly zero
                if verbose
                    fprintf('\nCJG succeeded: gradient zero.\n')
                end
                return
            end
            
            gam = dgg/gg;
            g = -xi;
            h = g + gam*h;
            xi = h;
        end
        
        break;
    catch
        if Retry_cnt == MAX_RETRY
           error('CJG: Max number of retries due to error reached. Aborting optimization.') ;
        end
        GSTEP = GSTEP/10;
        warning('CJG: An error occured during function evaluation, retrying with smaller GSTEP %g..\n',GSTEP);
        Retry_cnt = Retry_cnt+1;
    end
    
    if Retry_cnt>0
       warning('WARNING: GSTEP was reduced during run time due to errors during function evaluation. This forces the algorithm to restart the complete optimization. Choosing a smaller initial value may prevent this.') 
    end
    
    if verbose
        fprintf('\n');
    end
    % warning('CJG: Max iterations reached in conjugateGradient.');
    
end

end


function g = NUMERIC_GRAD(fname, x, varargin)
% Calculates the numerical (central) gradient of function with name/handle
% 'fname' [f = feval(fname, x, pars)] in the variables x with stepsize h.
% The struct pars must hold all parameters neccessary for evaluating
% 'fname' (e.g. pars.a, pars.b).

% minimum step size
H_MIN = 1e-8;

g = zeros(numel(x), 1);
for gradi = 1:numel(x)
    % Set step size
    h = max(H_MIN, abs(x(gradi))* 1e-4);
    
    %Buffer value of the gradi-th variable
    xi_tmp = x(gradi);
    
    % Compute f(x+h)
    x(gradi) = xi_tmp+h;
    f1 = feval(fname, x, varargin{:});
    
    % Compute f(x-h)
    x(gradi) = xi_tmp-h;
    f2 = feval(fname, x, varargin{:});
    
    % Central derivative
    g(gradi) = (f1-f2)/(2*h);
    
    %Reset value of the gradi-th variable
    x(gradi) = xi_tmp;
end

end



function [ax,bx,cx] = bracket(a,b, fname)
% Bracket the minimum of a 1d-function fname using the inital points a,b
% Given the output triple ax<bx<cx, the minimum is within (ax,cx)
% and f(b) < f(a), f(b) < f(c) holds.

% Parabolic extrapolation to the minimum is used for bracketing
% -> For parabolic functions this will already return the minimum as the
% midpoint.

% From numerical recipes 3rd edition
% Author: Simon Christoph Stein
% Date: 2015

GOLD = 1.618034; % Magnification from golden ratios
TINY = 1e-20; % Prevent division by zero
GLIMIT = 100; % Maximum allowed magnification per step.

ax = a;
bx = b;

fa = feval(fname,ax);
fb = feval(fname,bx);

% Swap so that a-> b is downhill
if (fb > fa)
    [ax,bx] = swap(ax,bx);
    [fa,fb] = swap(fa,fb);
end

% First guess for c
cx = bx + GOLD*(bx-ax);
fc = feval(fname, cx);

while(fb > fc) % For fc>fb we have bracketet the minimum
    % Parabolic extrapolation using a,b,c. Search for minimum and save in u
    r = (bx-ax)*(fb-fc);
    q = (bx-cx)*(fb-fa);
    if (q-r) < 0
        u = bx - ((bx-cx)*q-(bx-ax)*r) / (2*min(q-r,-TINY));
    else
        u = bx - ((bx-cx)*q-(bx-ax)*r) / (2*max(q-r,TINY));
    end
    ulim = bx + GLIMIT*(cx-bx); % max allowed value for u
    
    
    % Check various cases
    
    % Case 1:  u is between b and c
    if( (bx-u) * (u-cx) > 0 )
        fu = feval(fname, u);
        if( fu < fc) % Minimum between b and c!
            ax = bx;
            bx = u;
            fa = fb;
            fb = fu;
            return
        elseif (fu > fb) % Minimum between a and u (fb<fa always holds)
            cx = u;
            fc = fu;
            return
        end
        u = cx + GOLD*(cx-bx); % parabolic guess did not find minimum, use default magnification
        fu = feval(fname, u);
        
        % Case 2: u between c and limit
    elseif ( (cx-u) * (u-ulim) > 0)
        fu = feval(fname, u);
        if(fu < fc) % no minimum
            cx_tmp =cx;
            bx = cx;    cx = u;    u = u + GOLD*(u-cx_tmp); % update positions
            fb = fc;  fc = fu;  fu = feval(fname, u);
        end
        
        % Case 3: % u beyond the limit
    elseif ( (u-ulim)*(ulim-cx) >= 0)
        u = ulim;
        fu = feval(fname, u);
        
        % Case 4: u is somewhere it should not be. Reject the parabolic fit, use default magnification
    else
        u = cx + GOLD*(cx-bx);
        fu = feval(fname, u);
    end
    
    % Eliminate oldest point and continue
    ax = bx;    bx = cx;    cx = u;
    fa = fb;  fb = fc;  fc = fu;
end

val = sort([ax, bx, cx]);
ax = val(1);
bx = val(2);
cx = val(3);
end



function [a,b] = swap(ax,bx)
a = bx;
b = ax;
end



function [xmin, fmin] = brentsMinimumSearch( ax,bx,cx, fname, TOL)
% Usage: [xmin, fmin] = brentsMinimumSearch( ax,bx,cx, fname )
% Brents method uses parabolic interpolation to find the minimum of the
% function 'fname'. For "non-cooperative" functions it falls back to golden
% section search . The minimum is searched for in the interval (ax,cx)
% with f(bx) < f(ax), f(bx) < f(cx).

% From numerical recipes 3rd edition
% Author: Simon Christoph Stein
% Date: 2015

if nargin < 5
    % Decreasing the tolarance increases speed
    % Default can be about 1e-8 for double, 1e-4 for single precision
    % (square root of machine precision)
    TOL = 1e-8;
end

ITMAX = 100; % Maximum number of iterations
R = 0.61803399; CGOLD = 1.0-R; % Golden section ratios
ZEPS = eps*1e-3; % Protect against trying to achieve fractional accuracy for a minimum that happens to be exactly zero

a=[]; b=[]; d=0.0; eptemp=[]; fu=[]; fv=[]; fw=[]; fx=[];
p=[]; q=[]; r=[]; tol1=[]; tol2=[]; u=[]; v=[]; w=[]; x=[]; xm=[];
e = 0.0; % The distance moved from the last step

% Put a and b in ascending order
if (ax < cx), a=ax; else a=cx; end;
if (ax > cx), b=ax; else b=cx; end;

% initializations
x = bx;
w = bx;
v = bx;

fx = feval(fname, x);
fv = fx;
fw = fv;

for iter=1:ITMAX % Main loop
    xm = 0.5*(a+b);
    tol1 = TOL*abs(x)+ZEPS;
    tol2=2.0*tol1;
    
    % Test if we're done
    if (abs(x-xm) <= (tol2-0.5*(b-a)) )
        fmin = fx;
        xmin = x;
        return;
    end
    
    if(abs(e) > tol1) % Construct trial parabolic fit through points x,v,w
        r = (x-w)*(fx-fv);
        q = (x-v)*(fx-fw);
        p = (x-v)*q - (x-w)*r;
        q = 2.0*(q-r);
        
        
        if (q > 0.0), p = -p; end;
        q = abs(q);
        etemp = e;
        e = d;
        
        % Check if parabolic fit is unacceptable.
        % unacceptable -> golden section step
        % acceptable -> parabolic step
        if(  ( abs(p) >= abs(0.5*q*etemp) )  ||  ( p <= q*(a-x) )  ||  (p >= q*(b-x))  )
            if (x >= xm),  e= a-x;  else  e = b-x;  end;
            d = CGOLD*e;
        else
            d = p/q;
            u = x+d;
            if(  (u-a < tol2) || (b-u < tol2) )
                if (xm-x > 0), d = tol1;  else d = -tol1; end;
            end
        end
        
    else
        if (x >= xm), e = a-x; else e = b-x; end;
        d = CGOLD*e;
    end
    
    if( abs(d) >= tol1)
        u = x+d;
    else
        if (d>0),  u = x + tol1;  else  u = x-tol1; end;
    end
    fu = feval(fname, u); % One function value per iteration
    
    % Decoide what to do with our function evaluation
    if( fu <= fx)
        if( u >= x),  a=x;  else b=x;  end;
        v = w;  w = x;  x = u;
        fv = fw;  fw = fx;  fx = fu;
    else
        if ( u<x ), a=u;  else b = u;  end;
        if( (fu <= fw) ||  ( w == x) )
            v = w;
            w = u;
            fv = fw;
            fw = fu;
        elseif( (fu <= fx) || (v == x) || (v==w)  )
            v = u;
            fv = fu;
        end
        
    end
end

error('Too many iterations in brentsMinimumSearch');
end




function [xmin, fmin] = goldenSectionSearch( ax,bx,cx, fname, TOL)
% Usage: [xmin, fmin] = goldenSectionSearch( ax,bx,cx, fname )
% Golden section search to find a minimum of the function 'fname' in the
% interval (ax,cx) with f(bx) < f(ax), f(bx) < f(cx). The converges of this
% method is linear.

% From numerical recipes 3rd edition
% Author: Simon Christoph Stein
% Date: 2015

if nargin < 5
    % Decreasing the tolarance increases speed
    % Default can be about 1e-8 for double, 1e-4 for single precision
    % (square root of machine precision)
    TOL = 1e-8;
end
R = 0.61803399;  C = 1.0-R; % Golden section ratios

x0 = ax;
x1 = [];
x2 = [];
x3 = cx;

% x0 to x1 should be the smaller segment
% the other one will be bisected
if( abs(cx-bx) > abs(bx-ax) )
    x1 = bx;
    x2 = bx+C*(cx-bx); % the new point to be tried
else
    x1 = bx-C*(bx-ax); % the new point to be tried
    x2 = bx;
end

f1 = feval(fname, x1);
f2 = feval(fname, x2);

while( abs(x3-x0) > TOL*( abs(x1) + abs(x2) )) %check fractional width
    if(f2 < f1) % see which trial point is better, update boundary points
        x0 = x1;  x1 = x2;  x2 = R*x2 + C*x3;
        f1 = f2;  f2 = feval(fname, x2);
    else
        x3 = x2;  x2 = x1; x1 = R*x1+C*x0;
        f2 = f1;  f1 = feval(fname, x1);
    end
end

if (f1<f2)
    xmin = x1;
    fmin = f1;
else
    xmin = x2;
    fmin = f2;
end


end





