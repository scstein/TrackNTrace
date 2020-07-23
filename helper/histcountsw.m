function [n,edges,bin] = histcountsw(x, w, varargin)
%HISTCOUNTSW Weighted Histogram Bin Counts
% [n,edges,bin] = histcountsw(x, w, varargin)
% Like histcounts but with an additional weighting of the data x
% with w. For w = 1 they give identical output.
% Elements with the weight NaN are excluded.
% See HISTCOUNTS for more details.
%
% There are the following additional normaizaltion options:
%  'weights' (Default): Normalized to the weight of the counts.
%  'cumweight'        : Cummulative weights.
%  'weightdensity'    : Same as weights, but divded by the bin width.
% For the options 'count', 'cumcount' and 'countdensity' the weights are
% normalized to a mean of 1.
%
% Christoph Thiele, 2019

    if nargin<2 || isempty(w)
        w = 1;
    elseif ischar(w)
        varargin = [{w}, varargin];
        w = 1;
    end
    
    ind = isnan(w);
    if any(ind)
        w = w(~ind);
        x = x(~ind);
    end
%     norm = {'count', 'probability', 'countdensity', 'pdf', 'cumcount', 'cdf'};
    normidx = find(strcmpi(varargin,'Normalization'),1,'last');
    if ~isempty(normidx) 
        if numel(normidx)+1>numel(varargin)
            error(message('MATLAB:histcounts:ArgNameValueMismatch'));
        end
        norm = varargin{normidx+1};
        
        if contains(norm, 'weight')
            varargin{normidx+1} = strrep(norm,'weight','count');
        elseif contains(norm, 'count')
            w = w/sum(w(:))*numel(w);
        end
    else
        norm ='weights';
    end
    
    [n,edges,bin] = histcounts(x, varargin{:});
    ind = bin~=0; % zero indicates values out of the bounds or NaNs.
    
    % expand w if singleton
    if isscalar(w)
        w = w * ones(size(x));
    elseif numel(bin) ~= numel(w)
        error('Dimensions of values and weights do not match.');
    end
    if any(ind)
        bin = bin(:);
        w = w(:);
        n = accumarray(bin(ind),w(ind),size(n'))';% sum up w for each bin.
    end
    switch norm
        case {'countdensity', 'weightdensity'}
            n = n./double(diff(edges));
        case {'cumcount', 'cumweight'}
            n = cumsum(n);
        case 'probability'
            n = n / sum(w);
        case 'pdf'
            n = n/sum(w)./double(diff(edges));
        case 'cdf'
            n = cumsum(n / sum(w));
    end
end

