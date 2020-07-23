function h = histogramw(varargin)
%HISTOGRAMW Plots a weighted histogram
% h = HISTOGRAM(X,W,...)
% Like histogram but with an additional weighting of the data x
% with w. For w = 1 they give identical output.
%
% h = HISTOGRAM(AX,__)
% 
% Supports all the options from histogram, except of the input types 
% datetime, duration and categorial. See HISTOGRAM for more details.
%
% The binning is done by histcountsw. See HISTCOUNTSW for details on the
% normalization of the histogram.
%
%
% Note: If the Data property of histogram object is changed the updated
% histogram will not be weighted anymore.
%
% Christoph Thiele, 2019

    % Select the input arguments for histcountsw. 
    [cax,args] = axescheck(varargin{:});
    countarglist = {'NumBins','BinEdges','BinLimits',...
        'BinWidth','Normalization','BinMethod'};
    countargsidx = cellfun(@(arg)find(strcmpi(arg,args)),countarglist,'UniformOutput',false);
    if numel([countargsidx{:}])>0
        % There are extra arguments for histcounts
        countargsidx = [countargsidx{:}]+[0; 1];
        countargsidx = [1:min([numel(args),find(cellfun(@ischar,args),1)-1]), countargsidx(:)'];
        if countargsidx(end)>numel(args)
            error(message('MATLAB:histogram:ArgNameValueMismatch'));
        end
    else
        % There are no extra arguments for histcounts
        countargsidx = 1:min([numel(args),find(cellfun(@ischar,args),1)-1]);
    end
    countarg = args(countargsidx);

    % All other arguments are passed to histogram.
    histargs = args(setxor(1:numel(args),countargsidx));
    if ~isempty(cax)
        histargs = [histargs, {'Parent', cax}];
    end

    [n,edges] = histcountsw(countarg{:});

    histObj = histogram('BinEdges',edges,'BinCounts',n,histargs{:});

    if nargout > 0
        h = histObj;
    end

end

