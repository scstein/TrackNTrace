function [ax,rgb] = im2rgb(varargin)
%IM2RGB Converts a colormaped imagesc (imagesc) into RGB.
% im2rgb uses the alphadata to blend the image with the background color of the
% axis. This is necessary when exporting, as the background color of the axis is
% lost on export.
%
% [ax,rgb] = im2rgb(ax)
% 
% ax  - axes or image to convert.
%
% ax  - converted axes
% rgb - generated rgb map
%
% im2rgb is restricted to single axes and single images inside that axes.
%
% (C) Christoph Thiele, 2020.


% get axes
ax = gobjects(0);
if nargin>0
    if isgraphics(varargin{1})
        hdl = ancestor(varargin{1},'axes');
        if isgraphics(hdl)
            ax = hdl;
        else
            hdl = ancestor(varargin{1},'figure');
            ax = hdl.CurrentAxes;
        end
    end
end
if isempty(ax) || ~isgraphics(ax,'Axes')
    ax = gca;
end
% get image
if nargin>0 && isgraphics(varargin{1},'image')
    im = varargin{1};
else
    im = findobj(ax,'Type','image');
end

if isempty(im) || ~isprop(im,'CData') || isempty(im.CData)
    rgb = [];
    return
end
% get rgb values
cdata = im.CData;
if size(cdata,3)==1
    cmap = ax.Colormap;
    cmap_len = size(cmap,1);
    
    % Adapted from  John Iversen (2020). freezeColors / unfreezeColors (https://www.mathworks.com/matlabcentral/fileexchange/7943-freezecolors-unfreezecolors), MATLAB Central File Exchange. Retrieved March 20, 2020.
    %convert cdata to indexes into colormap
    if strcmpi(im.CDataMapping,'scaled')
        cdatalim = ax.CLim;
        cdata = ceil( (double(cdata) - cdatalim(1)) / (cdatalim(2)-cdatalim(1)) * cmap_len);
    else %direct mapping
        cdata = floor(cdata);
    end

    rgb = ind2rgb(cdata,cmap); % ind2rgb automatically limits cdata to [1 cmap_len]
    
else
    rgb = cdata;
end

% get alpha values
adata = im.AlphaData;
if ~strcmpi(im.AlphaDataMapping,'none')
    % we need to apply the alphamap of the axes
    amap = ax.Alphamap;
    amap_len = numel(amap);
    adatalim = ax.ALim;
    if strcmpi(im.AlphaDataMapping,'scaled')
        adata = ceil( (double(adata) - adatalim(1)) / (adatalim(2)-adatalim(1)) * amap_len);
    else
        if isinteger(adata)
            adata = double(adata) + 1;
        else
            adata = floor(adata);
        end
    end
    adata = amap(max(1,min(amap_len,adata)));
end
adata = max(0,min(1,adata));

% combine using axes color as background
bgcolor = ax.Color;
if strcmpi(bgcolor,'none')
    warning('Axes background is transparent. Using white instead.');
    bgcolor = [1 1 1];
end
bgcolor = shiftdim(bgcolor(:),-2);

rgb = rgb .* adata + bgcolor .* (1-adata);

%%
im.CData = rgb;
im.AlphaData = 1;
im.AlphaDataMapping = 'none';

end



