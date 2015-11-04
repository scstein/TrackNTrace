function h = createOptionsGUI( name, par_text, par_type, par_value, par_tooltip)
%CREATEOPTIONSGUI Summary of this function goes here
%   Detailed explanation goes here

h = figure('NumberTitle','off');

num_pars = numel(par_text);
% set(h,'Name',[name, ' Options']);
set(h,'Units','points');
set(h,'Resize','off');
set(h,'ToolBar','none');
set(h,'MenuBar','none');

bg_color = get(h,'Color');

fontSize = 10;
row_spacing = 2*fontSize;
col_spacing = 1*fontSize;
win_heigth = (num_pars+1)*row_spacing+fontSize*2;

elemHeigth = 1.4*fontSize;
editSize = 6*fontSize; % Size of edit fields


% Title text
text = [name, ' Options'];
text_pos = [col_spacing, win_heigth-2*fontSize, 1,elemHeigth];
h_text = uicontrol('Parent',h, 'Units','points', 'Style','text',...
    'BackgroundColor',bg_color,'FontSize',fontSize,'FontWeight','bold','String',text);
% Resize text to full width
extent = get(h_text,'Extent');
text_pos(3) = extent(3);
set(h_text,'Position',text_pos);

max_width = text_pos(3)+2*col_spacing;

for iP = 1:num_pars
    text = par_text{iP};
    ptype = par_type{iP};
    value = par_value{iP};
    
    % Add text
    text_pos = [col_spacing, win_heigth-5*fontSize-(iP-1)*row_spacing, 1, elemHeigth];
    h_text = uicontrol('Parent',h, 'Units','points','Style','text', ...
        'BackgroundColor',bg_color,'FontSize',fontSize,'String',text);
    extent = get(h_text,'Extent');
    text_pos(3) = extent(3);
    set(h_text,'Position',text_pos);
    textwidth = text_pos(3);
    
    % Set Tooltip
    set(h_text, 'TooltipString',par_tooltip{iP});
    
    % Add uicontrol for value
    value_position = [col_spacing+textwidth+col_spacing, win_heigth-5*fontSize-(iP-1)*row_spacing, 1, 1];
    switch ptype
        case 'float'
            val_pos = [value_position(1), value_position(2), editSize, elemHeigth];
            h_val=uicontrol('Parent',h, 'Units','points', 'Position', val_pos, ...
                'Style','edit','FontSize',fontSize,'String',num2str(value));% ,'Callback', @changeText);
        case 'int'
            val_pos = [value_position(1), value_position(2), editSize, elemHeigth];
            h_val= uicontrol('Parent',h, 'Units','points', 'Position', val_pos, ...
                'Style','edit','FontSize',fontSize,'String',num2str(value));% ,'Callback', @changeText);
        case 'bool'
            val_pos = [value_position(1), value_position(2), fontSize, elemHeigth];
            h_val= uicontrol('Parent',h, 'Units','points', 'Position', val_pos, ...
                'Style','checkbox','FontSize',fontSize,'Value',value);% ,'Callback', @changeText);
        case 'list'
            % TODO: Determine length of largest string, set extent based on this
            val_pos = [value_position(1), value_position(2), editSize, elemHeigth];
            h_val = uicontrol('Parent',h, 'Units','points', 'Position', val_pos, ...
                'Style','popupmenu','FontSize',fontSize,'String',value);% ,'Callback', @changeText);
            extent = get(h_val,'Extent');
            val_pos(3) = extent(3)+16;
            set(h_val,'Position',val_pos);
        otherwise
            error('Unkown parameter type');
    end
    
    % Set Tooltip
    set(h_val, 'TooltipString',par_tooltip{iP});
    
    overall_width = 3*col_spacing + text_pos(3) + val_pos(3);
    if(overall_width> max_width)
        max_width = overall_width;
    end
end

set(h,'Position', [0,0, max_width, win_heigth]);
movegui(h,'center');

% h_text = uicontrol('Parent',h_window, 'Units','normalized', 'Position',[0. 0.95 1 0.05], ...
%     'Style','text','FontSize',16,'String','Hello world.');
% h_button = uicontrol('Parent',h_window, 'Units','normalized', 'Position',[0.4 0.85 0.2 0.05], ...
%     'Style','pushbutton','String','Colorbutton','Callback',@changeColor);
% h_popup = uicontrol('Parent',h_window, 'Units','normalized', 'Position',[0.4 0.75 0.2 0.05], ...
%     'Style','popupmenu','String',{'Choice 1','Choice 2', 'Choice 3'},'Callback',@printPopupChoice);
% h_slider = uicontrol('Parent',h_window, 'Units','normalized', 'Position',[0.3 0.65 0.4 0.05], ...
%     'Style','slider','Min',0,'Max',100,'SliderStep',[0.05 0.2],'Callback', @printValue);
% h_edit = uicontrol('Parent',h_window, 'Units','normalized', 'Position',[0.3 0.55 0.4 0.05], ...
%     'Style','edit','String',get(h_text,'String') ,'Callback', @changeText);
% h_axes = axes('Parent',h_window, 'Units','normalized', 'Position',[0.225 0.05 0.5 0.45], ...
%     'ButtonDownFcn',@printPosition,'NextPlot','replacechildren');


end

