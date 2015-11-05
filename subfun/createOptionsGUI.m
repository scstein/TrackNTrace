function h = createOptionsGUI( plugin_name, par_name, par_type, par_defaultValue, par_tooltip, inputOptions)
%CREATEOPTIONSGUI Summary of this function goes here
%   Detailed explanation goes here

% All paramters are saved within the options struct. Variable names are
% derived from the given parameter name by erasing all white spaces.
if nargin < 6 || isempty(inputOptions)
    options = [];
    options.plugin_name = plugin_name;
else
    if(plugin_name == inputOptions.plugin_name)
        options = inputOptions;
    else
        error('createOptionsGUI: Tried to create options for plugin ''%s'' with inputOptions for ''%s''.');
    end
end

% strip white spaces from paramter names
struct_varNames = cell(numel(par_name));
for iP = 1:numel(par_name);
   stripped_name =  par_name{iP};
   stripped_name(stripped_name == ' ') = [];
   struct_varNames{iP} = stripped_name;
end


h = figure('NumberTitle','off');

num_pars = numel(par_name);
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
plugin_name = [plugin_name, ' Options'];
text_pos = [col_spacing, win_heigth-2*fontSize, 1,elemHeigth];
h_text = uicontrol('Parent',h, 'Units','points', 'Style','text',...
    'BackgroundColor',bg_color,'FontSize',fontSize,'FontWeight','bold','String',plugin_name);
% Resize text to full width
extent = get(h_text,'Extent');
text_pos(3) = extent(3);
set(h_text,'Position',text_pos);

max_width = text_pos(3)+2*col_spacing;

% Create controls for all parameters
for iP = 1:num_pars
    pName = par_name{iP};
    pStructVarName = struct_varNames{iP};
    pType = par_type{iP};
    
    % Check if value for parameter was given
    if(isfield(options,pStructVarName))
        pValue = options.(pStructVarName);
    else
        pValue = par_defaultValue{iP};
    end
    
    % Add text
    text_pos = [col_spacing, win_heigth-5*fontSize-(iP-1)*row_spacing, 1, elemHeigth];
    h_text = uicontrol('Parent',h, 'Units','points','Style','text', ...
        'BackgroundColor',bg_color,'FontSize',fontSize,'String',pName);
    extent = get(h_text,'Extent');
    text_pos(3) = extent(3);
    set(h_text,'Position',text_pos);
    textwidth = text_pos(3);
    
    % Set Tooltip
    set(h_text, 'TooltipString',par_tooltip{iP});
    
    % Add uicontrol for value
    value_position = [col_spacing+textwidth+col_spacing, win_heigth-5*fontSize-(iP-1)*row_spacing, 1, 1];
    switch pType
        case 'float'
            val_pos = [value_position(1), value_position(2), editSize, elemHeigth];
            h_val=uicontrol('Parent',h, 'Units','points', 'Position', val_pos, ...
                'Style','edit','FontSize',fontSize,'String',num2str(pValue));
        case 'int'
            val_pos = [value_position(1), value_position(2), editSize, elemHeigth];
            h_val= uicontrol('Parent',h, 'Units','points', 'Position', val_pos, ...
                'Style','edit','FontSize',fontSize,'String',num2str(pValue));
        case 'bool'
            val_pos = [value_position(1), value_position(2), 11.5, elemHeigth];
            h_val= uicontrol('Parent',h, 'Units','points', 'Position', val_pos, ...
                'Style','checkbox','FontSize',fontSize,'Value',pValue,'String','');
        case 'list'
            val_pos = [value_position(1), value_position(2), editSize, elemHeigth];
            h_val = uicontrol('Parent',h, 'Units','points', 'Position', val_pos, ...
                'Style','popupmenu','FontSize',fontSize,'String',par_defaultValue{iP});
            % Determine neccessary width of popup based on choosable options
            maxCharacters = 0;
            for iOption = 1:numel(par_defaultValue{iP})
                characters = numel(par_defaultValue{iP}{iOption});
                if(characters>maxCharacters); maxCharacters = characters; end;                    
            end
            set(h_val,'Value', 1);
            val_pos(3) = maxCharacters*fontSize/2+16; % Note: fontSize/2 is an (arbitrary) approximation of the font width
            set(h_val,'Position',val_pos);
            
            % Special case for lists: The given default value is the list 
            % options. If no selected option was given, the first item will
            % be automatically selected. If a selected option was given, 
            % we have to set the popup accordingly.
            if(isfield(options,pStructVarName))
                setPopup(h_val, pValue);
            end
            
        otherwise
            error('createOptionsGUI: Unknown parameter type');
    end
    
    % Set Tooltip
    set(h_val, 'TooltipString',par_tooltip{iP});
    
    % Set callback
    set(h_val, 'Callback', {@callback_saveOnChange,pStructVarName,pType});
    % Execute callback once to save the parameter starting value
    callback_saveOnChange(h_val, [], pStructVarName, pType);
    
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


    function callback_saveOnChange(uiObj, eventdata, structVarName, pType)       
        switch pType
            case 'float'
                pValue = str2double(get(uiObj,'String'));
                set(uiObj, 'String', num2str(pValue)); % Synchronise text field
                options.(structVarName) = pValue;
            case 'int'
                pValue = round(str2double(get(uiObj,'String'))); % get integer value
                set(uiObj, 'String', num2str(pValue)); % Synchronise text field
                options.(structVarName) = pValue;
            case 'bool'
                pValue = logical(get(uiObj,'Value'));
                options.(structVarName) = pValue;
            case 'list'
                choices = get(uiObj,'String');
                optionString = choices{get(uiObj,'Value')};
                options.(structVarName) = optionString;
            otherwise
                error('createOptionsGUI:callback_saveOnChange: Unknown parameter type');                              
        end
        setappdata(h,'options',options); % Update stored options
    end


% Sets a popup menu to the option with name 'optionString'
    function setPopup(hObj,optionString)
        choices = get(hObj,'String');
        for idx = 1:length(choices)
            if strcmp(choices{idx}, optionString)
                set(hObj,'Value',idx);
                return
            end
        end
        error('invalid option ''%s'' to setPopup tag: %s',optionString,get(hObj,'Tag'));
    end
end

