function createOptionsPanel( h_panel, plugin_name, par_name, par_type, par_defaultValue, par_tooltip, inputOptions)


% All paramters are saved within the options struct. Variable names are
% derived from the given parameter name by erasing all white spaces.
if nargin < 7 || isempty(inputOptions)
    options = [];
    options.plugin_name = plugin_name;
else
    if strcmp(plugin_name,inputOptions.plugin_name)
        options = inputOptions;
    else
        error('createOptionsPanel: Tried to create options for plugin ''%s'' with inputOptions for ''%s''.',plugin_name,inputOptions.plugin_name);
    end
end

% Check if properties are defined for all parameters
nElem = [numel(par_name), numel(par_type), numel(par_defaultValue), numel(par_tooltip)];
if(numel(unique(nElem)) >1)
   error('createOptionsPanel: Error: par_name, par_type, par_defaultValue, par_tooltip must be specified for every parameter to create a plugin!.');
end

% Number of parameters
num_pars = numel(par_name);

% Remove all existing uicontrols from the panel (in case it was used before)
child_handles = get(h_panel,'Children');
for iC = 1:numel(child_handles)
    delete(child_handles(iC));
end

% Compute the variable name for every parameter
% Replace white space ' ' by underscore '_' & remove dots '.'
struct_varNames = cell(num_pars);
for iP = 1:numel(par_name);
    stripped_name =  par_name{iP};
    stripped_name(stripped_name == ' ') = '_';
    stripped_name(stripped_name == '.') = '';
    struct_varNames{iP} = stripped_name;
end


set(h_panel,'Title',[plugin_name, ' Options']);
set(h_panel,'Units','Points');


% These values are fixed!
panel_pos = get(h_panel,'Position');
panel_width = panel_pos(3);
bg_color = get(h_panel,'BackgroundColor');

% These variables influence the appearence of the panel
fontSize = 10; % Note: Almost everything scales with the font size
row_gap = 1*fontSize; % Gap between two rows
col_gap = 1*fontSize; % Gap between controls for two parameters
text_gap = 0.5*fontSize; % Gap between text and value field
elemHeight = 1.4*fontSize; % Height of uicontrols
editWidth = 6*fontSize; % width of edit fields
maxElemPerRow = 2; % Maximum number of elements in a row

% 
elemRowIdx = 1;
left_pos = text_gap;
bott_pos = -0.5*fontSize;
newRow();

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
    text_pos = [1, 1, 1, elemHeight];
    h_text = uicontrol('Parent',h_panel, 'Units','points','Style','text', ...
        'BackgroundColor',bg_color,'FontSize',fontSize,'String',pName);
    extent = get(h_text,'Extent');
    text_pos(3) = extent(3);
    set(h_text,'Position',text_pos);
    textwidth = text_pos(3);
    
    % Set Tooltip
    set(h_text, 'TooltipString',par_tooltip{iP});
    
    % Add uicontrol for value
    value_position = [1, 1, 1, 1];
    switch pType
        case 'float'
            val_pos = [value_position(1), value_position(2), editWidth, elemHeight];
            h_val=uicontrol('Parent',h_panel, 'Units','points', 'Position', val_pos, ...
                'Style','edit','BackgroundColor', 'w','FontSize',fontSize,'String',num2str(pValue));
        case 'int'
            val_pos = [value_position(1), value_position(2), editWidth, elemHeight];
            h_val= uicontrol('Parent',h_panel, 'Units','points', 'Position', val_pos, ...
                'Style','edit','BackgroundColor', 'w','FontSize',fontSize,'String',num2str(pValue));
        case 'bool'
            val_pos = [value_position(1), value_position(2), 11.5, elemHeight];
            h_val= uicontrol('Parent',h_panel, 'Units','points', 'Position', val_pos, ...
                'Style','checkbox','FontSize',fontSize,'Value',pValue,'String','');
        case 'list'
            val_pos = [value_position(1), value_position(2), editWidth, elemHeight];
            h_val = uicontrol('Parent',h_panel, 'Units','points', 'Position', val_pos, ...
                'Style','popupmenu','BackgroundColor', 'w','FontSize',fontSize,'String',par_defaultValue{iP});
            % Determine neccessary width of popup based on choosable options
            maxCharacters = 0;
            for iOption = 1:numel(par_defaultValue{iP})
                characters = numel(par_defaultValue{iP}{iOption});
                if(characters>maxCharacters); maxCharacters = characters; end;
            end
            set(h_val,'Value', 1); % Select first item by default
            val_pos(3) = (maxCharacters+1)*fontSize/2+16; % Note: fontSize/2 is an (arbitrary) approximation of the font width
            set(h_val,'Position',val_pos);
            
            % Special case for lists: The given default value is the list of
            % options. If no selected option was given, the first item will
            % be automatically selected. If a selected option was given,
            % we have to set the popup accordingly.
            if(isfield(options,pStructVarName))
                setPopup(h_val, pValue);
            end
            
        otherwise
            error('createOptionsPanel: Unknown parameter type');
    end
    valwidth = val_pos(3);
    
    % Set Tooltip
    set(h_val, 'TooltipString',par_tooltip{iP});
    
    % Set callback
    set(h_val, 'Callback', {@callback_saveOnChange,pStructVarName,pType});
    % Execute callback once to save the parameter starting value
    callback_saveOnChange(h_val, [], pStructVarName, pType);
    
    % Width of 'text value'
    overall_width = textwidth + text_gap + valwidth;
    
    % Check if element should be placed in this row
    elemFitsInRow = (left_pos + overall_width < panel_width);
    if(~elemFitsInRow || elemRowIdx > maxElemPerRow)
        % Advance to new row
        newRow();
        placeElement();
    else
        placeElement();
        % Advance to next column
        elemRowIdx = elemRowIdx+1;
         left_pos = left_pos + overall_width + col_gap;
    end
end

panel_newheight = abs(bott_pos)+0.5*fontSize;

% Offset bottom positions of all uicontrols by the panel height
child_handles = get(h_panel,'Children');
for iC = 1:numel(child_handles)
   pos = get(child_handles(iC), 'Position');
   pos(2) = pos(2) + panel_newheight;
   set(child_handles(iC), 'Position',pos);
end

% Set Panel position to keep the top left corner at the same position as
% before this function call 
panel_pos(2) = panel_pos(2) + panel_pos(4) - panel_newheight;
panel_pos(4) = panel_newheight;
set(h_panel,'Position',panel_pos);



    function newRow()
        left_pos = text_gap;
        bott_pos = bott_pos-row_gap-elemHeight;
        elemRowIdx = 1;
    end

    function placeElement()
                % Place element
        text_pos = get(h_text,'Position');
        text_pos(1:2) = [left_pos,bott_pos];
        set(h_text,'Position',text_pos);

        val_pos = get(h_val, 'Position');
        val_pos(1:2) = [left_pos+textwidth+text_gap, bott_pos];
        set(h_val,'Position',val_pos);
    end

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
                error('createOptionsPanel:callback_saveOnChange: Unknown parameter type');
        end
        setappdata(h_panel,'options',options); % Update stored options
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

