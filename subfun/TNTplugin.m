classdef TNTplugin < handle % Inherit from handle class
    %TNTPLUGIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        param_specification % This saves the names and information about all parameters
        name
        type
        info
        initFunc
        mainFunc
        postFunc
    end
    
    properties(Access = private)
        options
    end
    
    methods
        % Constructor
        function obj = TNTplugin(name, type, mainFunc)
            % Check input
            if isempty(name)
                error('Empty plugin name not allowed!');
            end
            
            switch type
                case 1 % Candidate method
                case 2 % Fitting method
                case 3 % Tracking method
                otherwise
                    warning('Invalid plugin type ''%i'' when constructing plugin ''%s''', type, name);
            end
            
            if isempty(mainFunc) || ~isa(mainFunc,'function_handle')
                error('Invalid or empty mainFunc handle when constructing plugin ''%s''', name);
            end
            
            % Internal
            obj.param_specification = cell(0,4);
            
            % Mandatory
            obj.name = name;
            obj.type = type;
            obj.mainFunc = mainFunc;
            
            % Optional
            obj.info = '';
            obj.initFunc = [];
            obj.postFunc = [];
        end
        
        % Add parameter to this plugin
        function add_param(obj, par_name, par_type, par_settings, par_tooltip )
            %   -------------- HOW TO ADD PARAMETERS --------------
            % Add parameters from inside a plugins M-file using:
            %   add_param(par_name, par_type, par_settings, par_tooltip)
            %  Example:     add_param('particle_radius', 'int', {3, 0, inf}, 'Approximate spot radius');
            %
            % par_name: Name of parameter
            %   These translate to the names of variables (>NAME<) inside the options struct this plugin
            %   saves to by replacing all white spaces ' ' with underscores '_' and removing dots '.'.
            %   Use options.>NAME< within your plugins function to address these parameters
            %
            % par_type: Type of parameter
            %   One of 'float', 'int', 'bool','string','list'
            %
            % par_settings: Settings needed to setup parameter, depends on type
            %   'float':  Cell array {defaultValue, lowerBound, upperBound}
            %   'int':    Cell array {defaultValue, lowerBound, upperBound}
            %   'bool':   The default value true/false
            %   'string': The default value (a string)
            %   'list':   A cell array string list of possible choices for 'list' (first entry is default)
            %   'filechooser': A cell array string list {'default directory','filterEnding'}
            % 
            % par_tooltip: Tooltip shown when hovering over the parameter with the mouse
            
            
            %   -------------- Input checking --------------
            if( isempty(par_name))
                error('add_plugin_param: Missing parameter name!')
            end
            if( isempty(par_type))
                error('add_plugin_param: Missing parameter type!')
            end
            
            switch par_type
                case 'float'
                    if(~iscell(par_settings) || length(par_settings) ~= 3)
                        error('add_plugin_param: Failed adding param ''%s'' of type ''%s''. Settings need to be {defaultValue, lowerBound, upperBound}!', par_name, par_type)
                    end
                case 'int'
                    if(~iscell(par_settings) || length(par_settings) ~= 3)
                        error('add_plugin_param: Failed adding param ''%s'' of type ''%s''. Settings need to be {defaultValue, lowerBound, upperBound}!', par_name, par_type)
                    end
                case 'bool'
                    if(~islogical(par_settings) || length(par_settings)~= 1)
                        error('add_plugin_param: Failed adding param ''%s'' of type ''%s''. Settings need to be [defaultValue] (true/false)!', par_name, par_type)
                    end
                case 'string'
                    if(~ischar(par_settings))
                        error('add_plugin_param: Failed adding param ''%s'' of type ''%s''. Settings need to be [''defaultValue'']!', par_name, par_type)
                    end
                case 'list'
                    if(~iscell(par_settings))
                        error('add_plugin_param: Failed adding param ''%s'' of type ''%s''. Settings need to be {''choice1'',''choice2'',...})!', par_name, par_type)
                    end
                case 'filechooser'
                    if(~iscell(par_settings) || length(par_settings) ~= 2)
                        error('add_plugin_param: Failed adding param ''%s'' of type ''%s''. Settings need to be {''default directory'',''filterEnding''}!', par_name, par_type)
                    end
                otherwise
                    error('add_plugin_param: Unknown parameter type ''%s'' of parameter ''%s''.', par_type, par_name);
            end
            
            % --- Add parameter information ---
            obj.param_specification = vertcat(obj.param_specification, {par_name, par_type, par_settings, par_tooltip});
        end
        
        % Add data to the options that have nothing to do with the GUI and its parameters.
        function addInternalsToOptions(obj)
            % save plugin name
            obj.options.plugin_name = obj.name;
            
            % Save the functions in options
            obj.options.initFunc = obj.initFunc;
            obj.options.mainFunc = obj.mainFunc;
            obj.options.postFunc = obj.postFunc;
        end
        
        % Retrieves the options for this plugin (function handles, all parameter values etc.)
        function options = getOptions(obj)
            % Add non-GUI parameters to options
            obj.addInternalsToOptions();
            
            % Return the options
            options = obj.options;
        end
        
        % Sets the options for this plugin
        function setOptions(obj, inputOptions)
            % All paramters are the options struct. Variable names are
            % derived from the given parameter name by erasing all white spaces.
            if isempty(inputOptions)
                obj.options = [];
            else
                if ~isfield(inputOptions,'plugin_name')
                    error('createOptionsPanel: Tried to create options for plugin ''%s'' with inputOptions carrying no field plugin_name.',plugin_name);
                end
                
                if strcmp(obj.name,inputOptions.plugin_name)
                    obj.options = inputOptions;
                else
                    error('createOptionsPanel: Tried to create options for plugin ''%s'' with inputOptions for ''%s''.',plugin_name,inputOptions.plugin_name);
                end
            end
            
            obj.addInternalsToOptions();
        end
        
        % Adjust the input panel to a graphical representation of this plugin in its current 
        % state (values stored inside obj.options), where its parameters can be manipulated.
        function createOptionsPanel(obj, h_panel)            
            % unpack parameter specification
            par_name = obj.param_specification(:,1);
            par_type = obj.param_specification(:,2);
            par_settings = obj.param_specification(:,3);
            par_tooltip = obj.param_specification(:,4);
            
            % Check if properties are defined for all parameters
            nElem = [numel(par_name), numel(par_type), numel(par_settings), numel(par_tooltip)];
            if(numel(unique(nElem)) >1)
                error('createOptionsPanel: Error: par_name, par_type, par_settings, par_tooltip must be specified for every parameter to create plugin ''%s''!.',plugin_name);
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
            
            set(h_panel,'Title',[obj.name, ' Options']);
            set(h_panel,'Units','Points');
            
            
            % These values are fixed!
            panel_pos = get(h_panel,'Position');
            panel_width = panel_pos(3);
            bg_color = get(h_panel,'BackgroundColor');
            
            % These variables influence the appearence of the panel
            fontSize = 8; % Note: Almost everything scales with the font size
            row_gap = 1*fontSize; % Gap between two rows
            col_gap = 1*fontSize; % Gap between controls for two parameters
            text_gap = 0.5*fontSize; % Gap between text and value field
            elemHeight = 1.4*fontSize; % Height of uicontrols
            editWidth = 6*fontSize; % width of edit fields
            maxElemPerRow = 10; % Maximum number of elements in a row
            
            %
            elemRowIdx = 1;
            left_pos = text_gap;
            bott_pos = -8;
            newRow();
            
            % Create controls for all parameters
            for iP = 1:num_pars
                pName = par_name{iP};
                pStructVarName = struct_varNames{iP};
                pType = par_type{iP};
                
                % Check if value for parameter was given
                if(isfield(obj.options,pStructVarName))
                    pValue = obj.options.(pStructVarName);
                else
                    % Currently the default value is either the value of par_settings itself
                    % or the first value inside the settings cell array
                    if(iscell(par_settings{iP}))
                        pValue = par_settings{iP}{1};
                    else
                        pValue = par_settings{iP};
                    end
                end
                
                pSettings = par_settings{iP};
                
                
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
                    case 'string'
                        val_pos = [value_position(1), value_position(2), 2*editWidth, elemHeight];
                        h_val= uicontrol('Parent',h_panel, 'Units','points', 'Position', val_pos, ...
                            'Style','edit','BackgroundColor', 'w','FontSize',fontSize,'String',pValue);
                    case 'list'
                        val_pos = [value_position(1), value_position(2), editWidth, max(elemHeight,17)]; % popupmenus are at least 17 points high to be displayed right
                        h_val = uicontrol('Parent',h_panel, 'Units','points', 'Position', val_pos, ...
                            'Style','popupmenu','BackgroundColor', 'w','FontSize',fontSize,'String',par_settings{iP});
                        % Determine neccessary width of popup based on choosable options
                        maxCharacters = 0;
                        for iOption = 1:numel(par_settings{iP})
                            characters = numel(par_settings{iP}{iOption});
                            if(characters>maxCharacters); maxCharacters = characters; end;
                        end
                        set(h_val,'Value', 1); % Select first item by default
                        val_pos(3) = (maxCharacters+1)*fontSize/2+16; % Note: fontSize/2 is an (arbitrary) approximation of the font width
                        set(h_val,'Position',val_pos);
                        
                        % For lists the given settings value is the list of
                        % options. If no selected option was given, the first item will
                        % be automatically selected. If a selected option was given,
                        % we have to set the popup accordingly.
                        if(isfield(obj.options,pStructVarName))
                            setPopup(h_val, pValue);
                        end
                    case 'filechooser'
                        val_pos = [value_position(1), value_position(2), 3*editWidth, elemHeight];
                        h_val= uicontrol('Parent',h_panel, 'Units','points', 'Position', val_pos, ...
                            'Style','edit','BackgroundColor', 'w','FontSize',fontSize,'String',pValue);
                        
                        % Create a filechooser button
                        button_pos = [value_position(1), value_position(2), editWidth, elemHeight*1.1];
                        h_button= uicontrol('Parent',h_panel, 'Units','points', 'Position', button_pos, ...
                            'Style','pushbutton','FontSize',fontSize,'String','Select');
                        set(h_button, 'Callback', {@callback_filechooserButton, h_val, pStructVarName, pSettings});
                        buttonwidth = button_pos(3);
                        
                        % Set Tooltip
                        set(h_button, 'TooltipString',par_tooltip{iP});
                    otherwise
                        error('createOptionsPanel: Unknown parameter type %s for parameter %s!',pType, pName);
                end
                valwidth = val_pos(3);
                
                % Set Tooltip
                set(h_val, 'TooltipString',par_tooltip{iP});
                
                % Set callback
                set(h_val, 'Callback', {@callback_saveOnChange,pStructVarName,pType, pSettings});
                % Execute callback once to save the parameter starting value
                callback_saveOnChange(h_val, [], pStructVarName, pType, pSettings);
                               
                % Place element in panel
                placeElement();
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
                elemRowIdx = 1;
                left_pos = text_gap;
                bott_pos = bott_pos-row_gap-elemHeight;
            end
            
            % Place element in the panel
            function placeElement()
                % Filechoosers always occupy their own row
                if strcmp(pType,'filechooser')
                    % Width of 'text value'
                    overall_width = textwidth + text_gap + valwidth + text_gap + buttonwidth;
                    
                    % Check if element should be placed in this row
                    elemFitsInRow = (left_pos + overall_width < (panel_width -text_gap) );
                    if(~elemFitsInRow)
                        newRow();
                    end
                    
                    text_pos = get(h_text,'Position');
                    text_pos(1:2) = [left_pos,bott_pos];
                    set(h_text,'Position',text_pos);
                    
                    val_pos = get(h_val, 'Position');
                    val_pos(1:2) = [left_pos+textwidth+text_gap, bott_pos];
                    set(h_val,'Position',val_pos);
                    
                    button_pos = get(h_button, 'Position');
                    button_pos(1:2) = [left_pos+textwidth+text_gap+text_gap+valwidth+text_gap, bott_pos];
                    set(h_button,'Position',button_pos);
                    
                    left_pos = left_pos + overall_width + col_gap; % Next column
                    
                    % If next element exceeds the max elements per row, begin new row
                    % (Except for the last element)
                    if( elemRowIdx > maxElemPerRow && iP~=num_pars)
                        newRow();
                    end
                else % For all other parameter types:
                    % Width of 'text value'
                    overall_width = textwidth + text_gap + valwidth;
                    
                    % Check if element should be placed in this row
                    elemFitsInRow = (left_pos + overall_width < (panel_width -text_gap) );
                    if(~elemFitsInRow)
                        newRow();
                    end
                    
                    text_pos = get(h_text,'Position');
                    text_pos(1:2) = [left_pos,bott_pos];
                    set(h_text,'Position',text_pos);
                    
                    val_pos = get(h_val, 'Position');
                    % Note: +(elemHeight-val_pos(4))/2 is needed to center Elements that are note elemHeight high (like poupmenus which need at least 17 points in height.
                    val_pos(1:2) = [left_pos+textwidth+text_gap, bott_pos+(elemHeight-val_pos(4))/2]; 
                    set(h_val,'Position',val_pos);
                    
                    left_pos = left_pos + overall_width + col_gap; % Next column
                    elemRowIdx = elemRowIdx+1;
                    % If next element exceeds the max elements per row, begin new row
                    % (Except for the last element)
                    if( elemRowIdx > maxElemPerRow && iP~=num_pars)
                        newRow();
                    end
                end
            end
            
            % Opens a filechooser and saves the result into the h_edit
            % edit field as well as options(structVarName)
            function callback_filechooserButton(uiObj, eventdata, h_edit, structVarName, pSettings)
                filterEnding = pSettings{2}; % file ending to filter for
                
                % Opens a file chooser dialog to select the dark movie
                lastPath = get(h_edit,'String');
                path = [];
                if ~isempty(lastPath)
                    [path,~,~] = fileparts(lastPath);
                end
                if isempty(filterEnding)
                    [filename, path] = uigetfile([path,filesep]);
                else
                    [filename, path] = uigetfile([path,filesep,'*.',filterEnding]);
                end
                if( isfloat(filename)); return; end; % User pressed cancel.

                value = [path,filename];
                set(h_edit,'String',value);
                obj.options.(structVarName) = value; % Save the parameter value
            end
            
            function callback_saveOnChange(uiObj, eventdata, structVarName, pType, pSettings)
                switch pType
                    case 'float'
                        parValue = str2double(get(uiObj,'String'));
                        % Obey upper/lower limits
                        parValue = max(pSettings{2},parValue);
                        parValue = min(pSettings{3},parValue);
                        
                        set(uiObj, 'String', num2str(parValue)); % Synchronise text field
                        obj.options.(structVarName) = parValue;
                    case 'int'
                        parValue = str2double(get(uiObj,'String')); % get integer value
                        % Obey upper/lower limits
                        parValue = max(pSettings{2},parValue);
                        parValue = min(pSettings{3},parValue);
                        % Round to int
                        parValue = round(parValue);
                        
                        set(uiObj, 'String', num2str(parValue)); % Synchronise text field
                        obj.options.(structVarName) = parValue;
                    case 'bool'
                        parValue = logical(get(uiObj,'Value'));
                        obj.options.(structVarName) = parValue;
                    case 'string'
                        parValue = get(uiObj,'String');
                        obj.options.(structVarName) = parValue;
                    case 'list'
                        choices = get(uiObj,'String');
                        optionString = choices{get(uiObj,'Value')};
                        obj.options.(structVarName) = optionString;
                    case 'filechooser'
                        parValue = get(uiObj,'String');
                        obj.options.(structVarName) = parValue;
                        
                    otherwise
                        error('createOptionsPanel:callback_saveOnChange: Unknown parameter type ''%s''!', pType);
                end
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
        
        
    end
    
end

