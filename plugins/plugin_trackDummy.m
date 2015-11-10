function [plugin_name, plugin_type] = plugin_trackDummy(h_panel, inputOptions)
    if nargin < 2
        inputOptions = [];
    end

    % Name of the component these options are for
    plugin_name = 'Placeholder tracking';

    % Type of plugin.
    % 1: Candidate detection
    % 2: Spot fitting
    % 3: Tracking
    plugin_type = 3;
    
    % Calling the plugin function without arguments just returns its name and type
    if (nargin == 0); return; end
    
    % Enter names of the parameters
    % These translate to the names of variables inside options struct this plugin
    % outputs by removing all white spaces.
    par_name  = {'Param A','Param B'};

    % Enter type of the parameters
    % possible: 'float', 'int', 'bool','list'
    par_type  = {'bool','float'};

    % Default value for parameters
    % Should be a number for 'float'/'int', true/false for 'bool'
    % or a cell array string list of possible choices for 'list' (first entry is default)
    par_defaultValue = {true,3.2};

    % Tooltip for the parameters
    par_tooltip = {'Tooltip A','Tooltip B'};

    createOptionsPanel(h_panel, plugin_name, par_name, par_type, par_defaultValue, par_tooltip,inputOptions);

    % Save handle of the plugins function
    options = getappdata(h_panel,'options');
    options.functionHandle = @() 1;
    setappdata(h_panel,'options',options);
end


