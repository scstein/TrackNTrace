function plugin_crossCorrelation()

% Name of the component these options are for
plugin_name = 'Cross correlation';

% Enter names of the parameters
% These translate to the names of variables inside options struct this plugin
% outputs by removing all white spaces.
par_name  = {'Param a','bool param','NiceParam List','Name'};

% Enter type of the parameters
% possible: 'float', 'int', 'bool','list'
par_type  = {'int','bool','list','float'};

% Default value for parameters
% Should be a number for 'float'/'int', true/false for 'bool'
% or a cell array string list of possible choices for 'list' (first entry is default)
par_defaultValue = {1.263,true,{'Option A','Option B','Option Dideldumdai'},3.2};

% Tooltip for the parameters
par_tooltip = {'This is a very long sentence for a tooltip.', 'does', 'help','blub'};

h = createOptionsGUI(plugin_name, par_name, par_type, par_defaultValue, par_tooltip); %, inputOptions);
%options = getappdata(h,'options')

options = getappdata(h,'options');
options.functionHandle = @crossCorrelation;
setappdata(h,'options',options);



    function crossCorrelation(img,options)
       disp('h'); 
       a = options;
    end

end

