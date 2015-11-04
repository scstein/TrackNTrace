% Name of the component these options are for
name = 'Cross correlation';

% Enter names of the parameters
par_text  = {'param a','bool param','NiceParam List','name'};

% Enter type of the parameters
% possible: 'float', 'int', 'bool','list'
par_type  = {'float','bool','list','float'};

% Default value for parameters
% Should be a number for 'float'/'int', true/false for 'bool'
% or a cell array string list of possible choices for 'list' (first entry is default)
par_value = {1.263,true,{'Option A','Option B','Option Jan'},3.2};

% Tooltip for the parameters
par_tooltip = {'This is a very long sentence for a tooltip.', 'does', 'help','blub'};

h = createOptionsGUI(name, par_text, par_type, par_value, par_tooltip);