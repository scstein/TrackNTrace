function [ param_specification ] = add_plugin_param( param_specification, par_name, par_type, par_settings, par_tooltip )
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
             error('add_plugin_param: Failed adding param ''%s'' of type ''%s''. Settings need to be {defaultValue, lowerBound, upperBound)!', par_name, par_type) 
         end
    case 'int'
        if(~iscell(par_settings) || length(par_settings) ~= 3)
             error('add_plugin_param: Failed adding param ''%s'' of type ''%s''. Settings need to be {defaultValue, lowerBound, upperBound)!', par_name, par_type) 
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
    otherwise
        error('add_plugin_param: Unknown parameter type ''%s'' of parameter ''%s''.', par_type, par_name);
end
        
% --- Add parameter information ---
param_specification = vertcat(param_specification, {par_name, par_type, par_settings, par_tooltip});
end

