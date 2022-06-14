function [plugin] = plugin_tracking_prefilter()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Filter and ...';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
type = 3;

% The functions this plugin implements
mainFunc =  @filterAndRunTracking;

% Description of output parameters
outParamDescription = {' '}; % is set in mainFunc adapting to input data

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Function to execute before frame-by-frame processing starts
plugin.initFunc = @initFunc;
plugin.postFunc = @postFunc;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Applies the given filter to the refined localisations and passes the remaining localisations to the selected tracking plugin.\n',...
    'Use the TNTVisulizer to test filters.'];

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('Filter',...
    'string',...
    'Amp_Peak_>0',...
    'Filter string to apply to all localisations.');
plugin.add_param('Tracking plugin name',...
    'list',...
    @getTrackingPlugins,... % It is important that getTrackingPlugins is only called when the panel is built and not already in this constructor as this would lead to recursive calls.
    'Select tracking plugin.');
plugin.add_param('Select',...
    'button',...
    @getTrackingPluginOptions,...
    'Confirm selection to set the tracking plugin options.');
plugin.newRow();
plugin.add_text('Tracking plugin options');
plugin.newRow();

end


%   -------------- User functions --------------

function [trackingData,trackingOptions] = filterAndRunTracking(refinementData,trackingOptions)
    % Filter refinement data (code from TNTVisulizer)
    global refinementOptions % Needed as we drag along names specified here
    getFilterFun = @(paramsNames)str2func(['@(posData)true(size(posData,1),1)&' regexprep(vectorize(['( ' trackingOptions.Filter ' )']),...
        strcat('(?<!\w)',matlab.lang.makeValidName(paramsNames),'(?!\w)'),... % Replace spaces with underscore and makes sure the match is not within a name.
        cellfun(@(n)sprintf('posData(:,%i)',n),num2cell(1:numel(paramsNames)),'UniformOutput',false),...
        'ignorecase')]);
    
    try
        empty_ind = cellfun(@numel,refinementData)==0;
        filter_ind = cellfun(getFilterFun(refinementOptions.outParamDescription(:)),refinementData(~empty_ind),'UniformOutput',false);
        refinementData(~empty_ind) = cellfun(@(d,ind)d(ind,:),refinementData(~empty_ind),filter_ind,'UniformOutput',false);

    catch err
        disp( getReport( err, 'extended', 'hyperlinks', 'on' ) );
        warning('Error when applying filter: ''%s''. Proceeding with unfiltered localisations.',trackingOptions.Filter);
    end
        
    % Pass to selected tracking plugin
    switch nargout(trackingOptions.trackingPlugin.mainFunc)
        case 1 % Only tracking data assigned during call
            trackingData = trackingOptions.trackingPlugin.mainFunc(refinementData,trackingOptions);
        case 2 % Tracking data assigned + options changed
            [trackingData, trackingOptions] = trackingOptions.trackingPlugin.mainFunc(refinementData,trackingOptions);
        otherwise
            warning('Too many output arguments in mainFunc of tracking plugin ''%s''. Ignoring additional outputs', trackingOptions.plugin_name);
    end
end

function trackingOptions = initFunc(trackingOptions)
    if ~isempty(trackingOptions.trackingPlugin.initFunc)
        trackingOptions = trackingOptions.initFunc(trackingOptions);
    end
end
function [trackingData,trackingOptions] = postFunc(trackingData,trackingOptions)
    if ~isempty(trackingOptions.trackingPlugin.postFunc)
        [trackingData,trackingOptions] = trackingOptions.trackingPlugin.postFunc(trackingData,trackingOptions);
    end
end
%% Callback
function [opts,obj] = getTrackingPluginOptions(opts,obj)
    tracking_plugins = loadTrackingPlugins();
    pluginNames = {tracking_plugins.name};
    selectedPlugin = tracking_plugins(strcmp(opts.Tracking_plugin_name,pluginNames));
    obj.param_specification = [obj.param_specification(1:5,:);selectedPlugin.param_specification];
    opts.trackingPlugin = selectedPlugin;
    if isempty(selectedPlugin.param_specification)
        obj.param_specification{5,3}{1} = ['''' selectedPlugin.name ''' has no options'];        
    else
        obj.param_specification{5,3}{1} = [selectedPlugin.name ' Options'];        
    end
end
function pluginNames = getTrackingPlugins(opts)
    tracking_plugins = loadTrackingPlugins();
    pluginNames = {tracking_plugins.name};
    pluginNames = pluginNames(~strcmpi(pluginNames,opts.plugin_name));
end
%%
% Loads all tracking plugins
function [tracking_plugins] = loadTrackingPlugins()
    fullPathToThisFile = mfilename('fullpath');
    [path,~,~] = fileparts(fileparts(fullPathToThisFile));

    plugin_path = [path filesep 'plugins'];
    if ~isonpath(plugin_path)
        addpath(plugin_path);
    end
    
    plugin_files = dir([plugin_path filesep 'plugin_*.m']);
    nr_plugins = numel(plugin_files);

    tracking_plugins = [];

    for iPlug = 1:nr_plugins
        try %% Try loading this plugin
            [~, plugin_constructor, ~] = fileparts(plugin_files(iPlug).name);
            plugin_constructor = str2func(plugin_constructor); % Convert string to function handle
            plugin = plugin_constructor();

            switch plugin.type
                case 3 % Tracking
                    tracking_plugins = [tracking_plugins, plugin];
                otherwise
            end
        catch err
            warning('TrackNTrace: Failed to load plugin file ''%s''. \n  Error: %s', plugin_files(iPlug).name, err.message);
        end
    end

    found_tracking_plugin = numel(tracking_plugins)>0;

    if not(found_tracking_plugin)
        error('No file tracking plugin detected.');
    end
end

function onPath = isonpath(Folder)
pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(Folder, pathCell));
else
  onPath = any(strcmp(Folder, pathCell));
end
end