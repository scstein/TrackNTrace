% Loads all import plugins and returns a list of readable formats in
% the syntax of uigetfile.
% Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de, 2020
function [import_plugins,formats] = loadImportPlugins()
    fullPathToThisFile = mfilename('fullpath');
    [path,~,~] = fileparts(fileparts(fullPathToThisFile));

    plugin_path = [path filesep 'plugins'];
    if ~isonpath(plugin_path)
        addpath(plugin_path);
    end
    
    plugin_files = dir([plugin_path filesep 'plugin_*.m']);
    nr_plugins = numel(plugin_files);

    import_plugins = [];

    for iPlug = 1:nr_plugins
        try %% Try loading this plugin
            [~, plugin_constructor, ~] = fileparts(plugin_files(iPlug).name);
            plugin_constructor = str2func(plugin_constructor); % Convert string to function handle
            plugin = plugin_constructor();

            switch plugin.type
                case {1,2,3,4,6}
                    % Ignore
                case 5 % Import
                    import_plugins = [import_plugins, plugin];
                otherwise
                    warning('Detected unknown plugin of type %i',plugin.type);
            end
        catch err
            warning('TrackNTrace: Failed to load plugin file ''%s''. \n  Error: %s', plugin_files(iPlug).name, err.message);
        end
    end

    found_import_plugin = numel(import_plugins)>0;

    if not(found_import_plugin)
        error('No file import plugin detected.');
    end
    if nargout>1
        formats = {};
        for iPlug = 1:numel(import_plugins)
            formats = [formats; import_plugins(iPlug).info.supportedFormats]; %#ok<AGROW>
        end
        formats = [{strjoin(formats(:,1),';'), 'Compatible formats'};formats];
        %{'*.tif;*.ptu','Compatible formats';'*.tif','TIFF';'*.ptu','PTU'}
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