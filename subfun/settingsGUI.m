function [globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIreturns] = settingsGUI(globalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs)
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2015
%

% struct for communication of the outside world to the GUI
if nargin < 5 || isempty(GUIinputs)
    GUIinputs.titleText = '';
    GUIinputs.fileText = '';
    GUIinputs.singleFileMode = false;
end

% struct for communication of the GUI to the outside world
GUIreturns.useSettingsForAll = false;
GUIreturns.userExit = false;
GUIreturns.globalOptionsChanged = false;
GUIreturns.candidateOptionsChanged = false;
GUIreturns.fittingOptionsChanged = false;
GUIreturns.trackingOptionsChanged  = false;
GUIreturns.previewIntervalChanged = false;
GUIreturns.previewMode = false;

% Save options at startup (to check later if options changed)
globalOptions_atStartup = globalOptions;
candidateOptions_atStartup = candidateOptions;
fittingOptions_atStartup = fittingOptions;
trackingOptions_atStartup = trackingOptions;


% -- Preparing the GUI --
h_main = openfig('settingsGUI_Layout.fig');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % For cleanup

h_all = guihandles(h_main);

% Setup GUI specific elements
set(h_all.text_title, 'String', GUIinputs.titleText);
[~, filename, fileext] = fileparts(GUIinputs.fileText);
set(h_all.edit_title, 'String', [filename,fileext]);
set(h_all.button_save, 'Callback', @callback_saveSettings);
set(h_all.button_load, 'Callback', @callback_loadSettings);
set(h_all.button_continue, 'Callback',@callback_continue);
set(h_all.button_continueForAll, 'Callback', @callback_continueForAll);
set(h_all.button_preview, 'Callback', @callback_preview);

% % General options
% edit_darkMovie
set(h_all.edit_firstFrame,'Callback',{@callback_IntEdit,1,inf});
set(h_all.edit_lastFrame,'Callback',{@callback_IntEdit,1,inf});
set(h_all.edit_firstFrameTesting,'Callback',{@callback_IntEdit,1,inf});
set(h_all.edit_lastFrameTesting,'Callback',{@callback_IntEdit,1,inf});
set(h_all.button_addMovies, 'Callback', @callback_addMovies);
set(h_all.button_removeMovie, 'Callback', @callback_removeMovie);
set(h_all.button_darkMovie, 'Callback', @callback_selectDarkMovie);
% set(h_all.listbox_movieList, 'Callback', @callback_resetValueWhenEmpty);

%Photon conversion
set(h_all.cbx_usePhotonConv, 'Callback', @callback_updateMainGUIstate);
set(h_all.edit_photonBias, 'Callback', {@callback_IntEdit,0,inf});
set(h_all.edit_photonSensitivity, 'Callback', {@callback_FloatEdit,1.0,inf});
set(h_all.edit_photonGain, 'Callback', {@callback_IntEdit,1,1000});

% % Candidate plugin
set(h_all.popup_candidateMethod, 'Callback', @callback_updatePlugins);

% % Fitting plugin
set(h_all.popup_fittingMethod, 'Callback', @callback_updatePlugins);

% % Tracking plugin
set(h_all.cbx_enableTracking, 'Callback', @callback_updateMainGUIstate);
set(h_all.popup_trackingMethod, 'Callback', @callback_updatePlugins);
% cbx_verbose

% % Help buttons
set(h_all.button_candidateHelp, 'Callback', @callback_helpButtons);
set(h_all.button_fittingHelp, 'Callback', @callback_helpButtons);
set(h_all.button_trackingHelp, 'Callback', @callback_helpButtons);

% % Container for plugins
candidate_plugins = [];
fitting_plugins = [];
tracking_plugins = [];

% Save last selected plugin per category
selected_candidate_plugin = -1;
selected_fitting_plugin = -1;
selected_tracking_plugin = -1;

% Load the plugins
if ~GUIinputs.singleFileMode % Show warnings only on startup
    fprintf('TNT: Loading plugins ...\n')
    loadPlugins();
    fprintf('TNT: Successfully imported %i plugins (%i candidate detection, %i fitting, %i tracking).\n',numel(candidate_plugins)+numel(fitting_plugins)+numel(tracking_plugins),numel(candidate_plugins),numel(fitting_plugins),numel(tracking_plugins));
else
    warning off;
    loadPlugins();
    warning on;
end
movegui(h_main,'center');

% % GUI main
setGUIBasedOnOptions();
uiwait(h_main);
drawnow; % makes figure disappear instantly (otherwise it looks like it is existing until script finishes)


% Load plugins from the plugins folder
% Note: This loads the plugin data, but does not build their panels
    function loadPlugins()
        fullPathToThisFile = mfilename('fullpath');
        [folderPath,~,~] = fileparts(fullPathToThisFile);
        plugin_files = dir([folderPath filesep '..' filesep 'plugins' filesep 'plugin_*.m']);
        %         plugin_files = plugin_files(~cellfun(@(var) strcmp(var,'plugin_default.m'), {plugin_files(:).name}).');
        nr_plugins = numel(plugin_files);
        
        for iPlug = 1:nr_plugins
            try %% Try loading this plugin
                [~, plugin_constructor, ~] = fileparts(plugin_files(iPlug).name);
                plugin_constructor = str2func(plugin_constructor); % Convert string to function handle
                plugin = plugin_constructor();
                
                switch plugin.type
                    case 1 % Candidate method
                        candidate_plugins = [candidate_plugins, plugin];
                    case 2 % Fitting method
                        fitting_plugins  = [fitting_plugins, plugin];
                    case 3 % Tracking method
                        tracking_plugins = [tracking_plugins, plugin];
                    otherwise
                        warning('Detected unknown plugin of type %i',plugin.type);
                end
            catch err
                warning('TrackNTrace: Failed to load plugin file ''%s''. \n  Error: %s', plugin_files(iPlug).name, err.message);
            end
        end
        
        found_candidate_plugin = numel(candidate_plugins)>0;
        found_fitting_plugin = numel(fitting_plugins)>0;
        found_tracking_plugin = numel(tracking_plugins)>0;
        
        if not(found_candidate_plugin)
            error('No candidate detection plugin detected.');
        end
        if not(found_fitting_plugin)
            error('No fitting plugin detected.');
        end
        if not(found_tracking_plugin)
            error('No tracking plugin detected.')
        end
        
        % Set popup choices
        set(h_all.popup_candidateMethod, 'String', {candidate_plugins(:).name});
        set(h_all.popup_fittingMethod, 'String', {fitting_plugins(:).name});
        set(h_all.popup_trackingMethod, 'String', {tracking_plugins(:).name});
    end

% Select plugins based on the current candidateOptions, fittingOptions,
% trackingOptions. After selection, the plugin panels are constructed
    function selectPluginsBasedOnOptions()
        
        % Candidate detection
        if ~isempty(candidateOptions)
            selected_candidate_plugin = -1;
            % Search for a loaded plugin with the same name
            for iPlug = 1:numel(candidate_plugins)
                if strcmp(candidateOptions.plugin_name,candidate_plugins(iPlug).name)
                    selected_candidate_plugin = iPlug;
                end
            end
            
            if selected_candidate_plugin < 0
                error('Plugin ''%s'' not found.',candidateOptions.plugin_name) ;
            end
        else
            % Select the default plugin if it is found
            default_plugin_index = strfind({candidate_plugins(:).name},GUIinputs.TNToptions.defaultCandidatePlugin);
            default_plugin_index = find(~cellfun(@isempty,default_plugin_index));
            if isempty(default_plugin_index)
                warning off backtrace
                warning('Default candidate plugin ''%s'' not found. Selecting first one.',GUIinputs.TNToptions.defaultCandidatePlugin);
                warning on backtrace
                selected_candidate_plugin = 1;
            else
                selected_candidate_plugin = default_plugin_index;
            end
        end
        
        % Fitting
        if ~isempty(fittingOptions)
            selected_fitting_plugin = -1;
            % Search for a loaded plugin with the same name
            for iPlug = 1:numel(fitting_plugins)
                if strcmp(fittingOptions.plugin_name,fitting_plugins(iPlug).name)
                    selected_fitting_plugin = iPlug;
                end
            end
            
            if selected_fitting_plugin < 0
                error('Plugin ''%s'' not found.',fittingOptions.plugin_name) ;
            end
        else
            % Select the default plugin if it is found
            default_plugin_index = strfind({fitting_plugins(:).name},GUIinputs.TNToptions.defaultFittingPlugin);
            default_plugin_index = find(~cellfun(@isempty,default_plugin_index));
            if isempty(default_plugin_index)
                warning off backtrace
                warning('Default fitting plugin ''%s'' not found. Selecting first one.',GUIinputs.TNToptions.defaultFittingPlugin);
                warning on backtrace
                selected_fitting_plugin = 1;
            else
                selected_fitting_plugin = default_plugin_index;
            end
        end
        
        % Tracking
        if ~isempty(trackingOptions)
            selected_tracking_plugin = -1;
            % Search for a loaded plugin with the same name
            for iPlug = 1:numel(tracking_plugins)
                if strcmp(trackingOptions.plugin_name,tracking_plugins(iPlug).name)
                    selected_tracking_plugin = iPlug;
                end
            end
            
            if selected_tracking_plugin < 0
                error('Plugin ''%s'' not found.',trackingOptions.plugin_name) ;
            end
        else
            % Select the default plugin if it is found
            default_plugin_index = strfind({tracking_plugins(:).name},GUIinputs.TNToptions.defaultTrackingPlugin);
            default_plugin_index = find(~cellfun(@isempty,default_plugin_index));
            if isempty(default_plugin_index)
                warning off backtrace
                warning('Default tracking plugin ''%s'' not found. Selecting first one.',GUIinputs.TNToptions.defaultTrackingPlugin);
                warning on backtrace
                selected_tracking_plugin = 1;
            else
                selected_tracking_plugin = default_plugin_index;
            end
        end
        
        % Set popups to correct plugin name
        set(h_all.popup_candidateMethod, 'Value', selected_candidate_plugin);
        set(h_all.popup_fittingMethod, 'Value', selected_fitting_plugin);
        set(h_all.popup_trackingMethod, 'Value', selected_tracking_plugin);
        
        % Build panels by invoking the selected plugins function
        candidate_plugins( selected_candidate_plugin).setOptions(candidateOptions);
        fitting_plugins( selected_fitting_plugin).setOptions(fittingOptions);
        tracking_plugins( selected_tracking_plugin).setOptions(trackingOptions);
        
        % Create the panels
        candidate_plugins( selected_candidate_plugin).createOptionsPanel(h_all.panel_candidate);
        fitting_plugins( selected_fitting_plugin).createOptionsPanel(h_all.panel_fitting);
        tracking_plugins( selected_tracking_plugin).createOptionsPanel(h_all.panel_tracking);
        
        updatePanelPositions();
    end

% Update GUI based on currently set values
% Note: This does not update the plugins!
    function callback_updateMainGUIstate(hObj, event)
        
        % Enable/disable photon conversion
        if get(h_all.cbx_usePhotonConv, 'Value')
            set(h_all.edit_photonBias, 'Enable','on');
            set(h_all.edit_photonSensitivity, 'Enable','on');
            set(h_all.edit_photonGain, 'Enable','on');
        else
            set(h_all.edit_photonBias, 'Enable','off');
            set(h_all.edit_photonSensitivity, 'Enable','off');
            set(h_all.edit_photonGain, 'Enable','off');
        end
        
        % In single file mode we disable choosing the movie list
        if(GUIinputs.singleFileMode)
            set(h_all.listbox_movieList,'Enable', 'off');
            set(h_all.button_addMovies,'Enable', 'off');
            set(h_all.button_removeMovie,'Enable', 'off');
            %             set(h_all.button_continueForAll, 'Visible','off');
        else
            set(h_all.button_preview,'Enable', 'off');
            set(h_all.text_firstFrameTesting,'Enable', 'off');
            set(h_all.edit_firstFrameTesting,'Enable', 'off');
            set(h_all.text_lastFrameTesting,'Enable', 'off');
            set(h_all.edit_lastFrameTesting,'Enable', 'off');
        end
        
        % Enable/disable tracking panel
        if get(h_all.cbx_enableTracking, 'Value')
            set(h_all.popup_trackingMethod,'Enable','on');
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','on');
        else
            set(h_all.popup_trackingMethod,'Enable','off');
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','off');
        end
    end

% This function is called if the selection of plugins changes
% It saves the options for the previously selected plugin and builds
% the panel for the newly selected one.
    function callback_updatePlugins(hObj, event)
        selected_candidate_plugin =  get(h_all.popup_candidateMethod,'Value');
        selected_fitting_plugin = get(h_all.popup_fittingMethod,'Value');
        selected_tracking_plugin = get(h_all.popup_trackingMethod,'Value');
        
        % Update panels to display currently selected plugin
        candidate_plugins( selected_candidate_plugin).createOptionsPanel(h_all.panel_candidate);
        fitting_plugins( selected_fitting_plugin).createOptionsPanel(h_all.panel_fitting);
        tracking_plugins( selected_tracking_plugin).createOptionsPanel(h_all.panel_tracking);
        
        updatePanelPositions();
        
        % Enable/disable tracking panel
        if get(h_all.cbx_enableTracking, 'Value')
            set(h_all.popup_trackingMethod,'Enable','on');
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','on');
        else
            set(h_all.popup_trackingMethod,'Enable','off');
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','off');
        end
    end

% Used to resize the GUI after selecting a different plugin
    function updatePanelPositions()
        above_panel_spacing = 0.5;
        topic_spacing = 2.5;
        units = 'characters';
        
        set(h_all.panel_candidate,'Units',units);
        panel_candidate_pos = get(h_all.panel_candidate,'Position');
        
        % Relative to candidate detection panel
        set(h_all.label_fittingMethod, 'Units',units);
        label_fittingMethod_pos = get(h_all.label_fittingMethod, 'Position');
        label_fittingMethod_pos(2) = panel_candidate_pos(2)-topic_spacing-label_fittingMethod_pos(4)/2;
        set(h_all.label_fittingMethod, 'Position',label_fittingMethod_pos);
        
        set(h_all.popup_fittingMethod, 'Units',units);
        pos = get(h_all.popup_fittingMethod, 'Position');
        pos(2) = panel_candidate_pos(2)-topic_spacing-pos(4)/2;
        set(h_all.popup_fittingMethod, 'Position',pos);
        
        set(h_all.button_fittingHelp, 'Units',units);
        pos = get(h_all.button_fittingHelp, 'Position');
        pos(2) = panel_candidate_pos(2)-topic_spacing-pos(4)/2;
        set(h_all.button_fittingHelp, 'Position',pos);
        
        % Relative to label_fitting_Method
        set(h_all.panel_fitting, 'Units',units);
        panel_fitting_pos = get(h_all.panel_fitting, 'Position');
        panel_fitting_pos(2) = label_fittingMethod_pos(2)-above_panel_spacing-panel_fitting_pos(4);
        set(h_all.panel_fitting, 'Position',panel_fitting_pos);
        
        % Relative to fitting panel
        set(h_all.label_enableTracking, 'Units',units);
        label_enableTracking_pos = get(h_all.label_enableTracking, 'Position');
        label_enableTracking_pos(2) = panel_fitting_pos(2)-topic_spacing-label_enableTracking_pos(4)/2;
        set(h_all.label_enableTracking, 'Position',label_enableTracking_pos);
        
        set(h_all.cbx_enableTracking, 'Units',units);
        pos = get(h_all.cbx_enableTracking, 'Position');
        pos(2) = panel_fitting_pos(2)-topic_spacing-pos(4)/2;
        set(h_all.cbx_enableTracking, 'Position',pos);
        
        set(h_all.label_trackingMethod, 'Units',units);
        pos = get(h_all.label_trackingMethod, 'Position');
        pos(2) = panel_fitting_pos(2)-topic_spacing-pos(4)/2;
        set(h_all.label_trackingMethod, 'Position',pos);
        
        set(h_all.popup_trackingMethod, 'Units',units);
        pos = get(h_all.popup_trackingMethod, 'Position');
        pos(2) = panel_fitting_pos(2)-topic_spacing-pos(4)/2;
        set(h_all.popup_trackingMethod, 'Position',pos);
        
        set(h_all.button_trackingHelp, 'Units',units);
        pos = get(h_all.button_trackingHelp, 'Position');
        pos(2) = panel_fitting_pos(2)-topic_spacing-pos(4)/2;
        set(h_all.button_trackingHelp, 'Position',pos);
        
        % Relative to label tracking
        set(h_all.panel_tracking, 'Units',units);
        panel_tracking_pos = get(h_all.panel_tracking, 'Position');
        panel_tracking_pos(2) = label_enableTracking_pos(2)-above_panel_spacing-panel_tracking_pos(4);
        set(h_all.panel_tracking, 'Position',panel_tracking_pos);
        
        % Relative to tracking panel
        button_spacing = 1;
        
        set(h_all.button_continueForAll, 'Units',units);
        pos = get(h_all.button_continueForAll, 'Position');
        pos(2) = panel_tracking_pos(2)-button_spacing-pos(4);
        set(h_all.button_continueForAll, 'Position',pos);
        
        set(h_all.button_continue, 'Units',units);
        pos = get(h_all.button_continue, 'Position');
        pos(2) = panel_tracking_pos(2)-button_spacing-pos(4);
        set(h_all.button_continue, 'Position',pos);
        
        set(h_all.button_preview, 'Units',units);
        pos = get(h_all.button_preview, 'Position');
        pos(2) = panel_tracking_pos(2)-button_spacing-pos(4);
        set(h_all.button_preview, 'Position',pos);
        
        button_preview_height = pos(4);
        
        set(h_all.text_firstFrameTesting, 'Units',units);
        pos = get(h_all.text_firstFrameTesting, 'Position');
        pos(2) = panel_tracking_pos(2)-button_spacing-pos(4);
        set(h_all.text_firstFrameTesting, 'Position',pos);
        
        set(h_all.edit_firstFrameTesting, 'Units',units);
        pos = get(h_all.edit_firstFrameTesting, 'Position');
        pos(2) = panel_tracking_pos(2)-button_spacing-pos(4);
        set(h_all.edit_firstFrameTesting, 'Position',pos);
        
        set(h_all.text_lastFrameTesting, 'Units',units);
        pos = get(h_all.text_lastFrameTesting, 'Position');
        pos(2) = panel_tracking_pos(2)-button_spacing-button_preview_height;
        set(h_all.text_lastFrameTesting, 'Position',pos);
        
        set(h_all.edit_lastFrameTesting, 'Units',units);
        pos = get(h_all.edit_lastFrameTesting, 'Position');
        pos(2) = panel_tracking_pos(2)-button_spacing-button_preview_height;
        set(h_all.edit_lastFrameTesting, 'Position',pos);
        
        % Rescale window based on last element
        set(h_main,'Units',units);
        win_pos = get(h_main,'Position');
        diff_height = pos(2)-button_spacing;
        
        % To resize the figure properly, we first need to move all objects
        % inside.. (Matlab ..)
        all_uiObjects = get(h_main,'Children');
        
        for iObj = 1:numel(all_uiObjects)
            set(all_uiObjects(iObj),'Units',units);
            pos = get(all_uiObjects(iObj),'Position');
            pos(2) = pos(2) - diff_height;
            set(all_uiObjects(iObj),'Position',pos);
        end
        
        win_pos(4) = win_pos(4)-diff_height;
        set(h_main,'Position', win_pos);
    end


% This will process all given movies with the current settings without
% individual figures showing up for each one.
    function callback_continueForAll(hObj,event)
        GUIreturns.useSettingsForAll = true;
        callback_continue(hObj,event);
    end

% Save the current settings to file
    function callback_saveSettings(hObj, event)
        storeOptions();
        
        [outfile,path] = uiputfile('.mat');
        if isfloat(outfile); return; end; % User clicked cancel
        
        outfile = [path,outfile];
        save(outfile,'globalOptions', 'candidateOptions','fittingOptions','trackingOptions');
    end

% Load settings from a file
    function callback_loadSettings(hObj,event)
        % In single file mode, the list of movies to process should not be
        % changed, we store it and restore after loading the settings
        if(GUIinputs.singleFileMode)
            filename_movies = globalOptions.filename_movies;
        end
        
        [infile, path] = uigetfile('.mat');
        if isfloat(infile);
            return;
        end; % User clicked cancel
        % Note: Loading has to be done this way, as variables "can not be
        % added to a static workspace" (e.g. the one of this GUI).
        allOptions = load([path,infile],'globalOptions', 'candidateOptions','fittingOptions','trackingOptions');
        globalOptions   = allOptions.globalOptions;
        if(GUIinputs.singleFileMode)
            globalOptions.filename_movies = filename_movies;
        end
        candidateOptions = allOptions.candidateOptions;
        fittingOptions   = allOptions.fittingOptions;
        trackingOptions  = allOptions.trackingOptions;
        
        setGUIBasedOnOptions();
    end

% Opens a file chooser dialog to choose multiple input (movie) files
% for processing. Note: filenames will be seperated by ';'
    function callback_addMovies(hObj,event)
        % Get current text field to set starting path of uigetfile
        listbox_entries = get(h_all.listbox_movieList,'String');
        
        % We have to do this, since repopulating a listbox does not
        % automatically reset its value..
        if numel(listbox_entries) == 0
            set(h_all.listbox_movieList,'Value',1);
        end
        
        path = [];
        if ~isempty(listbox_entries)
            [path,~,~] = fileparts(listbox_entries{end});
        end
        [movieList, path] = uigetfile([path,filesep,'*.tif'],'MultiSelect','on');
        if( isfloat(movieList) ); return; end; % User pressed cancel.
        
        % atach path to every entry in the list then add to listbox
        if iscell(movieList)
            for i=1:length(movieList)
                movieList{i} =  [path,movieList{i}];
            end
            % Add to listbox
            listbox_entries = [listbox_entries; movieList.'];
        elseif ischar(movieList)
            movieList = [path,movieList];
            % Add to listbox
            listbox_entries = [listbox_entries; {movieList}];
        end
        
        set(h_all.listbox_movieList,'String',listbox_entries);
    end

% Remove a movie from the movie list
    function callback_removeMovie(hObj,event)
        selected_entry = get(h_all.listbox_movieList,'Value');
        listbox_entries = get(h_all.listbox_movieList,'String');
        
        % When listbox is empty, do nothing
        if numel(listbox_entries) == 0
           return; 
        end
        
        % When last selected item is deleted, select the one before it
        if selected_entry == numel(listbox_entries)
            set(h_all.listbox_movieList,'Value',selected_entry-1);
        end
        
        listbox_entries(selected_entry) = [];
        set(h_all.listbox_movieList,'String',listbox_entries);
    end

% Opens a file chooser dialog to select the dark movie
    function callback_selectDarkMovie(hObj, event)
        darkMoviePath = get(h_all.edit_darkMovie,'String');
        path = [];
        if ~isempty(darkMoviePath)
            [path,~,~] = fileparts(darkMoviePath);
        end
        [darkMovie, path] = uigetfile([path,filesep,'*.tif']);
        if( isfloat(darkMovie)); return; end; % User pressed cancel.
        
        set(h_all.edit_darkMovie,'String',[path,darkMovie]);
    end

% Shows a message dialog with the info of the currently selected plugin
    function callback_helpButtons(hObj, event)
        switch hObj
            case h_all.button_candidateHelp
                name = candidate_plugins(selected_candidate_plugin).name;
                msg = candidate_plugins(selected_candidate_plugin).info;
            case h_all.button_fittingHelp
                name = fitting_plugins(selected_fitting_plugin).name;
                msg = fitting_plugins(selected_fitting_plugin).info;
            case h_all.button_trackingHelp
                name = tracking_plugins(selected_tracking_plugin).name;
                msg = tracking_plugins(selected_tracking_plugin).info;
            otherwise
                error('Unknown caller of callback. Something went wront =?!');
        end
        
        % Get title
        title = ['About ', name];
        % Show message box, the sprintf allows using commands like '\n' for a line break
        % in the message.
        msgbox(sprintf(msg), title, 'help','non-modal');
    end

% Callback for edit fields containing floats. Checks if a correct
% number was entered and restricts it to the given bounds.
    function callback_FloatEdit(hObj,event, minVal, maxVal)
        if nargin<3 || isempty(minVal);
            minVal=-inf;
        end
        if nargin<4 || isempty(maxVal);
            maxVal=inf;
        end
        
        value = str2num(get(hObj, 'String'));
        if isempty(value)
            set(hObj,'ForegroundColor','r');
            set(hObj,'String','INVALID');
            uicontrol(hObj);
        else
            value = max(minVal,value);
            value = min(maxVal,value);
            set(hObj,'ForegroundColor','k');
            set(hObj,'String',sprintf('%.2f',value));
        end
    end

% Callback for edit fields containing integer values. Checks if a correct
% number was entered and restricts it to the given bounds.
    function callback_IntEdit(hObj,event, minVal,maxVal)
        if nargin<3 || isempty(minVal);
            minVal=0;
        end
        if nargin<4 || isempty(maxVal);
            maxVal=inf;
        end
        
        value = round(str2num(get(hObj,'String')));
        if isempty(value)
            set(hObj,'ForegroundColor','r');
            set(hObj,'String','INVALID');
            uicontrol(hObj);
        else
            value = max(minVal,value);
            value = min(maxVal,value);
            set(hObj,'ForegroundColor','k');
            set(hObj,'String',sprintf('%i',value));
        end
    end

% Takes a cell array containing strings and concatenates them into one
% string seperated by ';'. If the input is not a cell but a string, the
% output is equal to the input string;
    function str = cell2str(cellObj, delimiter)
        if nargin<2
            delimiter = ';';
        end
        
        str = '';
        if(iscell(cellObj))
            str = cellObj{1};
            for i=2:length(cellObj)
                str = [str,delimiter,cellObj{i}];
            end
        elseif ischar(cellObj)
            str = cellObj;
        end
    end

% Gets the 'String' property from a uicontrol and returns the string
% converted to a number
    function value = getNum(hObj)
        value = str2double(get(hObj,'String'));
    end

% Sets the 'String' of an edit field to the input 'value'. Set
% 'isInteger' to true for setting/displaying integer values.
    function setNum(hObj,value,isInteger)
        %         value = num2str(value);
        if nargin<3 || isempty(isInteger)
            isInteger = false;
        end
        
        if isInteger
            set(hObj,'String',sprintf('%i',round(value)));
        else
            set(hObj,'String',sprintf('%.2f',value));
        end
    end

% Sets a popup-menu's 'value' to the choice defined by 'optionString'.
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

% Gets the string of a popup-menu matching the currently selected value.
    function optionString = getPopup(hObj)
        choices = get(hObj,'String');
        optionString = choices{get(hObj,'Value')};
    end

% Set all UI fields based on the current value of the options structs
% (globalOptions, candidateOptions, fittingOptions, trackingOptions)
    function setGUIBasedOnOptions()
        % % General Options
        
        % -- Listbox movie list --
        % We have to reset the listbox value, since repopulating a listbox does not do it automatically..
        if numel(get(h_all.listbox_movieList,'String')) == 0
            set(h_all.listbox_movieList,'Value',1);
        end
        % Convert single file name 'filepath' to cell. Usually filename_movies is always a cell, but the user can set a (char) path 'path' for the default global options.
        if (isempty(globalOptions.filename_movies)) 
            globalOptions.filename_movies = {};
        elseif (ischar(globalOptions.filename_movies))
            globalOptions.filename_movies = {globalOptions.filename_movies}; 
        end       
        set(h_all.listbox_movieList,'String', globalOptions.filename_movies);
        % -- ----------------- --
        
        set(h_all.edit_darkMovie,'String', globalOptions.filename_dark_movie);
        setNum(h_all.edit_firstFrame, globalOptions.firstFrame, true);
        setNum(h_all.edit_lastFrame, globalOptions.lastFrame, true);
        setNum(h_all.edit_firstFrameTesting, globalOptions.firstFrameTesting, true);
        setNum(h_all.edit_lastFrameTesting, globalOptions.lastFrameTesting, true);
        
        % Photon conversion
        set(h_all.cbx_usePhotonConv,'Value',globalOptions.usePhotonConversion)
        setNum(h_all.edit_photonBias,globalOptions.photonBias, true);
        setNum(h_all.edit_photonSensitivity,globalOptions.photonSensitivity);
        setNum(h_all.edit_photonGain,globalOptions.photonGain, true);
        
        % % Tracking
        set(h_all.cbx_enableTracking,'Value', globalOptions.enableTracking); % --
        
        % Update plugins
        selectPluginsBasedOnOptions();
        
        % Update the GUI
        callback_updateMainGUIstate();
    end

% Set all options structs (globalOptions, candidateOptions, fittingOptions, trackingOptions)
% based on the current state of the uicontrols. Also checks which
% options structs have changed compared to the startup values.
    function storeOptions()
        % % General Options
        globalOptions.filename_movies = get(h_all.listbox_movieList,'String');
        globalOptions.filename_dark_movie = get(h_all.edit_darkMovie,'String');
        globalOptions.firstFrame = getNum(h_all.edit_firstFrame);
        globalOptions.lastFrame =  getNum(h_all.edit_lastFrame);
        globalOptions.firstFrameTesting = getNum(h_all.edit_firstFrameTesting);
        globalOptions.lastFrameTesting = getNum(h_all.edit_lastFrameTesting);
        
        globalOptions.usePhotonConversion = logical(get(h_all.cbx_usePhotonConv,'Value'));
        globalOptions.photonBias = getNum(h_all.edit_photonBias);
        globalOptions.photonSensitivity = getNum(h_all.edit_photonSensitivity);
        globalOptions.photonGain = getNum(h_all.edit_photonGain);
        
        % Tracking
        globalOptions.enableTracking = logical(get(h_all.cbx_enableTracking,'Value')); % --
        
        % Store options from plugins
        candidateOptions = candidate_plugins(selected_candidate_plugin).getOptions();
        fittingOptions = fitting_plugins(selected_fitting_plugin).getOptions();
        trackingOptions = tracking_plugins(selected_tracking_plugin).getOptions();
        
        %Check if options were changed compared to intial ones
        GUIreturns.globalOptionsChanged   = ~isequaln(globalOptions_atStartup, globalOptions);
        GUIreturns.candidateOptionsChanged = ~isequaln(candidateOptions_atStartup, candidateOptions);
        GUIreturns.fittingOptionsChanged   = ~isequaln(fittingOptions_atStartup, fittingOptions);
        GUIreturns.trackingOptionsChanged  = ~isequaln(trackingOptions_atStartup, trackingOptions);
        %Check if preview window changed
        GUIreturns.previewIntervalChanged = (globalOptions_atStartup.firstFrameTesting ~= globalOptions.firstFrameTesting) || (globalOptions_atStartup.lastFrameTesting ~= globalOptions.lastFrameTesting);
        
        % Check if only change in globalOptions is enableTracking
        GUIreturns.globalOptionsChanged_ExcludingEnableTracking = false;
        if GUIreturns.globalOptionsChanged
            globalOptions_tmp = globalOptions;
            globalOptions_tmp.enableTracking = globalOptions_atStartup.enableTracking;
            
            GUIreturns.globalOptionsChanged_ExcludingEnableTracking = ~isequaln(globalOptions_atStartup, globalOptions_tmp);
        end
    end



    function callback_continue(hObj, event)
        storeOptions(); % Store the options before closing
        delete(h_main);
    end

% Cleanup function. This is neccessary to delete the timer!
    function callback_preview(hObj, event)
        storeOptions(); % Store the options before closing
        GUIreturns.previewMode = true;
        delete(h_main);
    end

% Called when closing the application via the 'X' button (or via close)
    function onAppClose(hObj, event)
        storeOptions();
        GUIreturns.userExit = true;
        delete(h_main);
    end

end
