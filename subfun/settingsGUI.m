% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
%     Copyright (C) 2016  Simon Christoph Stein, scstein@phys.uni-goettingen.de
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [globalOptions, candidateOptions,refinementOptions,trackingOptions, GUIreturns] = settingsGUI(globalOptions, candidateOptions,refinementOptions,trackingOptions, GUIinputs)
% TrackNTrace main GUI. Here plugins can be selected, their options 
% adjusted. This GUI is used for this management only, while the overall
% program flow is handled by RunTrackNTrace.m
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2016
%

titleText = 'Adjust options for movie:';
filename_movie = GUIinputs.filename_movie; % The movie we edit settings for

% struct for communication of the GUI to the outside world
GUIreturns.useSettingsForAll = false;
GUIreturns.userExit = false;
GUIreturns.globalOptionsChanged = false;
GUIreturns.candidateOptionsChanged = false;
GUIreturns.refinementOptionsChanged = false;
GUIreturns.trackingOptionsChanged  = false;
GUIreturns.previewIntervalChanged = false;
GUIreturns.previewMode = false;
GUIreturns.outputFolderSameAsMovie = false;
GUIreturns.outputFolder = [];

% Save options at startup (to check later if options changed)
globalOptions_atStartup = globalOptions;
candidateOptions_atStartup = candidateOptions;
refinementOptions_atStartup = refinementOptions;
trackingOptions_atStartup = trackingOptions;


% -- Preparing the GUI --
h_main = openfig('settingsGUI_Layout.fig');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % For cleanup

h_all = guihandles(h_main);

% Unix specific
if isunix
        uiObjects = get(h_main,'Children');        
        for iO = 1:numel(uiObjects)
            set(uiObjects(iO),'FontSize',8);
        end
        
        uiObjects = get(h_all.panel_general,'Children');        
        for iO = 1:numel(uiObjects)
            set(uiObjects(iO),'FontSize',8);
        end
end



% Setup GUI specific elements at top of GUI
set(h_all.text_title, 'String', titleText);
[~, fname,extension] = fileparts(filename_movie);
set(h_all.edit_title, 'String', [fname,extension]);
set(h_all.edit_title,'Tooltip',filename_movie);
set(h_all.button_save, 'Callback', @callback_saveSettings);
set(h_all.button_load, 'Callback', @callback_loadSettings);


% % General options
set(h_all.cbx_outputFolder, 'Callback', @callback_updateMainGUIstate);
set(h_all.cbx_outputFolder, 'Value', GUIinputs.outputFolderSameAsMovie);
set(h_all.button_outputFolder, 'Callback', @callback_selectOutputFolder);
set(h_all.edit_outputFolder, 'String', GUIinputs.outputFolder);
set(h_all.edit_outputFolder, 'Tooltip', GUIinputs.outputFolder);

set(h_all.button_darkMovie, 'Callback', @callback_selectDarkMovie);
% edit_darkMovie

set(h_all.edit_firstFrame,'Callback',{@callback_IntEdit,1,inf});
set(h_all.edit_lastFrame,'Callback',{@callback_IntEdit,1,inf});
%Photon conversion
set(h_all.cbx_usePhotonConv, 'Callback', @callback_updateMainGUIstate);
set(h_all.edit_photonBias, 'Callback', {@callback_IntEdit,0,inf});
set(h_all.edit_photonSensitivity, 'Callback', {@callback_FloatEdit,1.0,inf});
set(h_all.edit_photonGain, 'Callback', {@callback_IntEdit,1,1000});

% % Candidate plugin
set(h_all.popup_candidateMethod, 'Callback', @callback_updatePlugins);
set(h_all.cbx_candidate_loaded, 'Callback', @callback_updateMainGUIstate_loadedData);
if(GUIinputs.show_candidateData_fromFile_cbx)
 set(h_all.cbx_candidate_loaded, 'Value', GUIinputs.use_loaded_candidateData);
else
 set(h_all.cbx_candidate_loaded, 'Visible', 'off');
end

% % Refinement plugin
set(h_all.popup_refinementMethod, 'Callback', @callback_updatePlugins);
set(h_all.cbx_refinement_loaded, 'Callback', @callback_updateMainGUIstate_loadedData);
if(GUIinputs.show_refinementData_fromFile_cbx)
 set(h_all.cbx_refinement_loaded, 'Value', GUIinputs.use_loaded_refinementData);
else
 set(h_all.cbx_refinement_loaded, 'Visible', 'off');
end

% % Tracking plugin
set(h_all.cbx_enableTracking, 'Callback', @callback_updateMainGUIstate);
set(h_all.popup_trackingMethod, 'Callback', @callback_updatePlugins);
% cbx_verbose

% % Help buttons
set(h_all.button_candidateHelp, 'Callback', @callback_helpButtons);
set(h_all.button_refinementHelp, 'Callback', @callback_helpButtons);
set(h_all.button_trackingHelp, 'Callback', @callback_helpButtons);

% % Elements at bottom of GUI
set(h_all.button_preview, 'Callback', @callback_preview);
set(h_all.edit_firstFrameTesting,'Callback',{@callback_IntEdit,1,inf});
set(h_all.edit_lastFrameTesting,'Callback',{@callback_IntEdit,1,inf});
set(h_all.button_continueForAll, 'Callback', @callback_continueForAll);

set(h_all.button_continue, 'Callback',@callback_continue);
if(GUIinputs.lastMovieInList)
    set(h_all.button_continue, 'String','<html> <center> START <br> processing </center></html>');
    set(h_all.button_continue, 'ForegroundColor',[16,33,130]/255);
    set(h_all.button_continue, 'Tooltip','Start processing all movies.');
else
    set(h_all.button_continue, 'String','Next movie');
    set(h_all.button_continue, 'ForegroundColor',[0,0,0]);
    set(h_all.button_continue, 'Tooltip','Go to next movie in list to adjust its settings.');
end


% % Container for plugins
candidate_plugins = [];
refinement_plugins = [];
tracking_plugins = [];

% Save last selected plugin per category
selected_candidate_plugin = -1;
selected_refinement_plugin = -1;
selected_tracking_plugin = -1;

% Load the plugins
if GUIinputs.showStartupInformation % Show warnings only on startup
    fprintf('TNT: Loading plugins ...\n')
    loadPlugins();
    fprintf('TNT: Successfully imported %i plugins (%i candidate detection, %i refinement, %i tracking).\n',numel(candidate_plugins)+numel(refinement_plugins)+numel(tracking_plugins),numel(candidate_plugins),numel(refinement_plugins),numel(tracking_plugins));
else
    warning off;
    loadPlugins();
    warning on;
end
movegui(h_main,'center');

% % GUI main
setGUIBasedOnOptions();
callback_updateMainGUIstate_loadedData();

uiwait(h_main);
drawnow; % makes figure disappear instantly (otherwise it looks like it is existing until script finishes)


%% -------- Nested functions --------

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
                    case 1 % Candidate search
                        candidate_plugins = [candidate_plugins, plugin];
                    case 2 % Refinement/fitting
                        refinement_plugins  = [refinement_plugins, plugin];
                    case 3 % Tracking
                        tracking_plugins = [tracking_plugins, plugin];
                    otherwise
                        warning('Detected unknown plugin of type %i',plugin.type);
                end
            catch err
                warning('TrackNTrace: Failed to load plugin file ''%s''. \n  Error: %s', plugin_files(iPlug).name, err.message);
            end
        end
        
        found_candidate_plugin = numel(candidate_plugins)>0;
        found_refinement_plugin = numel(refinement_plugins)>0;
        found_tracking_plugin = numel(tracking_plugins)>0;
        
        if not(found_candidate_plugin)
            error('No candidate detection plugin detected.');
        end
        if not(found_refinement_plugin)
            error('No refinement plugin detected.');
        end
        if not(found_tracking_plugin)
            error('No tracking plugin detected.')
        end
        
        % Set popup choices
        set(h_all.popup_candidateMethod, 'String', {candidate_plugins(:).name});
        set(h_all.popup_refinementMethod, 'String', {refinement_plugins(:).name});
        set(h_all.popup_trackingMethod, 'String', {tracking_plugins(:).name});
    end

% Select plugins based on the current candidateOptions, refinementOptions,
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
        
        % Refinement
        if ~isempty(refinementOptions)
            selected_refinement_plugin = -1;
            % Search for a loaded plugin with the same name
            for iPlug = 1:numel(refinement_plugins)
                if strcmp(refinementOptions.plugin_name,refinement_plugins(iPlug).name)
                    selected_refinement_plugin = iPlug;
                end
            end
            
            if selected_refinement_plugin < 0
                error('Plugin ''%s'' not found.',refinementOptions.plugin_name) ;
            end
        else
            % Select the default plugin if it is found
            default_plugin_index = strfind({refinement_plugins(:).name},GUIinputs.TNToptions.defaultRefinementPlugin);
            default_plugin_index = find(~cellfun(@isempty,default_plugin_index));
            if isempty(default_plugin_index)
                warning off backtrace
                warning('Default fitting plugin ''%s'' not found. Selecting first one.',GUIinputs.TNToptions.defaultRefinementPlugin);
                warning on backtrace
                selected_refinement_plugin = 1;
            else
                selected_refinement_plugin = default_plugin_index;
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
        set(h_all.popup_refinementMethod, 'Value', selected_refinement_plugin);
        set(h_all.popup_trackingMethod, 'Value', selected_tracking_plugin);
        
        % Build panels by invoking the selected plugins function
        candidate_plugins( selected_candidate_plugin).setOptions(candidateOptions);
        refinement_plugins( selected_refinement_plugin).setOptions(refinementOptions);
        tracking_plugins( selected_tracking_plugin).setOptions(trackingOptions);
        
        % Create the panels
        candidate_plugins( selected_candidate_plugin).createOptionsPanel(h_all.panel_candidate);
        refinement_plugins( selected_refinement_plugin).createOptionsPanel(h_all.panel_fitting);
        tracking_plugins( selected_tracking_plugin).createOptionsPanel(h_all.panel_tracking);
        
        updatePanelPositions();
    end

% Update GUI based on currently set values
% Note: This does not update the plugins!
    function callback_updateMainGUIstate(hObj, event)        
        % Enable/disable output folder
        if get(h_all.cbx_outputFolder, 'Value')
            set(h_all.edit_outputFolder, 'Enable', 'off');
            set(h_all.button_outputFolder, 'Enable', 'off');
        else 
            set(h_all.edit_outputFolder, 'Enable', 'on');
            set(h_all.button_outputFolder, 'Enable', 'on');
        end
        
        % Enable/disable photon conversion
        if get(h_all.cbx_usePhotonConv, 'Value')
            set(h_all.text_bias, 'Enable', 'on');
            set(h_all.text_sensitivity, 'Enable', 'on');
            set(h_all.text_gain, 'Enable', 'on');
            set(h_all.edit_photonBias, 'Enable','on');
            set(h_all.edit_photonSensitivity, 'Enable','on');
            set(h_all.edit_photonGain, 'Enable','on');
        else
            set(h_all.text_bias, 'Enable', 'off');
            set(h_all.text_sensitivity, 'Enable', 'off');
            set(h_all.text_gain, 'Enable', 'off');
            set(h_all.edit_photonBias, 'Enable','off');
            set(h_all.edit_photonSensitivity, 'Enable','off');
            set(h_all.edit_photonGain, 'Enable','off');
        end
              
        
        % Enable/disable tracking panel
        if get(h_all.cbx_enableTracking, 'Value')
            set(h_all.popup_trackingMethod,'Enable','on');
            set(h_all.label_trackingMethod,'Enable','on');
            set(h_all.button_trackingHelp,'Enable','on');
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','on');
        else
            set(h_all.popup_trackingMethod,'Enable','off');
            set(h_all.label_trackingMethod,'Enable','off');
            set(h_all.button_trackingHelp,'Enable','off');
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','off');
        end
        
    end

% Update GUI based on if data loaded from a settings file should be used
    function callback_updateMainGUIstate_loadedData(hObj, event)    
        % Enable/disable panels depending on data loaded from a TNT file
        % Loaded fitting data can only be used if the loaded candidate data
        % is also used.
        if get(h_all.cbx_candidate_loaded, 'Value')
            candidateOptions = GUIinputs.candidateOptions_loaded;
        end
        
        if get(h_all.cbx_refinement_loaded, 'Value')
            refinementOptions = GUIinputs.refinementOptions_loaded;
        end
        
        if get(h_all.cbx_candidate_loaded, 'Value') || get(h_all.cbx_refinement_loaded, 'Value')
            % Update plugins
            selectPluginsBasedOnOptions();
            
            globalOptions.firstFrame = GUIinputs.firstFrame_lastFrame_loaded(1);
            globalOptions.lastFrame = GUIinputs.firstFrame_lastFrame_loaded(2);
            setNum(h_all.edit_firstFrame, globalOptions.firstFrame, true);
            setNum(h_all.edit_lastFrame, globalOptions.lastFrame, true);
            
            set(h_all.text_FrameInterval, 'Enable', 'off');
            set(h_all.text_firstFrame, 'Enable', 'off');
            set(h_all.text_lastFrame, 'Enable', 'off');
            set(h_all.edit_firstFrame, 'Enable','off');
            set(h_all.edit_lastFrame, 'Enable','off');
        else
            set(h_all.text_FrameInterval, 'Enable', 'on');
            set(h_all.text_firstFrame, 'Enable', 'on');
            set(h_all.text_lastFrame, 'Enable', 'on');
            set(h_all.edit_firstFrame, 'Enable','on');
            set(h_all.edit_lastFrame, 'Enable','on');
        end
        
        if get(h_all.cbx_candidate_loaded, 'Value')            
            set(h_all.cbx_refinement_loaded,'Enable','on');
            set(h_all.popup_candidateMethod,'Enable','off');
            set(findall(h_all.panel_candidate,'-property','Enable'), 'Enable','off');
        else
            set(h_all.cbx_refinement_loaded,'Enable','off');
            set(h_all.popup_candidateMethod,'Enable','on');
            set(findall(h_all.panel_candidate,'-property','Enable'), 'Enable','on');
        end
        
        if get(h_all.cbx_refinement_loaded, 'Value')            
            set(h_all.cbx_candidate_loaded,'Enable','off');
            set(h_all.popup_refinementMethod,'Enable','off');
            set(findall(h_all.panel_fitting,'-property','Enable'), 'Enable','off');
        else
            set(h_all.cbx_candidate_loaded,'Enable','on');
            set(h_all.popup_refinementMethod,'Enable','on');
            set(findall(h_all.panel_fitting,'-property','Enable'), 'Enable','on');
        end
        
        
    end

% This function is called if the selection of plugins changes
% It saves the options for the previously selected plugin and builds
% the panel for the newly selected one.
    function callback_updatePlugins(hObj, event)
        selected_candidate_plugin =  get(h_all.popup_candidateMethod,'Value');
        selected_refinement_plugin = get(h_all.popup_refinementMethod,'Value');
        selected_tracking_plugin = get(h_all.popup_trackingMethod,'Value');
        
        % Update panels to display currently selected plugin
        candidate_plugins( selected_candidate_plugin).createOptionsPanel(h_all.panel_candidate);
        refinement_plugins( selected_refinement_plugin).createOptionsPanel(h_all.panel_fitting);
        tracking_plugins( selected_tracking_plugin).createOptionsPanel(h_all.panel_tracking);
        
        % Disable candidate search / refinement panels if the box to use
        % data loaded from a file is checked (the user should not be able
        % to alter the parameters then.
                if get(h_all.cbx_candidate_loaded, 'Value')            
            set(h_all.cbx_refinement_loaded,'Enable','on');
            set(h_all.popup_candidateMethod,'Enable','off');
            set(findall(h_all.panel_candidate,'-property','Enable'), 'Enable','off');
        else
            set(h_all.cbx_refinement_loaded,'Enable','off');
            set(h_all.popup_candidateMethod,'Enable','on');
            set(findall(h_all.panel_candidate,'-property','Enable'), 'Enable','on');
        end
        
        if get(h_all.cbx_refinement_loaded, 'Value')            
            set(h_all.cbx_candidate_loaded,'Enable','off');
            set(h_all.popup_refinementMethod,'Enable','off');
            set(findall(h_all.panel_fitting,'-property','Enable'), 'Enable','off');
        else
            set(h_all.cbx_candidate_loaded,'Enable','on');
            set(h_all.popup_refinementMethod,'Enable','on');
            set(findall(h_all.panel_fitting,'-property','Enable'), 'Enable','on');
        end
        % --------------
        
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
        % -- Spacings --
        units = 'pixels';
        TOPIC_SPACING = 39; % Spacing between topic panels (candidate/tracking etc)
        ABOVE_PANEL_SPACING = 6.5; % Spacing topic panel and the UI elements above the panel (Refinement Method etc.)
        BUTTON_SPACING = 26; % Spacing between last panels and buttons at the bottom
        BOTTOM_SPACING = 9.75; % Spacing between buttons at the bottom and bottom of the GUI window        
        % -- -------- --
        
        set(h_all.panel_candidate,'Units',units);
        panel_candidate_pos = get(h_all.panel_candidate,'Position');
        
        % Relative to candidate detection panel
        set(h_all.label_refinementMethod, 'Units',units);
        label_refinementMethod_pos = get(h_all.label_refinementMethod, 'Position');
        label_refinementMethod_pos(2) = panel_candidate_pos(2)-TOPIC_SPACING-label_refinementMethod_pos(4)/2;
        set(h_all.label_refinementMethod, 'Position',label_refinementMethod_pos);
        
        set(h_all.popup_refinementMethod, 'Units',units);
        pos = get(h_all.popup_refinementMethod, 'Position');
        pos(2) = panel_candidate_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.popup_refinementMethod, 'Position',pos);
        
        set(h_all.button_refinementHelp, 'Units',units);
        pos = get(h_all.button_refinementHelp, 'Position');
        pos(2) = panel_candidate_pos(2)-TOPIC_SPACING-pos(4)/2 + 0.1; % Note the 0.1 is empirical and fits better
        set(h_all.button_refinementHelp, 'Position',pos);
        
        set(h_all.cbx_refinement_loaded, 'Units',units);
        pos = get(h_all.cbx_refinement_loaded, 'Position');
        pos(2) = panel_candidate_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.cbx_refinement_loaded, 'Position',pos);
        
        % Relative to label_fitting_Method
        set(h_all.panel_fitting, 'Units',units);
        panel_fitting_pos = get(h_all.panel_fitting, 'Position');
        panel_fitting_pos(2) = label_refinementMethod_pos(2)-ABOVE_PANEL_SPACING-panel_fitting_pos(4);
        set(h_all.panel_fitting, 'Position',panel_fitting_pos);
        
        % Relative to fitting panel
        set(h_all.label_enableTracking, 'Units',units);
        label_enableTracking_pos = get(h_all.label_enableTracking, 'Position');
        label_enableTracking_pos(2) = panel_fitting_pos(2)-TOPIC_SPACING-label_enableTracking_pos(4)/2;
        set(h_all.label_enableTracking, 'Position',label_enableTracking_pos);
        
        set(h_all.cbx_enableTracking, 'Units',units);
        pos = get(h_all.cbx_enableTracking, 'Position');
        pos(2) = panel_fitting_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.cbx_enableTracking, 'Position',pos);
        
        set(h_all.label_trackingMethod, 'Units',units);
        pos = get(h_all.label_trackingMethod, 'Position');
        pos(2) = panel_fitting_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.label_trackingMethod, 'Position',pos);
        
        set(h_all.popup_trackingMethod, 'Units',units);
        pos = get(h_all.popup_trackingMethod, 'Position');
        pos(2) = panel_fitting_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.popup_trackingMethod, 'Position',pos);
        
        set(h_all.button_trackingHelp, 'Units',units);
        pos = get(h_all.button_trackingHelp, 'Position');
        pos(2) = panel_fitting_pos(2)-TOPIC_SPACING-pos(4)/2 + 0.1; % Note the 0.1 is empirical and fits better
        set(h_all.button_trackingHelp, 'Position',pos);
        
        % Relative to label tracking
        set(h_all.panel_tracking, 'Units',units);
        panel_tracking_pos = get(h_all.panel_tracking, 'Position');
        panel_tracking_pos(2) = label_enableTracking_pos(2)-ABOVE_PANEL_SPACING-panel_tracking_pos(4);
        set(h_all.panel_tracking, 'Position',panel_tracking_pos);
        
        % Relative to tracking panel        
        set(h_all.button_continueForAll, 'Units',units);
        pos = get(h_all.button_continueForAll, 'Position');
        pos(2) = panel_tracking_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.button_continueForAll, 'Position',pos);
        
        set(h_all.button_continue, 'Units',units);
        pos = get(h_all.button_continue, 'Position');
        pos(2) = panel_tracking_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.button_continue, 'Position',pos);
        
        set(h_all.button_preview, 'Units',units);
        pos = get(h_all.button_preview, 'Position');
        pos(2) = panel_tracking_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.button_preview, 'Position',pos);
        
        button_preview_height = pos(4);
        
        
        set(h_all.text_firstFrameTesting, 'Units',units);
        pos = get(h_all.text_firstFrameTesting, 'Position');
        pos(2) = panel_tracking_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.text_firstFrameTesting, 'Position',pos);
        
        set(h_all.edit_firstFrameTesting, 'Units',units);
        pos = get(h_all.edit_firstFrameTesting, 'Position');
        pos(2) = panel_tracking_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.edit_firstFrameTesting, 'Position',pos);
        
        set(h_all.text_lastFrameTesting, 'Units',units);
        pos = get(h_all.text_lastFrameTesting, 'Position');
        pos(2) = panel_tracking_pos(2)-BUTTON_SPACING-button_preview_height;
        set(h_all.text_lastFrameTesting, 'Position',pos);
        
        set(h_all.edit_lastFrameTesting, 'Units',units);
        pos = get(h_all.edit_lastFrameTesting, 'Position');
        pos(2) = panel_tracking_pos(2)-BUTTON_SPACING-button_preview_height;
        set(h_all.edit_lastFrameTesting, 'Position',pos);
        
        % Rescale window based on last element
        set(h_main,'Units',units);
        win_pos = get(h_main,'Position');
        win_top_pos = win_pos(2)+win_pos(4);
        diff_height = pos(2)-BOTTOM_SPACING;
        
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
        win_pos(2) = win_top_pos-win_pos(4); % Keeps the top position constant
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
        
        [~,name,~] = fileparts(outfile); % Append _TNT to filename
        outfile = [path,name,'_TNT.mat'];
        save(outfile,'filename_movie','globalOptions', 'candidateOptions','refinementOptions');
        if(globalOptions.enableTracking) % Save tracking options only if tracking is desired
            save(outfile,'trackingOptions','-append');
        end
    end

% Load settings from a file
    function callback_loadSettings(hObj,event)   
        % Note: we don't load the filename of the movie, as this should be
        % unchanged when loading settings.
        
        [infile, path] = uigetfile({'*.mat','TNT settings'});
        if isfloat(infile);
            return;
        end; % User clicked cancel
        
        % Note: Loading has to be done this way, as variables "can not be
        % added to a static workspace" (e.g. the one of this GUI).
        warning off % We turn warnings off, as trackingOptions might not exist
        allOptions = load([path,infile],'globalOptions', 'candidateOptions','refinementOptions','trackingOptions');
        warning on
        globalOptions   = allOptions.globalOptions;
        candidateOptions = allOptions.candidateOptions;
        refinementOptions   = allOptions.refinementOptions;
        if isfield(allOptions,'trackingOptions')
            trackingOptions  = allOptions.trackingOptions;
        end
        
        setGUIBasedOnOptions();
    end

% Opens a file chooser dialog to select the output folder
    function callback_selectOutputFolder(~, ~)
        outputFolderPath = get(h_all.edit_darkMovie,'String');
        path = [];
        if ~isempty(outputFolderPath)
            [path,~,~] = fileparts(outputFolderPath);
        end
        [path] = uigetdir(path);
        if( isfloat(path)); return; end; % User pressed cancel.
        
        set(h_all.edit_outputFolder,'String',path);
        set(h_all.edit_outputFolder,'Tooltip',path);
    end

% Opens a file chooser dialog to select the dark movie
    function callback_selectDarkMovie(~, ~)
        darkMoviePath = get(h_all.edit_darkMovie,'String');
        path = [];
        if ~isempty(darkMoviePath)
            [path,~,~] = fileparts(darkMoviePath);
        end
        [darkMovie, path] = uigetfile([path,filesep,'*.tif']);
        if( isfloat(darkMovie)); return; end; % User pressed cancel.
        
        set(h_all.edit_darkMovie,'String',[path,darkMovie]);
        set(h_all.edit_darkMovie,'Tooltip',[path,darkMovie]);
    end

% Shows a message dialog with the info of the currently selected plugin
    function callback_helpButtons(hObj, event)
        switch hObj
            case h_all.button_candidateHelp
                name = candidate_plugins(selected_candidate_plugin).name;
                msg = candidate_plugins(selected_candidate_plugin).info;
            case h_all.button_refinementHelp
                name = refinement_plugins(selected_refinement_plugin).name;
                msg = refinement_plugins(selected_refinement_plugin).info;
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
        
        % Accept end as inf, as MATLAB users are used to end as the last element
        if(strcmp(get(hObj,'String'), 'end'))
            value = inf;
        else
            value = round(str2num(get(hObj,'String')));
        end
        
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

% Set all UI fields based on the current value of the options structs
% (globalOptions, candidateOptions, refinementOptions, trackingOptions)
    function setGUIBasedOnOptions()
        % % General Options
        set(h_all.edit_darkMovie,'String', globalOptions.filename_dark_movie);
        set(h_all.edit_darkMovie,'Tooltip',globalOptions.filename_dark_movie);
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

% Set all options structs (globalOptions, candidateOptions, refinementOptions, trackingOptions)
% based on the current state of the uicontrols. Also checks which
% options structs have changed compared to the startup values.
    function storeOptions()
        % % General Options
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
        refinementOptions = refinement_plugins(selected_refinement_plugin).getOptions();
        trackingOptions = tracking_plugins(selected_tracking_plugin).getOptions();
        
        % Store if loaded data should be used
        GUIreturns.use_loaded_candidateData = get(h_all.cbx_candidate_loaded,'Value');
        GUIreturns.use_loaded_refinementData = get(h_all.cbx_refinement_loaded,'Value');
        
        %Check if options were changed compared to intial ones
        GUIreturns.globalOptionsChanged   = ~isequaln(globalOptions_atStartup, globalOptions);
        GUIreturns.candidateOptionsChanged = ~isequaln(candidateOptions_atStartup, candidateOptions);
        GUIreturns.refinementOptionsChanged   = ~isequaln(refinementOptions_atStartup, refinementOptions);
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
        
        % Save output folder
        GUIreturns.outputFolder = get(h_all.edit_outputFolder,'String');
        GUIreturns.outputFolderSameAsMovie = logical(get(h_all.cbx_outputFolder,'Value'));        
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
