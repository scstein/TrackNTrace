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
function [globalOptions,importOptions,candidateOptions,refinementOptions,trackingOptions,postprocOptions, GUIreturns] = settingsGUI(globalOptions,importPlugin,candidateOptions,refinementOptions,trackingOptions, postprocOptions, GUIinputs)
% TrackNTrace main GUI. Here plugins can be selected, their options 
% adjusted. This GUI is used for this management only, while the overall
% program flow is handled by RunTrackNTrace.m
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2016
%
% Extended CT, 2019:
%	- added UI panel for importPlugin
%	- added UI framebinning
%	- added UI timegate and preview
%	- added UI panel for postprocessing
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

% Get importOptions (the import plugin is already loaded, but it's options
% might get changed).
importOptions = importPlugin.getOptions();
% Save options at startup (to check later if options changed)
globalOptions_atStartup = globalOptions;
importOptions_atStartup = importOptions;
candidateOptions_atStartup = candidateOptions;
refinementOptions_atStartup = refinementOptions;
trackingOptions_atStartup = trackingOptions;
postprocOptions_atStartup = postprocOptions;


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
set(h_all.edit_binFrame,'Callback',{@callback_IntEdit,1,inf});
%Photon conversion
set(h_all.cbx_usePhotonConv, 'Callback', @callback_updateMainGUIstate);
set(h_all.edit_photonBias, 'Callback', {@callback_IntEdit,0,inf});
set(h_all.edit_photonSensitivity, 'Callback', {@callback_FloatEdit,1.0,inf});
set(h_all.edit_photonGain, 'Callback', {@callback_IntEdit,1,1000});
%Time gate
if GUIinputs.hasTCSPC
    set(h_all.cbx_useTimegate, 'Enable', 'On');
    set(h_all.cbx_useTimegate, 'Callback', @callback_updateMainGUIstate);
    set(h_all.edit_tgStart, 'Callback', {@callback_IntEdit,0,inf});
    set(h_all.edit_tgEnd, 'Callback', {@callback_IntEdit,0,inf});
    set(h_all.button_showTCSPC, 'Callback', @callback_showTCSPC);
    set(h_all.button_showTCSPC, 'Enable', 'on');
else
    set(h_all.cbx_useTimegate, 'Enable', 'Off');
    set(h_all.cbx_useTimegate, 'Value', false);
    set(h_all.button_showTCSPC, 'Enable', 'Off');
end

% % Import plugin
if size(importPlugin.param_specification,1)>0 %if import plugin has options
    set(h_all.button_importOptions,'Enable','on');
    set(h_all.button_importOptions,'Callback',@callback_importOptionsPanel);
else
    set(h_all.button_importOptions,'Enable','off');
end
% % Candidate plugin
set(h_all.cbx_enableCandidate, 'Callback', @callback_updateMainGUIstate);
set(h_all.popup_candidateMethod, 'Callback', @callback_updatePlugins);
set(h_all.cbx_candidate_loaded, 'Callback', @callback_updateMainGUIstate_loadedData);
if(GUIinputs.show_candidateData_fromFile_cbx)
 set(h_all.cbx_candidate_loaded, 'Value', GUIinputs.use_loaded_candidateData);
else
 set(h_all.cbx_candidate_loaded, 'Visible', 'off');
end

% % Refinement plugin
set(h_all.cbx_enableRefinement, 'Callback', @callback_updateMainGUIstate);
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
set(h_all.cbx_tracking_loaded, 'Callback', @callback_updateMainGUIstate_loadedData);
if(GUIinputs.show_trackingData_fromFile_cbx)
 set(h_all.cbx_tracking_loaded, 'Value', GUIinputs.use_loaded_trackingData);
else
 set(h_all.cbx_tracking_loaded, 'Visible', 'off');
end

% % Postproc plugin
set(h_all.cbx_enablePostproc, 'Callback', @callback_updateMainGUIstate);
set(h_all.popup_postprocMethod, 'Callback', @callback_updatePlugins);
set(h_all.cbx_postproc_loaded, 'Callback', @callback_updateMainGUIstate_loadedData);
if(GUIinputs.show_postprocData_fromFile_cbx)
 set(h_all.cbx_postproc_loaded, 'Value', GUIinputs.use_loaded_postprocData);
else
 set(h_all.cbx_postproc_loaded, 'Visible', 'off');
end

% % Help buttons
set(h_all.button_candidateHelp, 'Callback', @callback_helpButtons);
set(h_all.button_refinementHelp, 'Callback', @callback_helpButtons);
set(h_all.button_trackingHelp, 'Callback', @callback_helpButtons);
set(h_all.button_postprocHelp, 'Callback', @callback_helpButtons);

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
postproc_plugins = [];

% Save last selected plugin per category
selected_candidate_plugin = -1;
selected_refinement_plugin = -1;
selected_tracking_plugin = -1;
selected_postproc_plugin = -1;

% Load the plugins
if GUIinputs.showStartupInformation % Show warnings only on startup
    fprintf('TNT: Loading plugins ...\n')
    loadPlugins();
    fprintf('TNT: Successfully imported %i plugins (%i candidate detection, %i refinement, %i tracking, %i postprocessing).\n',numel(candidate_plugins)+numel(refinement_plugins)+numel(tracking_plugins),numel(candidate_plugins),numel(refinement_plugins),numel(tracking_plugins),numel(postproc_plugins));
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
                    case 4 % Postprocessing
                        postproc_plugins = [postproc_plugins, plugin];
                    case {5,6} % Postprocessing, Import, Export
                        % Ignore
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
        found_postproc_plugin = numel(postproc_plugins)>0;
        
        if not(found_candidate_plugin)
            error('No candidate detection plugin detected.');
        end
        if not(found_refinement_plugin)
            error('No refinement plugin detected.');
        end
        if not(found_tracking_plugin)
            warning('No tracking plugin detected.')
        end
        if not(found_postproc_plugin)
            warning('No tracking plugin detected.')
        end
        
        % Set popup choices
        set(h_all.popup_candidateMethod, 'String', {candidate_plugins(:).name});
        set(h_all.popup_refinementMethod, 'String', {refinement_plugins(:).name});
        if found_tracking_plugin
            set(h_all.popup_trackingMethod, 'String', {tracking_plugins(:).name});
        end
        if found_postproc_plugin
            set(h_all.popup_postprocMethod, 'String', {postproc_plugins(:).name});
        end
    end

% Select plugins based on the current candidateOptions, refinementOptions,
% trackingOptions. After selection, the plugin panels are constructed
    function selectPluginsBasedOnOptions(loadCandidate,loadRefinement,loadTracking,loadPostproc)
        loadCandidate = ~exist('loadCandidate','var') || loadCandidate;
        loadRefinement = ~exist('loadRefinement','var') || loadRefinement;
        loadTracking = ~exist('loadPostproc','var') || loadTracking;
        loadPostproc = ~exist('loadPostproc','var') || loadPostproc;
        found_tracking_plugin = numel(tracking_plugins)>0;
        found_postproc_plugin = numel(postproc_plugins)>0;
        
        % Candidate detection
        if loadCandidate
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
                default_plugin_index = strcmp({candidate_plugins(:).name},GUIinputs.TNToptions.defaultCandidatePlugin);
                default_plugin_index = find(default_plugin_index);
                if isempty(default_plugin_index)
                    warning off backtrace
                    warning('Default candidate plugin ''%s'' not found. Selecting first one.',GUIinputs.TNToptions.defaultCandidatePlugin);
                    warning on backtrace
                    selected_candidate_plugin = 1;
                else
                    selected_candidate_plugin = default_plugin_index;
                end
            end
        end
        
        % Refinement
        if loadRefinement
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
                default_plugin_index = strcmp({refinement_plugins(:).name},GUIinputs.TNToptions.defaultRefinementPlugin);
                default_plugin_index = find(default_plugin_index);
                if isempty(default_plugin_index)
                    warning off backtrace
                    warning('Default fitting plugin ''%s'' not found. Selecting first one.',GUIinputs.TNToptions.defaultRefinementPlugin);
                    warning on backtrace
                    selected_refinement_plugin = 1;
                else
                    selected_refinement_plugin = default_plugin_index;
                end
            end
        end
        
        % Tracking
        if found_tracking_plugin && loadTracking
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
                default_plugin_index = strcmp({tracking_plugins(:).name},GUIinputs.TNToptions.defaultTrackingPlugin);
                default_plugin_index = find(default_plugin_index);
                if isempty(default_plugin_index)
                    warning off backtrace
                    warning('Default tracking plugin ''%s'' not found. Selecting first one.',GUIinputs.TNToptions.defaultTrackingPlugin);
                    warning on backtrace
                    selected_tracking_plugin = 1;
                else
                    selected_tracking_plugin = default_plugin_index;
                end
            end
        end
        
        % Postprocessing
        if found_postproc_plugin && loadPostproc
            if ~isempty(postprocOptions)
                selected_postproc_plugin = -1;
                % Search for a loaded plugin with the same name
                for iPlug = 1:numel(postproc_plugins)
                    if strcmp(postprocOptions.plugin_name,postproc_plugins(iPlug).name)
                        selected_postproc_plugin = iPlug;
                    end
                end
                
                if selected_postproc_plugin < 0
                    error('Plugin ''%s'' not found.',postprocOptions.plugin_name) ;
                end
            else
                % Select the default plugin if it is found
                default_plugin_index = strcmp({postproc_plugins(:).name},GUIinputs.TNToptions.defaultPostprocPlugin);
                default_plugin_index = find(default_plugin_index);
                if isempty(default_plugin_index)
                    warning off backtrace
                    warning('Default postproc plugin ''%s'' not found. Selecting first one.',GUIinputs.TNToptions.defaultPostprocPlugin);
                    warning on backtrace
                    selected_postproc_plugin = 1;
                else
                    selected_postproc_plugin = default_plugin_index;
                end
            end
        end
        % Set popups to correct plugin name
        set(h_all.popup_candidateMethod, 'Value', selected_candidate_plugin);
        set(h_all.popup_refinementMethod, 'Value', selected_refinement_plugin);
        if found_tracking_plugin
            set(h_all.popup_trackingMethod, 'Value', selected_tracking_plugin);
        end
        if found_postproc_plugin
            set(h_all.popup_postprocMethod, 'Value', selected_postproc_plugin);
        end
        
        % Build panels by invoking the selected plugins function
        if loadCandidate
            candidate_plugins( selected_candidate_plugin).setOptions(candidateOptions);
        end
        if loadRefinement
            refinement_plugins( selected_refinement_plugin).setOptions(refinementOptions);
        end
        if found_tracking_plugin && loadTracking
            tracking_plugins( selected_tracking_plugin).setOptions(trackingOptions);
        end
        if found_postproc_plugin && loadPostproc
            postproc_plugins( selected_postproc_plugin).setOptions(postprocOptions);
        end
        
        % Create the panels
        candidate_plugins( selected_candidate_plugin).createOptionsPanel(h_all.panel_candidate);
        refinement_plugins( selected_refinement_plugin).createOptionsPanel(h_all.panel_refinement);
        if found_tracking_plugin
            tracking_plugins( selected_tracking_plugin).createOptionsPanel(h_all.panel_tracking);
        end
        if found_postproc_plugin
            postproc_plugins( selected_postproc_plugin).createOptionsPanel(h_all.panel_postproc);
        end
        
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
        if get(h_all.cbx_usePhotonConv, 'Value') && strcmpi(get(h_all.cbx_usePhotonConv,'Enable'),'on')
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
        
        % Enable/disable timegate
        if get(h_all.cbx_useTimegate, 'Value')
            set(h_all.text_tgStart, 'Enable', 'on');
            set(h_all.text_tgEnd, 'Enable', 'on');
            set(h_all.edit_tgStart, 'Enable', 'on');
            set(h_all.edit_tgEnd, 'Enable', 'on');
        else
            set(h_all.text_tgStart, 'Enable', 'off');
            set(h_all.text_tgEnd, 'Enable', 'off');
            set(h_all.edit_tgStart, 'Enable', 'off');
            set(h_all.edit_tgEnd, 'Enable', 'off');
        end
        
        updatePanelEnable();
        
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
        
        if get(h_all.cbx_tracking_loaded, 'Value')
            trackingOptions = GUIinputs.trackingOptions_loaded;
        end
        
        if get(h_all.cbx_postproc_loaded, 'Value')
            postprocOptions = GUIinputs.postprocOptions_loaded;
        end
        
        if get(h_all.cbx_candidate_loaded, 'Value') && get(h_all.cbx_enableCandidate, 'Value')
            % Update plugins
            selectPluginsBasedOnOptions(...
                get(h_all.cbx_candidate_loaded, 'Value'), ...
                get(h_all.cbx_refinement_loaded, 'Value'),...
                get(h_all.cbx_tracking_loaded, 'Value'),...
                get(h_all.cbx_postproc_loaded, 'Value'));
            
            % Set and disable global options influencing the analysis
            set(h_all.edit_darkMovie,'String',globalOptions_atStartup.filename_dark_movie);
            
            setNum(h_all.edit_firstFrame, globalOptions_atStartup.firstFrame, true);
            setNum(h_all.edit_lastFrame, globalOptions_atStartup.lastFrame, true);
            setNum(h_all.edit_binFrame, globalOptions_atStartup.binFrame, true);

            set(h_all.cbx_usePhotonConv,'Value',globalOptions_atStartup.usePhotonConversion)
            setNum(h_all.edit_photonBias, globalOptions_atStartup.photonBias, false);
            setNum(h_all.edit_photonSensitivity, globalOptions_atStartup.photonSensitivity, false);
            setNum(h_all.edit_photonGain, globalOptions_atStartup.photonGain, false);
            
            set(h_all.cbx_useTimegate,'Value',globalOptions_atStartup.useTimegate)
            
            setNum(h_all.edit_tgStart, globalOptions_atStartup.tgStart, true);
            setNum(h_all.edit_tgEnd, globalOptions_atStartup.tgEnd, true);
                                
            set(h_all.text_FrameInterval, 'Enable', 'off');
            set(h_all.text_firstFrame, 'Enable', 'off');
            set(h_all.text_lastFrame, 'Enable', 'off');
            set(h_all.text_binFrame, 'Enable', 'off');
            set(h_all.edit_firstFrame, 'Enable','off');
            set(h_all.edit_lastFrame, 'Enable','off');
            set(h_all.edit_binFrame, 'Enable','off');
            
            
            set(h_all.text_darkMovie, 'Enable', 'off');
            set(h_all.edit_darkMovie, 'Enable', 'off');
            set(h_all.button_darkMovie, 'Enable', 'off');
            set(h_all.label_usePhotonConv, 'Enable','off');
            set(h_all.cbx_usePhotonConv, 'Enable','off');
            set(h_all.label_useTimegate, 'Enable','off');
            set(h_all.cbx_useTimegate, 'Enable','off');
        else
            set(h_all.text_FrameInterval, 'Enable', 'on');
            set(h_all.text_firstFrame, 'Enable', 'on');
            set(h_all.text_lastFrame, 'Enable', 'on');
            set(h_all.text_binFrame, 'Enable', 'on');
            set(h_all.edit_firstFrame, 'Enable','on');
            set(h_all.edit_lastFrame, 'Enable','on');
            set(h_all.edit_binFrame, 'Enable','on');
                        
            set(h_all.text_darkMovie, 'Enable', 'on');
            set(h_all.edit_darkMovie, 'Enable', 'on');
            set(h_all.button_darkMovie, 'Enable', 'on');
            set(h_all.label_usePhotonConv, 'Enable','on');
            set(h_all.cbx_usePhotonConv, 'Enable','on');
            
            set(h_all.label_useTimegate, 'Enable','on');
            if GUIinputs.hasTCSPC
                set(h_all.cbx_useTimegate, 'Enable', 'on');
            else
                set(h_all.cbx_useTimegate, 'Enable', 'off');
                set(h_all.cbx_useTimegate, 'Value', false);
            end
        end
        
        callback_updateMainGUIstate(); 
    end

% This function is called if the selection of plugins changes
% It saves the options for the previously selected plugin and builds
% the panel for the newly selected one.
    function callback_updatePlugins(hObj, event)
        selected_candidate_plugin =  get(h_all.popup_candidateMethod,'Value');
        selected_refinement_plugin = get(h_all.popup_refinementMethod,'Value');
        selected_tracking_plugin = get(h_all.popup_trackingMethod,'Value');
        selected_postproc_plugin = get(h_all.popup_postprocMethod,'Value');
        
        % Update panels to display currently selected plugin
        candidate_plugins( selected_candidate_plugin).createOptionsPanel(h_all.panel_candidate);
        refinement_plugins( selected_refinement_plugin).createOptionsPanel(h_all.panel_refinement);
        tracking_plugins( selected_tracking_plugin).createOptionsPanel(h_all.panel_tracking);
        postproc_plugins( selected_postproc_plugin).createOptionsPanel(h_all.panel_postproc);
               
        updatePanelPositions();
        updatePanelEnable();
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
        
        set(h_all.cbx_enableRefinement, 'Units',units);
        pos = get(h_all.cbx_enableRefinement, 'Position');
        pos(2) = panel_candidate_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.cbx_enableRefinement, 'Position',pos);
        
        set(h_all.label_enableRefinement, 'Units',units);
        pos = get(h_all.label_enableRefinement, 'Position');
        pos(2) = panel_candidate_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.label_enableRefinement, 'Position',pos);
        
        % Relative to label_fitting_Method
        set(h_all.panel_refinement, 'Units',units);
        panel_refinement_pos = get(h_all.panel_refinement, 'Position');
        panel_refinement_pos(2) = label_refinementMethod_pos(2)-ABOVE_PANEL_SPACING-panel_refinement_pos(4);
        set(h_all.panel_refinement, 'Position',panel_refinement_pos);
        
        % Relative to fitting panel
        set(h_all.label_enableTracking, 'Units',units);
        label_enableTracking_pos = get(h_all.label_enableTracking, 'Position');
        label_enableTracking_pos(2) = panel_refinement_pos(2)-TOPIC_SPACING-label_enableTracking_pos(4)/2;
        set(h_all.label_enableTracking, 'Position',label_enableTracking_pos);
        
        set(h_all.cbx_enableTracking, 'Units',units);
        pos = get(h_all.cbx_enableTracking, 'Position');
        pos(2) = panel_refinement_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.cbx_enableTracking, 'Position',pos);
        
        set(h_all.cbx_tracking_loaded, 'Units',units);
        pos = get(h_all.cbx_tracking_loaded, 'Position');
        pos(2) = panel_refinement_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.cbx_tracking_loaded, 'Position',pos);
        
        set(h_all.label_trackingMethod, 'Units',units);
        pos = get(h_all.label_trackingMethod, 'Position');
        pos(2) = panel_refinement_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.label_trackingMethod, 'Position',pos);
        
        set(h_all.popup_trackingMethod, 'Units',units);
        pos = get(h_all.popup_trackingMethod, 'Position');
        pos(2) = panel_refinement_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.popup_trackingMethod, 'Position',pos);
        
        set(h_all.button_trackingHelp, 'Units',units);
        pos = get(h_all.button_trackingHelp, 'Position');
        pos(2) = panel_refinement_pos(2)-TOPIC_SPACING-pos(4)/2 + 0.1; % Note the 0.1 is empirical and fits better
        set(h_all.button_trackingHelp, 'Position',pos);
        
        % Relative to label tracking
        set(h_all.panel_tracking, 'Units',units);
        panel_tracking_pos = get(h_all.panel_tracking, 'Position');
        panel_tracking_pos(2) = label_enableTracking_pos(2)-ABOVE_PANEL_SPACING-panel_tracking_pos(4);
        set(h_all.panel_tracking, 'Position',panel_tracking_pos);
        
        % Relative to tracking panel
        set(h_all.label_enablePostproc, 'Units',units);
        label_enablePostproc_pos = get(h_all.label_enablePostproc, 'Position');
        label_enablePostproc_pos(2) = panel_tracking_pos(2)-TOPIC_SPACING-label_enablePostproc_pos(4)/2;
        set(h_all.label_enablePostproc, 'Position',label_enablePostproc_pos);
        
        set(h_all.cbx_enablePostproc, 'Units',units);
        pos = get(h_all.cbx_enablePostproc, 'Position');
        pos(2) = panel_tracking_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.cbx_enablePostproc, 'Position',pos);
        
        set(h_all.cbx_postproc_loaded, 'Units',units);
        pos = get(h_all.cbx_postproc_loaded, 'Position');
        pos(2) = panel_tracking_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.cbx_postproc_loaded, 'Position',pos);
        
        set(h_all.label_postprocMethod, 'Units',units);
        pos = get(h_all.label_postprocMethod, 'Position');
        pos(2) = panel_tracking_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.label_postprocMethod, 'Position',pos);
        
        set(h_all.popup_postprocMethod, 'Units',units);
        pos = get(h_all.popup_postprocMethod, 'Position');
        pos(2) = panel_tracking_pos(2)-TOPIC_SPACING-pos(4)/2;
        set(h_all.popup_postprocMethod, 'Position',pos);
        
        set(h_all.button_postprocHelp, 'Units',units);
        pos = get(h_all.button_postprocHelp, 'Position');
        pos(2) = panel_tracking_pos(2)-TOPIC_SPACING-pos(4)/2 + 0.1; % Note the 0.1 is empirical and fits better
        set(h_all.button_postprocHelp, 'Position',pos);
        
        % Relative to label postproc
        set(h_all.panel_postproc, 'Units',units);
        panel_postproc_pos = get(h_all.panel_postproc, 'Position');
        panel_postproc_pos(2) = label_enablePostproc_pos(2)-ABOVE_PANEL_SPACING-panel_postproc_pos(4);
        set(h_all.panel_postproc, 'Position',panel_postproc_pos);
        
        % Relative to postproc panel        
        set(h_all.button_continueForAll, 'Units',units);
        pos = get(h_all.button_continueForAll, 'Position');
        pos(2) = panel_postproc_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.button_continueForAll, 'Position',pos);
        
        set(h_all.button_continue, 'Units',units);
        pos = get(h_all.button_continue, 'Position');
        pos(2) = panel_postproc_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.button_continue, 'Position',pos);
        
        set(h_all.button_preview, 'Units',units);
        pos = get(h_all.button_preview, 'Position');
        pos(2) = panel_postproc_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.button_preview, 'Position',pos);
        
        button_preview_height = pos(4);
        
        
        set(h_all.text_firstFrameTesting, 'Units',units);
        pos = get(h_all.text_firstFrameTesting, 'Position');
        pos(2) = panel_postproc_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.text_firstFrameTesting, 'Position',pos);
        
        set(h_all.edit_firstFrameTesting, 'Units',units);
        pos = get(h_all.edit_firstFrameTesting, 'Position');
        pos(2) = panel_postproc_pos(2)-BUTTON_SPACING-pos(4);
        set(h_all.edit_firstFrameTesting, 'Position',pos);
        
        set(h_all.text_lastFrameTesting, 'Units',units);
        pos = get(h_all.text_lastFrameTesting, 'Position');
        pos(2) = panel_postproc_pos(2)-BUTTON_SPACING-button_preview_height;
        set(h_all.text_lastFrameTesting, 'Position',pos);
        
        set(h_all.edit_lastFrameTesting, 'Units',units);
        pos = get(h_all.edit_lastFrameTesting, 'Position');
        pos(2) = panel_postproc_pos(2)-BUTTON_SPACING-button_preview_height;
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

% This function enables/disables the plugin GUI based on the state of the
% enable and 'use loaded data' checkbox.
    function updatePanelEnable()
        % If loaded data is used the popup and panel are disabled and the
        % loaded checkbox of the next panel is activated
                      
        % Enable/disable candidate panel
        enableState = 'off';
        if get(h_all.cbx_enableCandidate, 'Value')
            if ~get(h_all.cbx_candidate_loaded, 'Value')
                enableState = 'on';
                set(h_all.cbx_refinement_loaded, 'Value',false); % Unchecking use-loaded also unchecks the use loaded of the next step.
            end
            set(h_all.cbx_candidate_loaded,'Enable','on');
            set(h_all.button_candidateHelp,'Enable','on');
            set(h_all.cbx_enableRefinement,'Enable','on');
        else
            set(h_all.cbx_candidate_loaded,'Enable','off');
            set(h_all.button_candidateHelp,'Enable','off');
            set(h_all.cbx_enableRefinement,'Enable','off');
        end
        set(h_all.popup_candidateMethod,'Enable',enableState);
        set(h_all.label_candidateMethod,'Enable',enableState);
        set(findall(h_all.panel_candidate,'-property','Enable'), 'Enable',enableState);
              
        % Enable/disable refinement panel
        enableState = 'off';
        if get(h_all.cbx_enableRefinement, 'Value') && get(h_all.cbx_enableCandidate, 'Value') 
            if ~get(h_all.cbx_refinement_loaded, 'Value')
                enableState = 'on';
                set(h_all.cbx_tracking_loaded, 'Value',false); % Unchecking use-loaded also unchecks the use loaded of the next step.
            else
                set(h_all.cbx_candidate_loaded,'Enable','off'); % The checkbox use-loaded of the previous step is disabled if it is checked for this step.
            end
            if get(h_all.cbx_candidate_loaded, 'Value')
                set(h_all.cbx_refinement_loaded,'Enable','on');
            else
                set(h_all.cbx_refinement_loaded,'Enable','off');
            end
            set(h_all.button_refinementHelp,'Enable','on');
            set(h_all.cbx_enableTracking, 'Enable','on');
        else
            set(h_all.cbx_refinement_loaded,'Enable','off');
            set(h_all.button_refinementHelp,'Enable','off');
            set(h_all.cbx_enableTracking, 'Enable','off');
        end
        set(h_all.popup_refinementMethod,'Enable',enableState);
        set(h_all.label_refinementMethod,'Enable',enableState);
        set(findall(h_all.panel_refinement,'-property','Enable'), 'Enable',enableState);
        
        % Enable/disable tracking panel
        enableState = 'off';
        if get(h_all.cbx_enableTracking, 'Value') && get(h_all.cbx_enableRefinement, 'Value') && get(h_all.cbx_enableCandidate, 'Value')
            if ~get(h_all.cbx_tracking_loaded, 'Value')
                enableState = 'on';
                set(h_all.cbx_postproc_loaded, 'Value',false); % Unchecking use-loaded also unchecks the use loaded of the next step.
            else
                set(h_all.cbx_refinement_loaded,'Enable','off'); % The checkbox use-loaded of the previous step is disabled if it is checked for this step.
            end
            if get(h_all.cbx_refinement_loaded, 'Value')
                set(h_all.cbx_tracking_loaded,'Enable','on');
            else
                set(h_all.cbx_tracking_loaded,'Enable','off');
            end
            set(h_all.button_trackingHelp,'Enable','on');
            set(h_all.cbx_enablePostproc, 'Enable','on');
        else
            set(h_all.cbx_tracking_loaded,'Enable','off');
            set(h_all.button_trackingHelp,'Enable','off');
            set(h_all.cbx_enablePostproc, 'Enable','off');
        end
        set(h_all.popup_trackingMethod,'Enable',enableState);
        set(h_all.label_trackingMethod,'Enable',enableState);
        set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable',enableState);
        
        % Enable/disable postproc panel
        enableState = 'off';
        if get(h_all.cbx_enablePostproc, 'Value') && get(h_all.cbx_enableTracking, 'Value') && get(h_all.cbx_enableRefinement, 'Value') && get(h_all.cbx_enableCandidate, 'Value')
            if ~get(h_all.cbx_postproc_loaded, 'Value')
                enableState = 'on';
            else
                set(h_all.cbx_tracking_loaded,'Enable','off'); % The checkbox use-loaded of the previous step is disabled if it is checked for this step.
            end
            if get(h_all.cbx_tracking_loaded, 'Value')
                set(h_all.cbx_postproc_loaded,'Enable','on');
            else
                set(h_all.cbx_postproc_loaded,'Enable','off');
            end
            set(h_all.button_postprocHelp,'Enable','on');
        else
            set(h_all.cbx_postproc_loaded,'Enable','off');
            set(h_all.button_postprocHelp,'Enable','off');
        end
        set(h_all.popup_postprocMethod,'Enable',enableState);
        set(h_all.label_postprocMethod,'Enable',enableState);
        set(findall(h_all.panel_postproc,'-property','Enable'), 'Enable',enableState);
              
%         % Enable/disable tracking panel
%         if get(h_all.cbx_enableTracking, 'Value') && get(h_all.cbx_enableRefinement, 'Value') && get(h_all.cbx_enableCandidate, 'Value')
%             set(h_all.popup_trackingMethod,'Enable','on');
%             set(h_all.label_trackingMethod,'Enable','on');
%             set(h_all.button_trackingHelp,'Enable','on');
%             set(h_all.cbx_tracking_loaded,'Enable','on');
%             set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','on');
%             set(h_all.cbx_enablePostproc, 'Enable','on');
%         else
%             set(h_all.popup_trackingMethod,'Enable','off');
%             set(h_all.label_trackingMethod,'Enable','off');
%             set(h_all.button_trackingHelp,'Enable','off');
%             set(h_all.cbx_tracking_loaded,'Enable','off');
%             set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','off');
%             set(h_all.cbx_enablePostproc, 'Enable','off');
%         end
%               
%         % Enable/disable postproc panel
%         if get(h_all.cbx_enablePostproc, 'Value') && get(h_all.cbx_enableTracking, 'Value') && get(h_all.cbx_enableRefinement, 'Value') && get(h_all.cbx_enableCandidate, 'Value')
%             set(h_all.popup_postprocMethod,'Enable','on');
%             set(h_all.label_postprocMethod,'Enable','on');
%             set(h_all.button_postprocHelp,'Enable','on');
%             set(h_all.cbx_postproc_loaded,'Enable','on');
%             set(findall(h_all.panel_postproc,'-property','Enable'), 'Enable','on');
%         else
%             set(h_all.popup_postprocMethod,'Enable','off');
%             set(h_all.label_postprocMethod,'Enable','off');
%             set(h_all.button_postprocHelp,'Enable','off');
%             set(h_all.cbx_postproc_loaded,'Enable','off');
%             set(findall(h_all.panel_postproc,'-property','Enable'), 'Enable','off');
%         end
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
        
        % Save the globalOptions and Options of the enabled steps
        savevars = {'filename_movie','globalOptions','importOptions'};
        if(globalOptions.enableCandidate)
            savevars = [savevars,'candidateOptions'];
        end
        if(globalOptions.enableRefinement)
            savevars = [savevars,'refinementOptions'];
        end
        if(globalOptions.enableTracking)
            savevars = [savevars,'trackingOptions'];
        end
        if(globalOptions.enablePostproc)
            savevars = [savevars,'postprocOptions'];
        end
        
        save(outfile,savevars{:});
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
        allOptions = load([path,infile],'globalOptions', 'importOptions', 'candidateOptions','refinementOptions','trackingOptions','postprocOptions');
        warning on
        
        % replacefields avoids errors when options are added or their names are
        % changed and prevents data stored in the option struct to be transfered
        % to the new file.
        globalOptions = replacefields(globalOptions,allOptions.globalOptions);
        if ~GUIinputs.hasTCSPC
            globalOptions.useTimegate = false;
        end
        if isfield(allOptions,'importOptions') && strcmp(importOptions.plugin_name,allOptions.importOptions.plugin_name)
            % Only set importOptions if the import plugin matches.
            importPlugin.setOptions(allOptions.importOptions);
            importOptions = importPlugin.getOptions();
        end
        
        if isfield(allOptions,'candidateOptions')
            candidateOptions = candidate_plugins(selected_candidate_plugin).getOptions();
            if isfield(allOptions.candidateOptions,'plugin_name') && ~strcmpi(candidateOptions.plugin_name,allOptions.candidateOptions.plugin_name)
                % update plugin type to get the available options
                candidateOptions.plugin_name = allOptions.candidateOptions.plugin_name;
                selectPluginsBasedOnOptions(true,false,false,false);
                candidateOptions = candidate_plugins(selected_candidate_plugin).getOptions();
            end
            candidateOptions = replacefields(candidateOptions,allOptions.candidateOptions);
        end
        if isfield(allOptions,'refinementOptions')
            refinementOptions = refinement_plugins(selected_refinement_plugin).getOptions();
            if isfield(allOptions.refinementOptions,'plugin_name') && ~strcmpi(refinementOptions.plugin_name,allOptions.refinementOptions.plugin_name)
                % update plugin type to get the available options
                refinementOptions.plugin_name = allOptions.refinementOptions.plugin_name;
                selectPluginsBasedOnOptions(false,true,false,false);
                refinementOptions = refinement_plugins(selected_refinement_plugin).getOptions();
            end
            refinementOptions = replacefields(refinementOptions,allOptions.refinementOptions);
        end
        if isfield(allOptions,'trackingOptions')
            trackingOptions = tracking_plugins(selected_tracking_plugin).getOptions();
            if isfield(allOptions.trackingOptions,'plugin_name') && ~strcmpi(trackingOptions.plugin_name,allOptions.trackingOptions.plugin_name)
                % update plugin type to get the available options
                trackingOptions.plugin_name = allOptions.trackingOptions.plugin_name;
                selectPluginsBasedOnOptions(false,false,true,false);
                trackingOptions = tracking_plugins(selected_tracking_plugin).getOptions();
            end
            trackingOptions = replacefields(trackingOptions,allOptions.trackingOptions);
        end
        if isfield(allOptions,'postprocOptions')
            postprocOptions = postproc_plugins(selected_postproc_plugin).getOptions();
            if isfield(allOptions.postprocOptions,'plugin_name') && ~strcmpi(postprocOptions.plugin_name,allOptions.postprocOptions.plugin_name)
                % update plugin type to get the available options
                postprocOptions.plugin_name = allOptions.postprocOptions.plugin_name;
                selectPluginsBasedOnOptions(false,false,false,true);
                postprocOptions = postproc_plugins(selected_postproc_plugin).getOptions();
            end
            postprocOptions = replacefields(postprocOptions,allOptions.postprocOptions);
        end
        
        % Disable use of loaded data
        set(h_all.cbx_candidate_loaded, 'Value', false);
        set(h_all.cbx_refinement_loaded, 'Value', false);
        set(h_all.cbx_tracking_loaded, 'Value', false);
        set(h_all.cbx_postproc_loaded, 'Value', false);        
        
        setGUIBasedOnOptions();
        callback_updateMainGUIstate_loadedData();
    end

    function old = replacefields(old,new)
        % replaces fieldvalues in struct OLD with values of struct NEW.
        % No new fields are added. Fields not present in NEW are not changed.
        fns = fieldnames(old);
        for fn = fns(:)'
            if isfield(new,fn{1})
                old.(fn{1}) = new.(fn{1});
            end
        end        
    end

% Opens a file chooser dialog to select the output folder
    function callback_selectOutputFolder(~, ~)
        outputFolderPath = get(h_all.edit_outputFolder,'String');
        path = pwd;
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
        path = pwd;
        if ~isempty(darkMoviePath)
            [path,~,~] = fileparts(darkMoviePath);
        end
        [darkMovie, path] = uigetfile([path,filesep,'*.tif']);
        if( isfloat(darkMovie)); return; end; % User pressed cancel.
        
        set(h_all.edit_darkMovie,'String',[path,darkMovie]);
        set(h_all.edit_darkMovie,'Tooltip',[path,darkMovie]);
    end

% Show the time gates and the TCSPC histogram
    function callback_showTCSPC(~, ~)
        if isfield(importOptions.info,'getTCSPC') && isa(importOptions.info.getTCSPC,'function_handle')
            figure;
            showTimegate(filename_movie,importOptions,[str2double(get(h_all.edit_tgStart,'String')),str2double(get(h_all.edit_tgEnd,'String'))]);
        else
            fprintf('TNT: The import plugin %s does not support preview of the TCSPC.\n', importOptions.plugin_name)
        end
    end

% Create importOptions modal dialog
    function callback_importOptionsPanel(~, ~)
        pos = get(h_main,'Position');
        pos = [pos(1)+10 pos(2)+pos(4)-200 pos(3)-20 100];
        fig = figure('WindowStyle','modal','Name','Import Options', 'NumberTitle','off','Resize','off','Visible','off','Position',pos);
        imPanel = uipanel('Parent',fig,'Units','normalized','Position',[0.05 0.05 0.9 0.9]);
        set(imPanel,'Units','points');
        initPos = get(imPanel,'Position');
        importPlugin.createOptionsPanel(imPanel);
        newPos = get(imPanel,'Position');
        set(fig,'Position',get(fig,'Position')+[0 0 0 newPos(4) - initPos(4)]);
        set(imPanel,'Position',[initPos(1:2) newPos(3:4)]);
        % Diable options if loaded candidates are used (view only)
        if get(h_all.cbx_candidate_loaded, 'Value')
            set(findall(imPanel,'-property','Enable'), 'Enable','off');
        end
        set(fig,'Visible','on');
        uiwait(fig);
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
            case h_all.button_postprocHelp
                name = postproc_plugins(selected_postproc_plugin).name;
                msg = postproc_plugins(selected_postproc_plugin).info;
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
        setNum(h_all.edit_binFrame, globalOptions.binFrame, true);
        setNum(h_all.edit_firstFrameTesting, globalOptions.firstFrameTesting, true);
        setNum(h_all.edit_lastFrameTesting, globalOptions.lastFrameTesting, true);
        
        % Photon conversion
        set(h_all.cbx_usePhotonConv,'Value',globalOptions.usePhotonConversion)
        setNum(h_all.edit_photonBias,globalOptions.photonBias, true);
        setNum(h_all.edit_photonSensitivity,globalOptions.photonSensitivity);
        setNum(h_all.edit_photonGain,globalOptions.photonGain, true);
        
        % Timegate
        set(h_all.cbx_useTimegate,'Value',globalOptions.useTimegate)
        setNum(h_all.edit_tgStart,globalOptions.tgStart, true);
        setNum(h_all.edit_tgEnd,globalOptions.tgEnd, true);
        
        % % enable/disable Plugins
        set(h_all.cbx_enableCandidate,'Value', globalOptions.enableCandidate);
        set(h_all.cbx_enableRefinement,'Value', globalOptions.enableRefinement);
        set(h_all.cbx_enableTracking,'Value', globalOptions.enableTracking);
        set(h_all.cbx_enablePostproc,'Value', globalOptions.enablePostproc);
        
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
        globalOptions.binFrame =  getNum(h_all.edit_binFrame);
        globalOptions.firstFrameTesting = getNum(h_all.edit_firstFrameTesting);
        globalOptions.lastFrameTesting = getNum(h_all.edit_lastFrameTesting);
        
        globalOptions.usePhotonConversion = logical(get(h_all.cbx_usePhotonConv,'Value'));
        globalOptions.photonBias = getNum(h_all.edit_photonBias);
        globalOptions.photonSensitivity = getNum(h_all.edit_photonSensitivity);
        globalOptions.photonGain = getNum(h_all.edit_photonGain);
        
        globalOptions.useTimegate = logical(get(h_all.cbx_useTimegate,'Value'));
        globalOptions.tgStart = getNum(h_all.edit_tgStart);
        globalOptions.tgEnd = getNum(h_all.edit_tgEnd);
        
        % enable/disable Plugins
        globalOptions.enableCandidate = logical(get(h_all.cbx_enableCandidate,'Value')); 
        globalOptions.enableRefinement = logical(get(h_all.cbx_enableRefinement,'Value')); 
        globalOptions.enableTracking = logical(get(h_all.cbx_enableTracking,'Value')); 
        globalOptions.enablePostproc = logical(get(h_all.cbx_enablePostproc,'Value')); 
        
        % Store options from plugins
        importOptions = importPlugin.getOptions();
        candidateOptions = candidate_plugins(selected_candidate_plugin).getOptions();
        refinementOptions = refinement_plugins(selected_refinement_plugin).getOptions();
        trackingOptions = tracking_plugins(selected_tracking_plugin).getOptions();
        postprocOptions = postproc_plugins(selected_postproc_plugin).getOptions();
        
        % Store if loaded data should be used
        GUIreturns.use_loaded_candidateData  = get(h_all.cbx_candidate_loaded,'Value');
        GUIreturns.use_loaded_refinementData = get(h_all.cbx_refinement_loaded,'Value');
        GUIreturns.use_loaded_trackingData   = get(h_all.cbx_tracking_loaded,'Value');
        GUIreturns.use_loaded_postprocData   = get(h_all.cbx_postproc_loaded,'Value');
        
        %Check if options were changed compared to intial ones
        GUIreturns.globalOptionsChanged     = ~isequaln(globalOptions_atStartup, globalOptions);
        GUIreturns.importOptionsChanged     = ~isequaln(importOptions_atStartup, importOptions);
        GUIreturns.candidateOptionsChanged  = ~isequaln(candidateOptions_atStartup, candidateOptions);
        GUIreturns.refinementOptionsChanged = ~isequaln(refinementOptions_atStartup, refinementOptions);
        GUIreturns.trackingOptionsChanged   = ~isequaln(trackingOptions_atStartup, trackingOptions);
        GUIreturns.postprocOptionsChanged   = ~isequaln(postprocOptions_atStartup, postprocOptions);
        %Check if preview window changed
        GUIreturns.previewIntervalChanged = (globalOptions_atStartup.firstFrameTesting ~= globalOptions.firstFrameTesting) || (globalOptions_atStartup.lastFrameTesting ~= globalOptions.lastFrameTesting) || (globalOptions_atStartup.binFrame ~= globalOptions.binFrame);
        GUIreturns.timegateChanged = (globalOptions_atStartup.useTimegate ~= globalOptions.useTimegate) || (globalOptions_atStartup.tgStart ~= globalOptions.tgStart) || (globalOptions_atStartup.tgEnd ~= globalOptions.tgEnd);
        
        % Check if only change in globalOptions is enableTracking
        GUIreturns.globalOptionsChanged_ExcludingEnable = false;
        if GUIreturns.globalOptionsChanged
            globalOptions_tmp = globalOptions;
            globalOptions_tmp.enableCandidate = globalOptions_atStartup.enableCandidate;
            globalOptions_tmp.enableRefinement = globalOptions_atStartup.enableRefinement;
            globalOptions_tmp.enableTracking = globalOptions_atStartup.enableTracking;
            globalOptions_tmp.enablePostproc = globalOptions_atStartup.enablePostproc;
            
            GUIreturns.globalOptionsChanged_ExcludingEnable = ~isequaln(globalOptions_atStartup, globalOptions_tmp);
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

% Called when closing the application via the 'X' button (or via close).
% Use try to make sure that this is not blocking in case of errors.
    function onAppClose(hObj, event)
        try
            storeOptions();
            GUIreturns.userExit = true;
        catch err
            warning(err.message);
        end
        try
            delete(h_main);
        catch err
            warning(err.message);
        end
    end

end
