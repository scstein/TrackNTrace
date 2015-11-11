function [generalOptions, candidateOptions,fittingOptions,trackingOptions, GUIreturns] = settingsGUI(generalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs)
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
GUIreturns.generalOptionsChanged = false;
GUIreturns.candidateOptionsChanged = false;
GUIreturns.fittingOptionsChanged = false;
GUIreturns.trackingOptionsChanged  = false;

% Save options at startup (to check later if options changed)
generalOptions_atStartup = generalOptions;
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

% % General options
% edit_movieList
% edit_darkMovie
set(h_all.edit_firstFrame,'Callback',{@callback_IntEdit,1,inf});
set(h_all.edit_lastFrame,'Callback',{@callback_IntEdit,1,inf});
set(h_all.cbx_previewMode, 'Callback', @callback_updateGUIstate);
set(h_all.edit_firstFrameTesting,'Callback',{@callback_IntEdit,1,inf});
set(h_all.edit_lastFrameTesting,'Callback',{@callback_IntEdit,1,inf});
set(h_all.button_movieList, 'Callback', @callback_selectMovieList);
set(h_all.button_darkMovie, 'Callback', @callback_selectDarkMovie);

%Photon conversion
set(h_all.cbx_usePhotonConv, 'Callback', @callback_updateGUIstate);
set(h_all.edit_photonBias, 'Callback', {@callback_IntEdit,0,inf});
set(h_all.edit_photonSensitivity, 'Callback', {@callback_FloatEdit,1.0,inf});
set(h_all.edit_photonGain, 'Callback', {@callback_IntEdit,1,1000});

% % Candidate Options
set(h_all.popup_candidateMethod, 'Callback', @callback_updateGUIstate);
% popup_fitDirection
set(h_all.cbx_calcOnce, 'Callback', @callback_updateGUIstate);
set(h_all.edit_avgWinSize, 'Callback', {@callback_IntEdit,1,inf});

% % Tracking
set(h_all.cbx_enableTracking, 'Callback', @callback_updateGUIstate);
set(h_all.popup_trackingMethod, 'Callback', @callback_updateGUIstate);
% cbx_verbose

% % Parse plugins
% Save plugin info (function / name)
candidate_plugins = {};
fitting_plugins = {};
tracking_plugins = {};

% Save last selected plugin per category
selected_candidate_plugin = -1; 
selected_fitting_plugin = -1; 
selected_tracking_plugin = -1;

% Save options for every plugin
candidate_plugin_options = {};
fitting_plugin_options = {};
tracking_plugin_options = {};

% Load the plugins
loadPlugins();

% % GUI main
setGUIBasedOnOptions();
uiwait(h_main);
drawnow; % makes figure disappear instantly (otherwise it looks like it is existing until script finishes)

    % Load plugins from the plugins folder
    function loadPlugins()
        fullPathToThisFile = mfilename('fullpath');
        [folderPath,~,~] = fileparts(fullPathToThisFile);
        plugin_files = dir([folderPath filesep '..' filesep 'plugins' filesep 'plugin_*.m']);
        nr_plugins = numel(plugin_files);
        
        for iPlug = 1:nr_plugins
            [~, plugin_function, ~] = fileparts(plugin_files(iPlug).name);
            plugin_function = str2func(plugin_function); % Convert string to function handle
            [plugin_name, plugin_type] = plugin_function();
            
            switch plugin_type
                case 1 % Candidate method
                    candidate_plugins = [candidate_plugins; {plugin_function, plugin_name}];
                case 2 % Fitting method
                    fitting_plugins  = [fitting_plugins; {plugin_function, plugin_name}];
                case 3 % Tracking method
                    tracking_plugins = [tracking_plugins; {plugin_function, plugin_name}];
                otherwise
                    warning('Detected unknown plugin of type %i',plugin_type);
            end
        end
        
        found_candidate_plugin = size(candidate_plugins,1)>0;
        found_fitting_plugin = size(fitting_plugins,1)>0;
        found_tracking_plugin = size(tracking_plugins,1)>0;
        
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
        set(h_all.popup_candidateMethod, 'String', candidate_plugins(:,2));
        set(h_all.popup_fittingMethod, 'String', fitting_plugins(:,2));
        set(h_all.popup_trackingMethod, 'String', tracking_plugins(:,2)); 
        
        % Allocate size to store options
        candidate_plugin_options = cell(size(candidate_plugins,1),1);
        fitting_plugin_options = cell(size(fitting_plugins,1),1);
        tracking_plugin_options = cell(size(tracking_plugins,1),1);
        
        setPluginsBasedOnOptions();
    end

    function setPluginsBasedOnOptions()
       % Select plugin based on the input options
        if isfield(candidateOptions,'plugin_name')
           selected_candidate_plugin = -1; 
           for iPlug = 1:size(candidate_plugins,1)
                if strcmp(candidateOptions.plugin_name,candidate_plugins{iPlug,2})
                    selected_candidate_plugin = iPlug;
                end
           end
           
           if selected_candidate_plugin < 0
              error('Plugin ''%s'' not found.',candidateOptions.plugin_name) ;
           end
        else
            selected_candidate_plugin = 1; 
        end
        
        
        if isfield(fittingOptions,'plugin_name')
           selected_fitting_plugin = -1; 
           for iPlug = 1:size(fittingOptions,1)
                if strcmp(fittingOptions.plugin_name,fitting_plugins{iPlug,2})
                    selected_fitting_plugin = iPlug;
                end
           end
           
           if selected_fitting_plugin < 0
              error('Plugin ''%s'' not found.',fittingOptions.plugin_name) ;
           end
        else
            selected_fitting_plugin = 1; 
        end
        
        if isfield(trackingOptions,'plugin_name')
           selected_tracking_plugin = -1; 
           for iPlug = 1:size(trackingOptions,1)
                if strcmp(trackingOptions.plugin_name,tracking_plugins{iPlug,2})
                    selected_tracking_plugin = iPlug;
                end
           end
           
           if selected_tracking_plugin < 0
              error('Plugin ''%s'' not found.',trackingOptions.plugin_name) ;
           end
        else
            selected_tracking_plugin = 1; 
        end
        
        % Set popups to correct value
        set(h_all.popup_candidateMethod, 'Value', selected_candidate_plugin);
        set(h_all.popup_fittingMethod, 'Value', selected_fitting_plugin);
        set(h_all.popup_trackingMethod, 'Value', selected_tracking_plugin);
        
        % Build panels
        candidate_plugins{ selected_candidate_plugin,1}(h_all.panel_candidate, candidateOptions);
        fitting_plugins{ selected_fitting_plugin,1}(h_all.panel_fitting, fittingOptions);
        tracking_plugins{ selected_tracking_plugin,1}(h_all.panel_tracking, trackingOptions); 
    end

% Update GUI based on currently set values
    function callback_updateGUIstate(hObj, event)
        % Enable/disable testing fields
        if get(h_all.cbx_previewMode, 'Value');
            set(h_all.edit_firstFrameTesting, 'Enable','on');
            set(h_all.edit_lastFrameTesting, 'Enable','on');
        else
            set(h_all.edit_firstFrameTesting, 'Enable','off');
            set(h_all.edit_lastFrameTesting, 'Enable','off');
        end
        
        if get(h_all.cbx_calcOnce, 'Value')
            set(h_all.edit_avgWinSize, 'Enable','on');
        else
            set(h_all.edit_avgWinSize, 'Enable','off');
        end
        
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
            set(h_all.edit_movieList,'Enable', 'off');
            set(h_all.button_movieList,'Enable', 'off');
            set(h_all.button_continueForAll, 'Visible','off');
        end
        
        % % Plugins
        % Store options for last selected plugins
       candidate_plugin_options{selected_candidate_plugin} = getappdata(h_all.panel_candidate,'options');
       fitting_plugin_options{selected_fitting_plugin} = getappdata(h_all.panel_fitting,'options');
       tracking_plugin_options{selected_tracking_plugin} = getappdata(h_all.panel_tracking,'options');
        
        % Update panels to display currently selected plugin
        % This is done by executing the plugin function with the correct panel as the argument
        selected_candidate_plugin =  get(h_all.popup_candidateMethod,'Value');
        selected_fitting_plugin = get(h_all.popup_fittingMethod,'Value');
        selected_tracking_plugin = get(h_all.popup_trackingMethod,'Value');
        
        candidate_plugins{ selected_candidate_plugin,1}(h_all.panel_candidate, candidate_plugin_options{selected_candidate_plugin});
        fitting_plugins{ selected_fitting_plugin,1}(h_all.panel_fitting, fitting_plugin_options{selected_fitting_plugin});
        tracking_plugins{ selected_tracking_plugin,1}(h_all.panel_tracking, tracking_plugin_options{selected_tracking_plugin});
        
        updatePanelPositions();
                
        % Enable/disable tracking panel
        if get(h_all.cbx_enableTracking, 'Value')
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','on');
        else
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','off');
        end
    end

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
        save(outfile,'generalOptions', 'candidateOptions','fittingOptions','trackingOptions');
    end

% Load settings from a file
    function callback_loadSettings(hObj,event)
        % In single file mode, the list of movies to process should not be
        % changed, we store it and restore after loading the settings
        if(GUIinputs.singleFileMode)
            filename_movies = generalOptions.filename_movies;
        end
        
        [infile, path] = uigetfile('.mat');
        if isfloat(infile);
            return;
        end; % User clicked cancel
        % Note: Loading has to be done this way, as variables "can not be
        % added to a static workspace" (e.g. the one of this GUI).
        allOptions = load([path,infile],'generalOptions', 'candidateOptions','fittingOptions','trackingOptions');
        generalOptions   = allOptions.generalOptions;
        if(GUIinputs.singleFileMode)
            generalOptions.filename_movies = filename_movies;
        end
        candidateOptions = allOptions.candidateOptions;
        fittingOptions   = allOptions.fittingOptions;
        trackingOptions  = allOptions.trackingOptions;
        
        setGUIBasedOnOptions();
    end

% Opens a file chooser dialog to choose multiple input (movie) files
% for processing. Note: filenames will be seperated by ';'
    function callback_selectMovieList(hObj,event)
        % Get current text field to set starting path of uigetfile
        movieListString = get(h_all.edit_movieList,'String');
        path = [];
        if ~isempty(movieListString)
            movieList = strsplit(movieListString,';');
            [path,~,~] = fileparts(movieList{1});
        end
        [movieList, path] = uigetfile([path,filesep,'*.tif'],'MultiSelect','on');
        if( isfloat(movieList) ); return; end; % User pressed cancel.
        
        % atach path to every entry in the list
        if iscell(movieList)
            for i=1:length(movieList)
                movieList{i} =  [path,movieList{i}];
            end
        elseif ischar(movieList)
            movieList = [path,movieList];
        end
        
        set(h_all.edit_movieList,'String',cell2str(movieList));
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
% 'isInteger' to true for displaying integer values.
    function setNum(hObj,value,isInteger)
        %         value = num2str(value);
        if nargin<3 || isempty(isInteger)
            isInteger = false;
        end
        
        if isInteger
            set(hObj,'String',sprintf('%i',value));
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
% (generalOptions, candidateOptions, fittingOptions, trackingOptions)
    function setGUIBasedOnOptions()
        % % General Options
        set(h_all.edit_movieList,'String', cell2str(generalOptions.filename_movies));
        set(h_all.edit_darkMovie,'String', generalOptions.filename_dark_movie);
        setNum(h_all.edit_firstFrame, generalOptions.firstFrame, true);
        setNum(h_all.edit_lastFrame, generalOptions.lastFrame, true);
        set(h_all.cbx_previewMode, 'Value', generalOptions.previewMode);
        setNum(h_all.edit_firstFrameTesting, generalOptions.firstFrameTesting, true);
        setNum(h_all.edit_lastFrameTesting, generalOptions.lastFrameTesting, true);
        
        % Candidate Options
        if generalOptions.fitForward
            setPopup(h_all.popup_fitDirection,'Forward');
        else
            setPopup(h_all.popup_fitDirection,'Backward');
        end
        
        % Calc once
        set(h_all.cbx_calcOnce, 'Value', generalOptions.calculateCandidatesOnce);
        setNum(h_all.edit_avgWinSize, generalOptions.averagingWindowSize, true);
        
        % Photon conversion
        set(h_all.cbx_usePhotonConv,'Value',generalOptions.usePhotonConversion)
        set(h_all.edit_photonBias,'Value',generalOptions.photonBias);
        set(h_all.edit_photonSensitivity,'Value',generalOptions.photonSensitivity);
        set(h_all.edit_photonGain,'Value',generalOptions.photonGain);
        
        % % Tracking
        set(h_all.cbx_enableTracking,'Value', generalOptions.enableTracking); % --
        
        % Update plugins
        setPluginsBasedOnOptions();
        
        % Update the GUI
        callback_updateGUIstate();
    end

% Set all options structs (generalOptions, candidateOptions, fittingOptions, trackingOptions)
% based on the current state of the uicontrols. Also checks which
% options structs have changed compared to the startup values.
    function storeOptions()
        % % General Options
        generalOptions.filename_movies = strsplit(get(h_all.edit_movieList,'String'), ';');
        generalOptions.filename_dark_movie = get(h_all.edit_darkMovie,'String');
        generalOptions.firstFrame = getNum(h_all.edit_firstFrame);
        generalOptions.lastFrame =  getNum(h_all.edit_lastFrame);
        generalOptions.previewMode = get(h_all.cbx_previewMode, 'Value');
        generalOptions.firstFrameTesting = getNum(h_all.edit_firstFrameTesting);
        generalOptions.lastFrameTesting = getNum(h_all.edit_lastFrameTesting);
        
        generalOptions.fitForward = logical(strcmp(getPopup(h_all.popup_fitDirection), 'Forward'));
        generalOptions.calculateCandidatesOnce = logical(get(h_all.cbx_calcOnce, 'Value'));
        generalOptions.averagingWindowSize = getNum(h_all.edit_avgWinSize);
        
        generalOptions.usePhotonConversion = logical(get(h_all.cbx_usePhotonConv,'Value'));
        generalOptions.photonBias = getNum(h_all.edit_photonBias);
        generalOptions.photonSensitivity = getNum(h_all.edit_photonSensitivity);
        generalOptions.photonGain = getNum(h_all.edit_photonGain);
        
        % Tracking
        generalOptions.enableTracking = logical(get(h_all.cbx_enableTracking,'Value')); % --
        
        %Check if options were changed compared to intial ones
        GUIreturns.generalOptionsChanged   = ~isequaln(generalOptions_atStartup, generalOptions);
        GUIreturns.candidateOptionsChanged = ~isequaln(candidateOptions_atStartup, candidateOptions);
        GUIreturns.fittingOptionsChanged   = ~isequaln(fittingOptions_atStartup, fittingOptions);
        GUIreturns.trackingOptionsChanged  = ~isequaln(trackingOptions_atStartup, trackingOptions);
        %Check if preview window changed
        GUIreturns.testWindowChanged = (generalOptions_atStartup.firstFrameTesting ~= generalOptions.firstFrameTesting) || (generalOptions_atStartup.lastFrameTesting ~= generalOptions.lastFrameTesting);
        
        % Store options from plugins
        candidateOptions = getappdata(h_all.panel_candidate,'options');
        fittingOptions = getappdata(h_all.panel_fitting,'options');
        trackingOptions = getappdata(h_all.panel_tracking,'options');
    end


% Cleanup function. This is neccessary to delete the timer!
    function callback_continue(hObj, event)
        storeOptions(); % Store the options before closing
        delete(h_main);
    end

% Called when closing the application via the 'X' button (or via close)
    function onAppClose(hObj, event)
        storeOptions();
        GUIreturns.userExit = true;
        delete(h_main);
    end

end
