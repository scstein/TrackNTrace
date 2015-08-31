function [generalOptions, candidateOptions,fittingOptions,trackingOptions, GUIreturns] = settingsGUI(generalOptions, candidateOptions,fittingOptions,trackingOptions, GUIinputs)
% USAGE:
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
GUIreturns.trackingOptionsChanged  = false;

% Save options at startup
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
set(h_all.edit_firstFrame,'Callback',{@callback_positiveIntEdit,1,inf});
set(h_all.edit_lastFrame,'Callback',{@callback_positiveIntEdit,1,inf});
set(h_all.cbx_previewMode, 'Callback', @callback_updateGUIstate);
set(h_all.edit_firstFrameTesting,'Callback',{@callback_positiveIntEdit,1,inf});
set(h_all.edit_lastFrameTesting,'Callback',{@callback_positiveIntEdit,1,inf});
set(h_all.button_movieList, 'Callback', @callback_selectMovieList);
set(h_all.button_darkMovie, 'Callback', @callback_selectDarkMovie);

% % Candidate Options
set(h_all.popup_candidateMethod, 'Callback', @callback_updateGUIstate);
set(h_all.edit_stddev, 'Callback', {@callback_positiveFloatEdit,0,inf});
% popup_fitDirection
set(h_all.cbx_calcOnce, 'Callback', @callback_updateGUIstate);
set(h_all.edit_avgWinSize, 'Callback', {@callback_positiveIntEdit,1,inf});
% Only for correlation
set(h_all.edit_corrThreshold,'Callback', {@callback_positiveFloatEdit,0,1});
% Only for intensity filtering
set(h_all.edit_particleRadius,'Callback', {@callback_positiveIntEdit,1,inf});
set(h_all.edit_intensityThreshold,'Callback', {@callback_positiveFloatEdit,0,100});
set(h_all.edit_pTest,'Callback', {@callback_positiveFloatEdit,0,1});
set(h_all.edit_bgInterval,'Callback', {@callback_positiveIntEdit,1,inf});


% % Fitting options
% cbx_fitStddev
% cbx_usePixelIntegration
% cbx_useMleRefinement

% % Tracking
set(h_all.cbx_enableTracking, 'Callback', @callback_updateGUIstate);
set(h_all.popup_trackerMethod, 'Callback', @callback_updateGUIstate);
set(h_all.edit_trackerRadius,'Callback', {@callback_positiveFloatEdit,0,inf});
set(h_all.edit_maxGap,'Callback', {@callback_positiveIntEdit,0,inf});
% popup_linkingMethod % Only for simpletracker
set(h_all.edit_minTrackLength,'Callback',{@callback_positiveIntEdit,1,inf});
set(h_all.edit_splitMovieParts, 'Callback', {@callback_positiveIntEdit,1,inf});
% cbx_verbose

% GUI main
setOptions();
uiwait(h_main);
drawnow; % makes figure disappear instantly (otherwise it looks like it is existing until script finishes)

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
        
        % Show only candidate options matching the method
        if strcmp(getPopup(h_all.popup_candidateMethod), 'Cross-correlation')
            set(h_all.panel_intensityFiltering,'Visible','off')
            set(h_all.panel_correlation,'Visible','on');
        else
            set(h_all.panel_intensityFiltering,'Visible','on')
            set(h_all.panel_correlation,'Visible','off');
        end
        
        if get(h_all.cbx_calcOnce, 'Value')
            set(h_all.edit_avgWinSize, 'Enable','on');
        else
            set(h_all.edit_avgWinSize, 'Enable','off');
        end
        
        % Enable/disable tracking panel
        if get(h_all.cbx_enableTracking, 'Value')
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','on');
            
            % Disable linking method except for simpletracker
            if ~strcmp(getPopup(h_all.popup_trackerMethod), 'simpletracker')
                set(h_all.popup_linkingMethod,'Enable','off')
            end
        else
            set(findall(h_all.panel_tracking,'-property','Enable'), 'Enable','off');
        end
        
        % In single file mode we disable choosing the movie list
        if(GUIinputs.singleFileMode)
            set(h_all.edit_movieList,'Enable', 'off');
            set(h_all.button_movieList,'Enable', 'off');
            set(h_all.button_continueForAll, 'Visible','off');
        end

    end

    function callback_continueForAll(hObj,event)
        GUIreturns.useSettingsForAll = true;
        callback_continue(hObj,event); 
    end

    function callback_saveSettings(hObj, event)
        storeOptions();
        
        [outfile,path] = uiputfile('.mat');
        if isfloat(outfile); return; end; % User clicked cancel
        
        outfile = [path,outfile];
        save(outfile,'generalOptions', 'candidateOptions','fittingOptions','trackingOptions');
    end

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
       setOptions();
    end

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

    function callback_positiveFloatEdit(hObj,event, minVal, maxVal)
        if nargin<3 || isempty(minVal);
            minVal=0; 
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

    function callback_positiveIntEdit(hObj,event, minVal,maxVal)
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


    function value = getNum(hObj)
        value = str2num(get(hObj,'String'));
    end

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

    function optionString = getPopup(hObj)
        choices = get(hObj,'String');
        optionString = choices{get(hObj,'Value')};
    end

    function setOptions()
        % % General Options
        set(h_all.edit_movieList,'String', cell2str(generalOptions.filename_movies));
        set(h_all.edit_darkMovie,'String', generalOptions.filename_dark_movie);
        setNum(h_all.edit_firstFrame, generalOptions.firstFrame, true);
        setNum(h_all.edit_lastFrame, generalOptions.lastFrame, true);
        set(h_all.cbx_previewMode, 'Value', generalOptions.previewMode);
        setNum(h_all.edit_firstFrameTesting, generalOptions.firstFrameTesting, true);
        setNum(h_all.edit_lastFrameTesting, generalOptions.lastFrameTesting, true);
        
        % Candidate Options
        if candidateOptions.useCrossCorrelation
            setPopup(h_all.popup_candidateMethod,'Cross-correlation');
        else
            setPopup(h_all.popup_candidateMethod,'Intensity filtering');
        end
        
        if candidateOptions.fitForward
            setPopup(h_all.popup_fitDirection,'Forward');
        else
            setPopup(h_all.popup_fitDirection,'Backward');
        end
        setNum(h_all.edit_stddev, candidateOptions.sigma);
        set(h_all.cbx_calcOnce, 'Value', candidateOptions.calculateCandidatesOnce);
        setNum(h_all.edit_avgWinSize, candidateOptions.averagingWindowSize, true);
        % Only for correlation
        setNum(h_all.edit_corrThreshold, candidateOptions.corrThresh);
        % Only for intensity filtering
        setNum(h_all.edit_particleRadius, candidateOptions.particleRadius, true);
        setNum(h_all.edit_intensityThreshold,candidateOptions.intensityThreshold, true);
        setNum(h_all.edit_pTest, candidateOptions.intensityPtestVar);
        setNum(h_all.edit_bgInterval,candidateOptions.backgroundCalculationInterval, true);    %
        
        % Fitting options
        set(h_all.cbx_fitStddev,'Value', fittingOptions.fitSigma)
        set(h_all.cbx_usePixelIntegration,'Value', fittingOptions.usePixelIntegratedFit)
        set(h_all.cbx_useMleRefinement,'Value',fittingOptions.useMLErefine)
        
        % % Tracking
        set(h_all.cbx_enableTracking,'Value', trackingOptions.enableTracking); % --
        setPopup(h_all.popup_trackerMethod, trackingOptions.method);
        set(h_all.cbx_verbose, 'Value', trackingOptions.verbose);
        setNum(h_all.edit_trackerRadius, trackingOptions.maxRadius);
        setNum(h_all.edit_maxGap, trackingOptions.maxGap, true);
        setPopup(h_all.popup_linkingMethod, trackingOptions.linkingMatrix) % Only for simpletracker
        setNum(h_all.edit_minTrackLength, trackingOptions.minTrackLength, true);
        setNum(h_all.edit_splitMovieParts, trackingOptions.splitMovieParts, true);  
        
        % Update the GUI
        callback_updateGUIstate();
    end

 function storeOptions()
        % % General Options
        generalOptions.filename_movies = strsplit(get(h_all.edit_movieList,'String'), ';');
        generalOptions.filename_dark_movie = get(h_all.edit_darkMovie,'String');
        generalOptions.firstFrame = getNum(h_all.edit_firstFrame);
        generalOptions.lastFrame =  getNum(h_all.edit_lastFrame);
        generalOptions.previewMode = get(h_all.cbx_previewMode, 'Value');
        generalOptions.firstFrameTesting = getNum(h_all.edit_firstFrameTesting);
        generalOptions.lastFrameTesting = getNum(h_all.edit_lastFrameTesting);
        
        %             % % Candidate Options
        candidateOptions.useCrossCorrelation = strcmp(getPopup(h_all.popup_candidateMethod), 'Cross-correlation'); 
        candidateOptions.fitForward = strcmp(getPopup(h_all.popup_fitDirection), 'Forward'); 
        candidateOptions.sigma = getNum(h_all.edit_stddev);
        candidateOptions.calculateCandidatesOnce = get(h_all.cbx_calcOnce, 'Value');
        candidateOptions.averagingWindowSize = getNum(h_all.edit_avgWinSize);
        % Only for correlation
        candidateOptions.corrThresh = getNum(h_all.edit_corrThreshold);
        % Only for intensity filtering
        candidateOptions.particleRadius = getNum(h_all.edit_particleRadius);
        candidateOptions.intensityThreshold = getNum(h_all.edit_intensityThreshold);
        candidateOptions.intensityPtestVar = getNum(h_all.edit_pTest);
        candidateOptions.backgroundCalculationInterval = getNum(h_all.edit_bgInterval);
        
        % Fitting options
        fittingOptions.fitSigma = get(h_all.cbx_fitStddev,'Value');
        fittingOptions.usePixelIntegratedFit = get(h_all.cbx_usePixelIntegration,'Value');
        fittingOptions.useMLErefine = get(h_all.cbx_useMleRefinement,'Value');
        
        % % Tracking
        trackingOptions.enableTracking = get(h_all.cbx_enableTracking,'Value'); % --
        trackingOptions.method = getPopup(h_all.popup_trackerMethod);
        trackingOptions.verbose = logical(get(h_all.cbx_verbose, 'Value'));
        trackingOptions.maxRadius = getNum(h_all.edit_trackerRadius);
        trackingOptions.maxGap = getNum(h_all.edit_maxGap);
        trackingOptions.linkingMatrix = getPopup(h_all.popup_linkingMethod); % Only for simpletracker
        trackingOptions.minTrackLength = getNum(h_all.edit_minTrackLength);
        trackingOptions.splitMovieParts = getNum(h_all.edit_splitMovieParts);        
        
        
        %Check if options were changed compared to intial ones
        GUIreturns.generalOptionsChanged   = ~isequaln(generalOptions_atStartup, generalOptions);
        GUIreturns.candidateOptionsChanged = ~isequaln(candidateOptions_atStartup, candidateOptions);
        GUIreturns.fittingOptionsChanged   = ~isequaln(fittingOptions_atStartup, fittingOptions);
        GUIreturns.trackingOptionsChanged  = ~isequaln(trackingOptions_atStartup, trackingOptions);
        %Check if preview window changed
        GUIreturns.testWindowChanged = (generalOptions_atStartup.firstFrameTesting ~= generalOptions.firstFrameTesting) || (generalOptions_atStartup.lastFrameTesting ~= generalOptions.lastFrameTesting);
    end


% Cleanup function. This is neccessary to delete the timer!
    function callback_continue(hObj, event)
        storeOptions(); % Store the options before closing
        delete(h_main);
    end

    function onAppClose(hObj, event)
       storeOptions();
       GUIreturns.userExit = true;
       delete(h_main);
    end

end
