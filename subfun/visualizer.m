function [h_main] = visualizer(movie, candidateData, candidateParams, fitData, fitParams, trackingData, trackingParams, FPS, is_blocking)
% USAGE: visualizeFitDataGUI(movie, trackingData)
% [ Full USAGE: visualizeFitDataGUI(movie, fitData, FPS, use_bw, blockUI) ]
%
% Visualizer for TrackNTrace results.
%
% Input:
%   movie: 3D matrix (rows,cols,frames) of the analyzed movie
%   fitData: 3D array of Gaussian distribution parameters
%     (param,particle,frame) for all found particles. The line order is
%     [x;y;A;B;sigma;flag], the column order is the particle index and the
%     slice order is the movie frame. The positions x (column) and y (row)
%     are not corrected by a middle pixel shift, the center of the top left
%     pixel is [1.0,1.0]. A is the unnormalized amplitude of a Gaussian
%     distribution A*exp(...)+B with constant background B. flag is the
%     exit flag of the fitting routine, see psfFit_Image.m for details.
%   FPS: frames per second to play movie with | default: 30
%   use_bw: black/white image, otherwise colormap hot | default false
%   is_blocking: Blocks MATLAB execution while visualizer is open
%
%  Inputs (except movie, fitData) can be left empty [] for default values.
%
% Output:
%   h_main - Handle to the GUI figure
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2015
%

% Parse given inputs. For clarity we outsource this in a function.
if nargin<1 || isempty(movie)
    fprintf(' Need input movie!\n');
    return
end

% This variablea are computed on demand if the user wants to plot
% distributions of parameters
allFramesCandidateData = [];
allFramesFitData = [];

% Shared variables
use_bw = false;
mode = 'none';
traj_lifetime = 0;
n_colors = 20;
traj_displayLength = inf;

% Prepare tracking data
if isempty(trackingData)
    id_tracks = []; % Note: (in case track IDs go from 1 to N without missing numbers, the track index is identical to the tracks ID.
    n_tracks = 0;
else
    id_tracks = unique(trackingData(:,1)); % Note: (in case track IDs go from 1 to N without missing numbers, the track index is identical to the tracks ID.
    n_tracks = numel(id_tracks);
end

% Put data for each track into its own cell. This is much more efficient comfortable for plotting.
cell_traj = cell(n_tracks ,1);
cnt = 1;
for iTrack = 1:n_tracks
    cell_traj{iTrack} = trackingData( trackingData(:,1)== id_tracks(cnt) , 2:end);
    cnt = cnt+1;
end



% -- Preparing the GUI --
h_main = openfig('visualizer.fig');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % Executed on closing for cleanup
set(h_main,'Toolbar','figure');   % Add toolbar needed for zooming
set(h_main, 'DoubleBuffer', 'on') % Helps against flickering
movegui(h_main,'center');

h_all = guihandles(h_main); % Get handles of all GUI objects
axes(h_all.axes); % Select axis for drawing plots

dcm_obj = datacursormode(h_main);
defaultDatatipFunction = get(dcm_obj,'UpdateFcn');

parse_inputs_and_setup(nargin);


%  -- Setup UI elements --
% Text on top
set(h_all.toptext,'String',sprintf('Current frame: 1/%i',size(movie,3)));

% Buttons
set(h_all.button_candidateMode,'Callback', {@callback_changeMode,'candidate'});
set(h_all.button_fittingMode,'Callback', {@callback_changeMode,'fitting'});
set(h_all.button_trackingMode,'Callback', {@callback_changeMode,'tracking'});

set(h_all.but_play,'Callback',@playCallback);
set(h_all.but_contrast,'Callback',@contrastCallback);
set(h_all.but_autocontrast,'Callback',@autocontrastCallback);
set(h_all.but_distribution,'Callback',@distributionCallback);

% Slider
set(h_all.slider,'Value',1, 'Min',1,'Max',size(movie,3),'SliderStep',[1/size(movie,3) 1/size(movie,3)],'Callback', @sliderCallback);
hLstn = addlistener(h_all.slider,'ContinuousValueChange',@updateSlider); %#ok<NASGU> % Add event listener for continous update of the shown slider value

% Edit fields
set(h_all.edit_FPS,'String',sprintf('%i',FPS), 'Callback', @fpsCallback);
set(h_all.edit_distributionBins, 'Callback', {@callback_intEdit,1,inf});
setNum(h_all.edit_distributionBins, 50, true);
set(h_all.edit_distributionRange,'Callback',{@callback_FloatEdit,0,100});
setNum(h_all.edit_distributionRange,100);

% Checkbox
set(h_all.cb_bw, 'Value', use_bw, 'Callback',@bwCallback);

% Popupmenu
set(h_all.popup_distribution, 'String', fitParams);

%  -- Candidate UI elements --

%  -- Fitting UI elements --

%  -- Tracking UI elements --
set(h_all.edit_lifetime,'String',sprintf('%i',traj_lifetime), 'Callback', @lifetimeCallback);
set(h_all.edit_colors,'String',sprintf('%i',n_colors), 'Callback', @colorCallback);
set(h_all.edit_trajDisplayLength,'String',sprintf('%i', traj_displayLength), 'Callback', @dispLengthCallback);

% Timer -> this controls playing the movie
h_all.timer = timer(...
    'ExecutionMode', 'fixedDelay', ...    % Run timer repeatedly
    'Period', round(1/FPS*1000)/1000, ... % Initial period is 1 sec. Limited to millisecond precision
    'TimerFcn', @onTimerUpdate, ...
    'StartFcn', @onTimerStart, ...
    'StopFcn',  @onTimerStop); % Specify callback

% Store handles to the plot objects (which is faster)
% Handles are set on first use (mostly in plotFrame)
linehandles = -1*ones(n_tracks,1);
dothandle_fit = -1;
dothandle_cand = -1;
imagehandle = -1;

% Draw the marker color depending on background color
track_colors = [];
marker_color = [];
drawColors(n_colors);

% Plot first frame to get limits right
% Set x,y,color limits
xl = [0.5,size(movie,2)+0.5];
yl = [0.5,size(movie,1)+0.5];
firstImg = movie(:,:,1);
zl = [min(firstImg(:)), max(firstImg(:))];

plotFrame(1);

xlim(xl);
ylim(yl);
caxis(zl);

% Calling this creates handles for all tracks not plotted before.
% Although this takes some time, the visualizer will respond smoother
% afterwards. If you uncomment this function, the visualizer starts faster,
% but may show jerky behaviour during first play as handles for tracks that
% did not occur before are created on demand.
setUnitinializedTrackHandles();

% Variables for playback
timePerFrame = round(1/FPS*1000)/1000; % limit to millisecond precision
elapsed_time = 0;
frame = 1;

% Change into the right mode (candidate/fitting/tracking)
callback_changeMode();

% In case the RunAgain dialog should be displayed, we stop scripts/functions
% calling the GUI until the figure is closed
if(is_blocking)
    uiwait(h_main);
    drawnow; % makes figure disappear instantly (otherwise it looks like it is existing until script finishes)
end


% --- Nested Functions ---

    % Change the chosen mode 'candidate','fitting','tracking' and display
    % its relevant content. 
    % If "modus" input is given, the mode is set to "modus". Its
    % implemented this way to use one callback for all buttons selecting the modes.
    function callback_changeMode(hObj,event,modus)
        if nargin>2
            mode = modus;
        end
        DEFAULT_COLOR = [0.941,0.941,0.941]; % Default color of buttons.
        SELECTED_COLOR = [0.65, 0.9, 0]; % Color of selected button.
        
        %Reset button colors       
        set(h_all.button_candidateMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_fittingMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_trackingMode,'BackgroundColor', DEFAULT_COLOR);
                       
        % Reset handles
        if dothandle_fit ~= -1
            set(dothandle_fit,'xdata',[],'ydata',[]);
        end
        if dothandle_cand ~= -1
            set(dothandle_cand,'xdata',[],'ydata',[]);
        end
        for iTr = 1:n_tracks
            if(linehandles(iTr)~=-1)
                set(linehandles(iTr),'xdata',[],'ydata',[]);
            end
        end
        
        % Mode specific changes (setting datatip function, highlight
        % button.
        switch mode
            case 'candidate'
                set(dcm_obj,'UpdateFcn',defaultDatatipFunction);
                set(h_all.button_candidateMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', candidateParams);
            case 'fitting'
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_fittingMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', fitParams);
            case 'tracking'
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_trackingMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', trackingParams);
            otherwise
                error('Unkown mode ''%s''!', mode);
        end
        
        % Update size of GUI, show mode specific panels
        updateGUIforMode();
        
        % Replot
        updateFrameDisplay();
    end


    % Update size of GUI, show mode specific panels
    function updateGUIforMode()
        units = 'characters';
        BOTTOM_SPACING = 0.5;
        
        % Set all mode specific panels invisible
        set(h_all.panel_tracking,'Visible','off');
        
        % Rescale window based on last element
        switch mode
            case 'candidate'
                set(h_all.panel_histogram,'Units',units);
                pos = get(h_all.panel_histogram,'Position');
            case 'fitting'
                set(h_all.panel_histogram,'Units',units);
                pos = get(h_all.panel_histogram,'Position');
            case 'tracking'
                set(h_all.panel_tracking,'Units',units);
                pos = get(h_all.panel_tracking,'Position');
                set(h_all.panel_tracking,'Visible','on');
        end
        
        set(h_main,'Units',units);
        win_pos = get(h_main,'Position');
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
        set(h_main,'Position', win_pos);
        
        % Reset units back to normalized, so figure resizes "properly"
        % (cough..)
        set(h_main,'Units','normalized');
        for iObj = 1:numel(all_uiObjects)
            set(all_uiObjects(iObj),'Units','normalized');
        end
    end

% Function for datacursor
    function txt = modeSpecificDatatipFunction(~,event_obj)
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        graphObjHandle = get(event_obj,'Target'); % The target object (line/image) of the cursor
        I = get(event_obj, 'DataIndex');
        
        if(numel(I) == 1) % Plotted position is selected
            txt = {};
            switch mode
                case 'candidate'
                    % Plot all parameters available for that spot in the datacursor window
                    for iPar=1:numel(candidateParams)
                        txt = [txt, {[candidateParams{iPar},': ', num2str(candidateData{frame}(I,iPar))]}];
                    end
                case 'fitting'
                    % Plot all parameters available for that spot in the datacursor window
                    for iPar=1:numel(fitParams)
                        txt = [txt, {[fitParams{iPar},': ', num2str(fitData{frame}(I,iPar))]}];
                    end
                case 'tracking'
                    TrackNr = find(linehandles==graphObjHandle); % Find lineobject for the selected point
                    PointData = cell_traj{TrackNr}(I,:); % Data of the selected point
                    TrackID = sprintf('%i',id_tracks(TrackNr)); % Get track ID from its index (in case TracIDs go from 1 to N without missing numbers, TrackNr==TrackID)
                    % Plot all parameters available for that spot in the datacursor window
                    txt = [txt, {['TrackID: ', TrackID]}];
                    for iPar=2:numel(trackingParams)
                        txt = [txt, {[trackingParams{iPar},': ', num2str(PointData(iPar-1))]}];
                    end
                otherwise
                    error('Unsupported mode ''%s'' for datatip function.',mode)
            end
            
        elseif (numel(I) == 2) % Image is selected
            txt = {['X: ',num2str(pos(1))],...
                ['Y: ',num2str(pos(2))],...
                ['Value: ', num2str(movie(I(2),I(1),frame))]};
        end
    end

% Get numeric value of edit field
    function value = getNum(hObj)
        value = str2num(get(hObj,'String'));
    end

% Set numeric value of edit field
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

% Callback for an edit field taking only values between minVal and maxVal
% Causes the plot to be redrawn
    function callback_FloatEdit_Plus_Update(hObj,event, minVal, maxVal)
        if nargin<3 || isempty(minVal);
            minVal=-inf;
        end
        if nargin<4 || isempty(maxVal);
            maxVal=inf;
        end
        
        % Callback
        callback_FloatEdit(hObj,event,minVal,maxVal);
        
        % Update the display
        updateFrameDisplay();
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
        
        % Check if a valid number was entered
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
    function callback_intEdit(hObj,event, minVal,maxVal)
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


% The main function of the application. This plays the movie if the
% timer is running
    function onTimerUpdate(timer, event)
        % Progress frame counter, clip at length of movie and stop at last
        % frame.
        frame = frame+1;
        if(frame >= size(movie,3))
            frame = size(movie,3);
            updateTopText();
            
            stop(h_all.timer);
        end
        set(h_all.slider,'Value',frame);
        updateTopText()
        
        % Skip frame if computer is too slow drawing
        if elapsed_time > timePerFrame
            elapsed_time = elapsed_time - timePerFrame;
            return;
        end
        tic_start = tic;
        updateFrameDisplay();
        elapsed_time = elapsed_time + toc(tic_start)- timePerFrame;
    end

    function onTimerStart(timer, event)
        set(h_all.but_play,'String','Pause');
        elapsed_time = 0;
    end

    function onTimerStop(timer, event)
        set(h_all.but_play,'String','Play');
    end

% Used to display the current frame as selected by the 'frame' variable
% Also this sets and saves the axis states (e.g. for zooming);
    function updateFrameDisplay()
        % Needed to minimize interference with other figures the user
        % brings into focus. It can be that the images are not plotted to
        % the GUI then but to the selected figure window
        set(0,'CurrentFigure',h_main);
        
        xl = xlim;
        yl = ylim;
        
        plotFrame(frame);
        
        % save the axis limits in case the user zoomed
        xlim(xl);
        ylim(yl);
        caxis(zl);
        
        % Deactivate labels
        set(h_all.axes,'YTickLabel',[]);
        set(h_all.axes,'XTickLabel',[]);
        
        % Adjust contrast continously if shift key is pressed
        modifiers = get(gcf,'currentModifier');
        shiftIsPressed = ismember('shift',modifiers);
        if(shiftIsPressed)
            autocontrastCallback([],[]);
        end
        
        drawnow; % Important! Or Matlab will skip drawing for high FPS
    end

% Plots the frame with the input index
    function plotFrame(iF)
        % Plot the movie frame
        if imagehandle == -1
            imagehandle = imagesc(movie(:,:,iF)); axis image; colormap gray;
        else
            set(imagehandle,'CData',movie(:,:,iF));
        end
        if use_bw
            colormap gray;
        else
            colormap hot;
        end
        
        switch mode
            case 'candidate'
                if(isempty(candidateData{iF}));
                    if dothandle_cand ~= -1
                        set(dothandle_cand,'xdata',[],'ydata',[]);
                    end
                    return
                end % Jump empty frames
                
                hold on;
                if dothandle_cand == -1
                    dothandle_cand = plot(candidateData{iF}(:,1), candidateData{iF}(:,2), 'x','Color',marker_color,'MarkerSize',8);
                else
                    set(dothandle_cand,'xdata',candidateData{iF}(:,1),'ydata',candidateData{iF}(:,2));
                end
                hold off;
            case 'fitting'
                if(isempty(fitData{iF}));
                    if dothandle_fit ~= -1
                        set(dothandle_fit,'xdata',[],'ydata',[]);
                    end
                    return
                end % Jump empty frames
                
                hold on;
                if dothandle_fit == -1
                    dothandle_fit = plot(fitData{iF}(:,1), fitData{iF}(:,2), 'o','Color',marker_color);
                else
                    set(dothandle_fit,'xdata',fitData{iF}(:,1),'ydata',fitData{iF}(:,2));
                end
                hold off;
            case 'tracking'
                % Draw the tracks of currently visible particles
                hold on;
                for iTr = 1:n_tracks
                    % Don't draw tracks not yet visible (iF <...) or not visible any more
                    % (iF > ...). If traj_lifetime>0 the tracks are displayed for
                    % the given number of frames longer.
                    if iF < cell_traj{iTr}(1,1) || iF > cell_traj{iTr}(end,1) + traj_lifetime
                        mask_toPlot = false(size(cell_traj{iTr},1),1);
                    else
                        % Plot trajectories a) only the last traj_displayLength positoins AND  b) up to the current frame
                        mask_toPlot = ((cell_traj{iTr}(:,1)>iF-traj_displayLength) & cell_traj{iTr}(:,1)<=iF);
                    end
                    
                    % If this trajectory was already plotted before, we just set
                    % its data via its handle (which is fast!). If not, create a
                    % new lineseries by using the plot command.
                    if (linehandles(iTr) == -1)
                        if( sum(mask_toPlot) ~= 0) % only plot the first time we have actual data to display
                            linehandles(iTr) = plot(cell_traj{iTr}(mask_toPlot, 2), cell_traj{iTr}(mask_toPlot, 3), '.--','Color',track_colors(iTr,:));
                        end
                    else
                        set(linehandles(iTr),'xdata',cell_traj{iTr}(mask_toPlot, 2),'ydata', cell_traj{iTr}(mask_toPlot, 3));
                    end
                end
                hold off;
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end
    end

% Switch play/pause by button
    function playCallback(hObj, eventdata)
        if frame == size(movie,3)
            frame = 1;
        end
        if strcmp(get(h_all.timer, 'Running'), 'off')
            start(h_all.timer);
        else
            stop(h_all.timer);
        end
    end

% Stop playing, adjust contrast, continue
    function contrastCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        % We clip the visible display range before the contrast dialog
        % to prevent a warning dialog about display range outside data range
        axes(h_all.axes);
        currImg = movie(:,:,frame);
        zImg = [min(currImg(:)), max(currImg(:))];
        zl = [max(zImg(1),zl(1)), min(zImg(2),zl(2))];
        caxis(zl);
        
        % Show contrast dialog, update color axis
        him = imcontrast;
        uiwait(him);
        zl = caxis;
        
        if isTimerOn
            start(h_all.timer);
        end
    end

% Stop playing, set contrast to match image min/max values, continue
    function autocontrastCallback(hObj, eventdata)
        axes(h_all.axes);
        xl = xlim; % update the axis limits in case the user zoomed
        yl = ylim;
        
        %         currImg = movie(:,:,frame); % Take whole frame for autocontrast
        % Take visible image cutout for autocontrast
        visibleXRange = max(1,floor(xl(1))):min(size(movie,2),ceil(xl(2)));
        visibleYRange = max(1,floor(yl(1))):min(size(movie,1),ceil(yl(2)));
        currImg = movie(visibleYRange,visibleXRange,frame);
        
        % Adjust contrast to match min/max intensity
        zl = [min(currImg(:)), max(currImg(:))];
        caxis(zl);
    end

% Plots the distribution selected by popup_distribution
    function distributionCallback(hObj, eventdata)
        selected_parameter = get(h_all.popup_distribution,'Value');
        choices = get(h_all.popup_distribution,'String');
        selected_string = choices{selected_parameter};
        
        % distribution to plot, must be synchronized with Popup-Menu (h_all.popup_distribution)
        dataRange = getNum(h_all.edit_distributionRange); % Range of data to histogram
        figure;
        ylabel('frequency');
        switch mode
            case 'candidate'
                if(isempty(allFramesCandidateData)) % Concatenate all frames
                    allFramesCandidateData = vertcat(candidateData{:});
                end
                rangedHist(allFramesCandidateData(:,selected_parameter), getNum(h_all.edit_distributionBins),dataRange);
                xlabel(fitParams{selected_parameter});
            case 'fitting'
                if(isempty(allFramesFitData)) % Concatenate all frames
                    allFramesFitData = vertcat(fitData{:});
                end
                rangedHist(allFramesFitData(:,selected_parameter), getNum(h_all.edit_distributionBins),dataRange);
                xlabel(fitParams{selected_parameter});
            case 'tracking'
                rangedHist(trackingData(:,selected_parameter), getNum(h_all.edit_distributionBins),dataRange);
                xlabel(trackingParams{selected_parameter});
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end
    end

% Switch black-white and hot display mode
    function bwCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        use_bw = ~use_bw;
        
        drawColors(n_colors); % Recompute colors
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay()
        end
        
    end

% Recompute the colors based on the current background
    function drawColors(num_colors)
        if use_bw
            bg = {'k'}; % background color
        else
            bg = {'r'};
        end
        
        marker_color = distinguishable_colors(1, bg);
        if dothandle_cand ~= -1
            set(dothandle_cand,'Color',marker_color);
        end
        
        marker_color = distinguishable_colors(1, bg);
        if dothandle_fit ~= -1
            set(dothandle_fit,'Color',marker_color);
            
        end
        
        track_colors = repmat( distinguishable_colors(num_colors, bg), ceil(n_tracks/num_colors) ,1);
        track_colors = track_colors(1:n_tracks,:);
        
        % Set the line color for each track
        for iH = 1:length(linehandles)
            if(linehandles(iH) ~= -1)
                set(linehandles(iH),'Color',track_colors(iH,:));
            end
        end
    end



% Update the movie FPS
    function fpsCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        FPS = str2num(get(h_all.edit_FPS, 'String'));
        if isempty(FPS) || FPS<=0
            FPS = 30;
        end
        
        % Timer is limited to 1ms precision
        if FPS>1000
            FPS = 1000;
            warning('Max FPS is 1000 due to timer precision');
        end
        timePerFrame = round(1/FPS*1000)/1000; % limit to millisecond precision
        set(h_all.timer,'Period', timePerFrame);
        set(h_all.edit_FPS,'String',sprintf('%.1f',FPS));
        
        if isTimerOn
            start(h_all.timer);
        end
    end

% This is called after letting the slider go. We update the frame display
% once more, otherwise synchronisation issues can occur.
    function sliderCallback(hObj, eventdata)
        updateFrameDisplay();
        elapsed_time = 0;
    end

% This is called continously when dragging the slider
    function updateSlider(hObj,eventdata)
        % Round slider value, sync text
        frame = round(get(h_all.slider,'Value'));
        set(h_all.slider,'Value',round(frame));
        updateTopText();
        
        % Stop timer
        if strcmp(get(h_all.timer, 'Running'), 'on')
            stop(h_all.timer);
        end
        
        % Skip frame if computer is too slow
        if elapsed_time >timePerFrame
            elapsed_time = elapsed_time - timePerFrame;
            return;
        end
        
        tic_start = tic;
        updateFrameDisplay();
        elapsed_time = elapsed_time + toc(tic_start)- timePerFrame;
    end

% Sets top text according to the current frame
    function updateTopText()
        set(h_all.toptext,'String',[sprintf('Current frame: %i/%i',frame,size(movie,3))])
    end

% Parse input variables
    function parse_inputs_and_setup(num_argin)
        % Is candidate data available?
        if ~isempty(candidateData)
            candidateParams = checkParameterDescription(candidateData, candidateParams);
            mode = 'candidate';
        else
            set(h_all.button_candidateMode,'Enable','off');
        end
        
        % Is fitting data available?
        if ~isempty(fitData)
            fitParams = checkParameterDescription(fitData, fitParams);
            mode = 'fitting';
        else
            set(h_all.button_fittingMode,'Enable','off');
        end
        
        % Is tracking data available?
        if ~isempty(trackingData)
            trackingParams = checkParameterDescription(trackingData, trackingParams);
            mode = 'tracking';
        else
            set(h_all.button_trackingMode,'Enable','off');
        end
        
        % input parsing
        if isempty(FPS)
            FPS = 30;
        end
        
        if isempty(is_blocking)
            is_blocking = false;
        end
        
    end

% Check if the number of parameters in the data and their description match.
% If there are not enough descriptions, pad with '<Unknown>', if there are
% too much, truncate.
    function description = checkParameterDescription(data, description)
        % Find number of parameters in fitData
        if iscell(data)
            nrParams = 0;
            for iFrame = 1:size(data,1)
                if ~isempty(data{iFrame})
                    nrParams = size(data{iFrame},2);
                    break;
                end
            end
        else
            nrParams = size(data,2);
        end
        
        % Make parameter description the same size as the number of params
        if isempty(description)
            description = repmat({'<Unknown>'}, nrParams,1);
        else
            if numel(description) > nrParams
                description = description(1:nrParams);
            elseif numel(description) < nrParams
                tmp = description;
                description = cell(nrParams,1);
                description(:) = {'<Unknown>'};
                description(1:numel(tmp)) = tmp(1:numel(tmp));
            end
        end
        
    end

% Cleanup function. This is neccessary to delete the timer!
    function onAppClose(hObj, event)
        if strcmp(get(h_all.timer, 'Running'), 'on')
            stop(h_all.timer);
        end
        delete(h_all.timer);
        delete(h_main);
    end


%%  Tracking only functions
    function lifetimeCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        traj_lifetime = round(str2num(get(h_all.edit_lifetime,'String')));
        if traj_lifetime<=0 || isempty(traj_lifetime)
            traj_lifetime = 0;
        end
        set(h_all.edit_lifetime,'String',sprintf('%i%',traj_lifetime));
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end

% Update the displayed length of trajectories
    function dispLengthCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        traj_displayLength = round(str2num(get(h_all.edit_trajDisplayLength,'String')));
        if traj_displayLength<=0 || isempty(traj_displayLength)
            traj_displayLength = 0;
        end
        set(h_all.edit_trajDisplayLength,'String',sprintf('%i%',traj_displayLength));
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end

% The user entered a different number of colors -> update color pool
    function colorCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        n_colors = round(str2num(get(h_all.edit_colors,'String')));
        if n_colors<=0 || isempty(n_colors)
            n_colors = 1;
        end
        set(h_all.edit_colors,'String',sprintf('%i%',n_colors));
        
        drawColors(n_colors);
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end

% Creates a line handle for every track that was not plotted before
    function setUnitinializedTrackHandles()
        hold on;
        for iTr = 1:n_tracks
            if(linehandles(iTr)==-1)
                linehandles(iTr) = plot(cell_traj{iTr}(1, 2), cell_traj{iTr}(1, 3), '.--','Color',track_colors(iTr,:));
                set(linehandles(iTr),'xdata',[],'ydata',[]);
            end
        end
        hold off;
    end
end

% Function that cuts data from upper and lower tails of the distribution
% keeping at least 'percentOfData' percent of all values.
function [freq, centers] = rangedHist(data, nbins, percentOfData)
if(percentOfData<100)
    % Limits for the cumulative density function (which goes from 0 to 1)
    lower_limit = (1-percentOfData/100)/2; % below this we throw away
    upper_limit = 1-(1-percentOfData/100)/2; % above this we throw away
    
    % Increase histogram resolution until highest bin has at most 10% of the data
    % This is neccessary to get a sufficiently good cumulative distribution function
    currBins = nbins;
    [freq,centers] = hist(data,currBins);
    maxRelativeBinCount = max(freq)/numel(data);
    % If there is only one unique value (e.g. PSF size), the method would
    % fail, so we set the number of bins very small instead
    if maxRelativeBinCount<1
        while(maxRelativeBinCount>0.1)
            currBins = currBins*2;
            [freq,centers] = hist(data,currBins);
            maxRelativeBinCount = max(freq)/numel(data);
        end
    else
        nbins=3;
    end
    
    % Integrate to get cumulative density function
    fraction = cumsum(freq)/sum(freq);
    % Find lower cutoff
    minIdx = find(fraction>lower_limit,1,'first')-1;
    if minIdx == 0
        minval = min(data);
    else
        minval = centers(minIdx);
    end
    % Find upper cutoff
    maxIdx = find(fraction>upper_limit,1,'first');
    if maxIdx == length(data)
        maxval = max(data);
    else
        maxval = centers(maxIdx);
    end
    
    data = data(data>=minval & data<=maxval);
end

freq = [];
centers = [];
switch(nargout)
    case 0
        hist(data,nbins);
    case 1
        freq = hist(data,nbins);
    case 2
        [freq,centers] = hist(data,nbins);
end

end
