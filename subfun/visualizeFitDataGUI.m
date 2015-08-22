function [h_main, run_again] = visualizeFitDataGUI(movie, fitData, FPS, use_bw, show_RunAgain)
% USAGE: visualizeFitDataGUI(movie, trajectoryData)
% [ Full USAGE: visualizeFitDataGUI(movie, fitData, FPS, use_bw, show_RunAgain) ]
%
% Visualizer for fit results (fitData).
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
%   show_RunAgain: Display run again dialog after closing. Used by tracker
%                  preview mode.
%
%  Inputs (except movie, fitData) can be left empty [] for default values.
%
% Output:
%   h_main - Handle to the GUI figure
%   run_again - Used for preview mode. Returns if the user selected to
%               run to software again in the onAppClose dialog
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2015
%

% Parse given inputs. For clarity we outsource this in a function.
if nargin<2 || isempty(movie) || isempty(fitData)
    fprintf(' Need input (movie,fitData)!\n');
    return
end
parse_inputs(nargin);
run_again = false;


% -- Preparing the GUI --
h_main = openfig('visualizeFitDataGUI_Layout.fig');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % Executed on closing for cleanup
set(h_main,'Toolbar','figure');   % Add toolbar needed for zooming
set(h_main, 'DoubleBuffer', 'on') % Helps against flickering

h_all = guihandles(h_main); % Get handles of all GUI objects
axes(h_all.axes); % Select axis for drawing plots

% Text on top
set(h_all.toptext,'String',sprintf('frame = 1/%i',size(movie,3)));

% Buttons
set(h_all.but_play,'Callback',@playCallback);
set(h_all.but_contrast,'Callback',@contrastCallback);

% Slider
set(h_all.slider,'Value',1, 'Min',1,'Max',size(movie,3),'SliderStep',[1/size(movie,3) 1/size(movie,3)],'Callback', @sliderCallback);
hLstn = handle.listener(h_all.slider,'ActionEvent',@updateSlider); % Add event listener for continous update of the shown slider value

% Edit fields
set(h_all.edit_FPS,'String',sprintf('%i',FPS), 'Callback', @fpsCallback);
set(h_all.edit_ampThresh, 'Callback', @callback_FloatEdit_Plus_Update);
setNum(h_all.edit_ampThresh, 0);
set(h_all.edit_snrThresh, 'Callback', @callback_FloatEdit_Plus_Update);
setNum(h_all.edit_snrThresh, 0);

% Checkbox
set(h_all.cb_bw, 'Value', use_bw, 'Callback',@bwCallback);

% Timer -> this controls playing the movie
h_all.timer = timer(...
    'ExecutionMode', 'fixedDelay', ...    % Run timer repeatedly
    'Period', round(1/FPS*1000)/1000, ... % Initial period is 1 sec. Limited to millisecond precision
    'TimerFcn', @onTimerUpdate, ...
    'StartFcn', @onTimerStart, ...
    'StopFcn',  @onTimerStop); % Specify callback



% -- Preparation for plotting --
% Convert data into format better for plotting
% (for Matlabs column major memory layout)
fitData_xCoords = squeeze(fitData(1,:,:));
fitData_yCoords = squeeze(fitData(2,:,:));
fitData_Amplitude = squeeze(fitData(3,:,:));
fitData_Background = squeeze(fitData(4,:,:));
fitData_SNR = fitData_Amplitude./fitData_Background;

% Draw the marker color depending on background color
marker_color = [];
drawColors();


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

% Variables for playback
timePerFrame = round(1/FPS*1000)/1000; % limit to millisecond precision
elapsed_time = 0;
frame = 1;

% In case the RunAgain dialog should be displayed, we stop scripts/functions
% calling the GUI until the figure is closed
if(show_RunAgain)
    uiwait(h_main);
end


% --- Nested Functions ---

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
        
        updateFrameDisplay();
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
        xl = xlim;
        yl = ylim;
        
        plotFrame(frame);
        
        % save the axis limits in case the user zoomed
        xlim(xl);
        ylim(yl);
        caxis(zl);
        
        drawnow; % Important! Or Matlab will skip drawing for high FPS
    end

    % Plots the frame with the input index
    function plotFrame(iF)
        imagesc(movie(:,:,iF)); axis image; colormap gray;
        if use_bw
            colormap gray;
        else
            colormap hot;
        end
        %     title(sprintf('Frame %i/%i',iF,size(movie,3)));
        
        % By default all particels are drawn
        ampMask = ones(size(fitData_Amplitude,1),1);
        snrMask = ones(size(fitData_Amplitude,1),1);
        
        % If thresholds are specified we filter the list of particles
        ampThresh = getNum(h_all.edit_ampThresh);
        snrThresh = getNum(h_all.edit_snrThresh);
        if( ~isempty(ampThresh) && ~(ampThresh==0) )
            ampMask = fitData_Amplitude(:,iF)>ampThresh;
        end
        if( ~isempty(ampThresh) && ~(snrThresh==0) )
            snrMask = fitData_SNR(:,iF)>snrThresh;
        end
        toPlotMask = ampMask & snrMask; % Particles with high enough amplitude AND signal-to-noise

        % Draw all particles above the given tresholds
        hold on;            
        plot(fitData_xCoords(toPlotMask,iF), fitData_yCoords(toPlotMask, iF), 'o','Color',marker_color);
        hold off;        
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

    % Switch black-white and hot display mode
    function bwCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        use_bw = ~use_bw;
        
        drawColors(); % Recompute colors
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay()
        end
        
    end

    % Recompute the colors based on the current background
    function drawColors()
        if use_bw
            bg = {'k'}; % background color
        else
            bg = {'r'};
        end
        
        marker_color = distinguishable_colors(1, bg);
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
        set(h_all.toptext,'String',[sprintf('frame =  %i/%i',frame,size(movie,3))])
    end

% Parse input variables
    function parse_inputs(num_argin)
        % input parsing
        if num_argin <3 || isempty(FPS)
            FPS = 30;
        end
        
        if num_argin <4 || isempty(use_bw)
            use_bw = false;
        end
        
        if num_argin < 5 || isempty(show_RunAgain)
            show_RunAgain = false;
        end
        
    end

% Cleanup function. This is neccessary to delete the timer!
    function onAppClose(hObj, event)
        if strcmp(get(h_all.timer, 'Running'), 'on')
            stop(h_all.timer);
        end
        delete(h_all.timer);
        
        % Dialog used in DEMO mode for getting return values.
        if show_RunAgain
            d = dialog('Position',[300 300 220 100],'Name','Run again?','WindowStyle','normal');
            
            txt = uicontrol('Parent',d,...
                'Style','text',...
                'Position',[5 40 210 40],...
                'String',sprintf('Run again to adjust settings?'));
            
            btn_yes = uicontrol('Parent',d,...
                'Position',[30 10 70 25],...
                'String','Yes',...
                'Callback',@buttonPress);
            
            btn_no = uicontrol('Parent',d,...
                'Position',[120 10 70 25],...
                'String','No',...
                'Callback',@buttonPress);
            run_again = false;
            uiwait(d);
        end
        delete(h_main);
        
        function buttonPress(hObj, event)
            run_again = strcmp(get(hObj,'String'),'Yes');
            delete(d);
        end
    end

end