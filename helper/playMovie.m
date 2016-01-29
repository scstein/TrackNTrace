function [h_main] = playMovie(movie, FPS, use_bw)
% USAGE: playMovie(movie, FPS, use_bw)
%
% Player for movies saved as 3D matrices.
%
% Input:
%   movie: 3D matrix (rows,cols,frames) movie
% Inputs also adjustable by GUI:
%   FPS: frames per second to play movie with | default: 30
%   use_bw: black/white image, otherwise colormap hot | default false
% Output:
%   h_main - Handle to the GUI figure
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2015
%

if nargin<1
    fprintf(' Need input movie!\n');
    return
end
num_argin = nargin;
parse_inputs();

% -- Preparing the GUI --
h_main = openfig('playMovie.fig');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % For cleanup
set(h_main,'Toolbar','figure'); % Toolbar needed for zooming
set(h_main, 'DoubleBuffer', 'on') % Helps against flickering
movegui(h_main,'center');

h_all = guihandles(h_main);
axes(h_all.axes);

% Text on top
set(h_all.toptext,'String',sprintf('frame = 1/%i',size(movie,3)));

% Buttons
set(h_all.but_play,'Callback',@playCallback);
set(h_all.but_contrast,'Callback',@contrastCallback);
set(h_all.but_autocontrast,'Callback',@autocontrastCallback);

% Slider
set(h_all.slider,'Value',1, 'Min',1,'Max',size(movie,3),'SliderStep',[1/size(movie,3) 1/size(movie,3)],'Callback', @sliderCallback);
hLstn = addlistener(h_all.slider,'ContinuousValueChange',@updateSlider); % For continous update of the shown slider value
set(h_all.slider,'ButtonDownFcn',@onSliderClick);

% Edit fields
set(h_all.edit_FPS,'String',sprintf('%i',FPS), 'Callback', @fpsCallback);

% Checkbox
set(h_all.cb_bw, 'Value', use_bw, 'Callback',@bwCallback);

% Timer -> this controls playing the movie
warning('off','MATLAB:TIMER:RATEPRECISION');
h_all.timer = timer(...
    'ExecutionMode', 'fixedDelay', ...   % Run timer repeatedly
    'Period', 1/FPS, ...                % Initial period is 1 sec.
    'TimerFcn', @timed_update, ...
    'StartFcn', @onTimerStart, ...
    'StopFcn',  @onTimerStop); % Specify callback
warning('on','MATLAB:TIMER:RATEPRECISION');


% Store handle to the image object (plotting is faster this way)
% Handles are set on first use (mostly in plotFrame)
imagehandle = -1;

% Plot first frame
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

% --- Nested Functions (mostly callbacks) ---

    % The main function of the application. This plays the movie if the
    % timer is running
    function timed_update(timer, event)
        frame = frame+1;
        if(frame >= size(movie,3))
            frame = size(movie,3);
            updateTopText();
            
            stop(h_all.timer);
        end
        set(h_all.slider,'Value',frame);
        updateTopText()
        
        % Skip frame if computer is too slow
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
        
        % Adjust contrast continously if shift key is pressed
        modifiers = get(gcf,'currentModifier');
        shiftIsPressed = ismember('shift',modifiers);
        if(shiftIsPressed)
           autocontrastCallback([],[]); 
        end
        
        drawnow;
    end

    % Plots the frame with the input index
    function plotFrame(iF)        
        if imagehandle == -1
            imagesc(movie(:,:,iF)); axis image; colormap gray;
        else
            set(imagehandle,'CData',movie(:,:,iF));
        end
        if use_bw
            colormap gray;
        else
            colormap hot;
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
        
        %This is done to prevent a warning dialog about display range outside data range
        axes(h_all.axes);
        currImg = movie(:,:,frame);
        zImg = [min(currImg(:)), max(currImg(:))];
        zl = [max(zImg(1),zl(1)), min(zImg(2),zl(2))];
        caxis(zl);
        
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

    % Switch black-white and hot display mode
    function bwCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        use_bw = ~use_bw;
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay()
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


    % This is called after letting the slider go
    function sliderCallback(hObj, eventdata)
        %
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
    function parse_inputs()
        % input parsing
        if num_argin <2 || isempty(FPS)
            FPS = 30;
        end
        
        if num_argin <3 || isempty(use_bw)
            use_bw = false;
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

end


