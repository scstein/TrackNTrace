% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
%     Copyright (C) 2016  Simon Christoph Stein, scstein@phys.uni-goettingen.de
% 	  Copyright (C) 2020, Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de
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
function [h_main, movie] = TNTvisualizer(movieOrTNTfile, candidateDataOrTNTfile, candidateParams, refinementData, refinementParams, trackingData, trackingParams, postprocData, postprocParams, FPS, is_blocking, firstFrame)
% Visualizer for movies and TrackNTrace results.
% 
% COMMON USAGEs: 
%    (x) TNTvisualizer();
%        Opens a dialog to choose a movie or TNT file from a completed TNT
%        run. When choosing a TNT file, the corresponding movie is loaded
%        automatically along all results from the evaluation.
%    (x) [~,movie] = TNTvisualizer();
%        Same as above, but the loaded movie is returned.
%    (x) [~,movie] = TNTvisualizer(PathToTNTfile);
%        Load TNT file. This loads the corresponding movie and all results from the evaluation.
%    (x) TNTvisualizer(movie);
%        Visualize movie already loaded in MATLAB.
%    (x) [~,movie] = TNTvisualizer(PathToMovie);
%        Same as above, but movie is loaded from disk first.
%    (x) TNTvisualizer(movie, PathToTNTfileOrStruct);
%        Visualize already loaded movie along with evaluated data from
%        specified TNT file. Make sure the movie fits to the evaluated data.
%    (x) [~,movie] = TNTvisualizer(PathToMovie, PathToTNTfileOrStruct); 
%        Same as above, but movie is loaded from disk first.

% [ Full USAGE: [GUIhandle, movie] = TNTvisualizer(movieOrTNTfile, candidateDataOrTNTfile, candidateParams, refinementData, refinementParams, trackingData, trackingParams, FPS, is_blocking, firstFrame) ]
%
% Input:
%   movieOrTNTfile: 
%     EITHER: 3D matrix (rows,cols,frames) with intensity values.
%     OR: FLIM movie: {intensity (rows,cols,frames), lifetime (rows,cols,frames)}
%     OR: Filepath to TIF movie -> Loads movie
%     OR: Filepath to TNT file. -> Loads movie and evaluation
%     OR: Struct like TNT file. -> Loads movie and evaluation
%   candidateDataOrTNTfile: 
%     EITHER: 1D cell array, where cell{F} is a 2D matrix with dimensions 
%       N(F)x(2+PC) saving the data of frame F. N(F) is the number
%       of candidates detected in frame F and the columns are 
%       'x','y' plus PC arbitrary additional parameters.
%     OR: Filepath to TNT file. -> Loads evaluated data 
%         (movie must have been given as the first parameter).
%     OR: Struct like TNT file.
%   candidateParams:
%     1D cell array with (2+PC) strings containing the names of the parameters
%     (columns of each cell) of candidateData (above).
%   refinementData: 
%     1D cell array, where cell{F} is a 2D matrix with dimensions
%     N(F)x(3+PF) saving the data of frame F. N(F) is the number
%     of fitted emitters detected in frame F and the columns are 
%     'x','y','z' position plus PF arbitrary additional parameters.
%   refinementParams:
%     1D cell array with (3+PF) strings containing the names of the parameters
%     (columns of each cell) of refinementData (above).
%   trackingData: 
%     2D matrix with dimensions Nx(5+PT). N is the number of all positions
%     tracked positions. The columns are 'TrackID','frame','x','y','z'
%     plus PT arbitrary additional parameters.
%   trackingParams:
%     1D cell array with (5+PT) strings containing the names of the parameters
%     (columns) of trackingData (above).
%   FPS:
%     Frames per second the visualizer plays set to use for playback on
%     startup.
%   is_blocking:
%     If true, the GUI blocks MATLAB execution until it closes.
%   firstFrame:
%     Index of first frame. This is used for displaying the true index
%     (with respect to the full movie) of a frame  if only a certain interval
%     of a movie was processed by TrackNTrace.
%     Example: If frame 201 to 250 from some movie were processed. Setting
%     firstFrame=201 shows 'Fr. 10/50 (210 in movie)' if the 10th processed
%     frame is shown.
%
%  Inputs can be left empty [] for default values. If the first or second
%  argument is a struct it is checked for the fields FPS, is_blocking,
%  and title.
%
% Output:
%   h_main:
%     Handle to the visualizer figure / GUI.
%   movie:
%     The visualized movie. Useful if movie was loaded by the visualizer.
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2016
% Author: Jan Christoph Thiele
% Date: 2020
%
% CT, 2019:
% - Added handeling of FLIM data
% - Added file export
% - Added tabbed GUI
% CT, 2020:
% - file preview
% - advanced filtering
% - reconstruction
% - drift correction
%

% Add all paths required to run TNT
addPathsVisualizer();

% Check MATLAB version
MATLABversion = strsplit(version,'.');
if(str2double(MATLABversion(1))>=8 && str2double(MATLABversion(2))>=6) % Use 'drawback nocallbacks' for MATLAB 2015b (8.6.x) and later
    MATLAB_2015b_or_newer = true; 
else
    MATLAB_2015b_or_newer = false; 
end

% Show filechooser dialog to choose movie / TNT file if started without input arguments.
importMode = false;
importFormats = {};
importPlugin = [];
if nargin==0
   [filename, path,filterindex] = uigetfile({'*.mat;*.tif','Visualizer files';'*.*','All files'},'Select movie or TNT file to visualize.');
    if isfloat(filename)
        % User clicked cancel
        importMode = true;
        movieOrTNTfile = '';
    else 
        [~,~,ext] = fileparts(filename);
        if filterindex == 1
            switch ext
                case '.tif'
                    movieOrTNTfile = [path,filename];
                case '.mat'
                    movieOrTNTfile = [path,filename];
                otherwise
                    error('Only movies or TNT files are supported for loading!');
            end
        else
            importMode = true;
            movieOrTNTfile = [path,filename];
        end
    end
end


% -- Preparing the GUI --
h_main = openfig(['subfun',filesep,'TNTvisualizer_Layout.fig'],'new','invisible');
set(h_main,'Renderer','painters');
% set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'handleVisibility','off');
set(h_main,'CloseRequestFcn',@onAppClose); % Executed on closing for cleanup
set(h_main,'Toolbar','figure');   % Add toolbar needed for zooming
set(h_main,'DoubleBuffer', 'on') % Helps against flickering
set(h_main,'FileName',''); % Set FileName to empty. This prevents that the user replaces TNTvisualizer_Layout.fig when pressing save.

h_all = guihandles(h_main); % Get handles of all GUI objects

% -- Turn off axes-specific toolbar and enable figure toolbar (for 2019a)
if ~verLessThan('matlab','9.5')
    addToolbarExplorationButtons(h_main) % Add the axes controls back to the figure toolbar 
    h_all.axes.Toolbar.Visible = 'off'; % Hide the integrated axes toolbar
    h_all.axes.Toolbar = []; % Remove the axes toolbar data
    try %#ok<TRYNC>
        tbar_removeList = {'DataManager.Linking','Standard.OpenInspector','Exploration.Rotate','Exploration.Brushing','Standard.EditPlot'};
        tbar = findall( h_main, 'tag', 'FigureToolBar');
        tbarbuttons = findall(tbar);
        delete(tbarbuttons(arrayfun(@(el)isgraphics(el)&&isprop(el,'tag')&&any(strcmp(el.Tag,tbar_removeList)),tbarbuttons)));
    end
end
% -- Datacursor --
dcm_obj = datacursormode(h_main);
% defaultDatatipFunction = get(dcm_obj,'UpdateFcn');
set(dcm_obj,'Enable','on');
% -- Scalebar --
try
    sb_obj = TNTscalebar(h_all.axes,'Visible','off');
    TNTscalebar.addToogle(h_main);
catch
    % TNTscalebar requires a newer matlab version
    sb_obj = [];
end


% -- Setup and input parsing --
% These variables are accessed/modified across different functions of the GUI
allFramesCandidateData = []; % Computed on demand for histogram plots
allFramesRefinementData = []; % Computed on demand for histogram plots

% use_bw = false; % Is visualization black/white or colored?
use_flim = false;                                                   % Plot FLIM image
cm_maps = {@cmap_heat,'hot','gray','jet','hsv',@cmap_isoluminant65,@cmap_isoluminant75,@cmap_rainbow_bgyrm};% List of available colormaps
cm_names = {'heat','hot','B/W','jet','hsv','iso65','iso75','Rainbow'};
cm_ind = 1;                                                         % Default colormap
cm_ind_inactive = 7;                                                % Default colormap for FLIM
cm_invert = false;                                                  % invert colormap
img_gamma = 1; % Gamma factor for plotting
mode = 'movie'; % 'movie' 'candidate'/'refinement'/'tracking'
traj_lifetime = 0; % Trajectories are kept for display for traj_lifetime after the particles last appearance
nr_track_colors = 20; % Size of colorpool to color tracks
traj_displayLength = inf; % Tracking data is displayed for the last traj_displayLength frames only (with respect to the current frame).

movie = [];   % Set by parse_inputs_and_setup depending on 'movieOrTNTfile'
movieLT = []; % Set by parse_inputs_and_setup depending on 'movieOrTNTfile'
candidateData = {}; % Set by parse_inputs_and_setup depending on 'candidateDataOrTNTfile'
id_tracks = []; % Unique ids for the tracks. If tracks are numbered continously from 1 to N, this is identical to the track index.
cell_traj = {}; % Each cell saves the data for one track. Is initialized by parse_inputs_and_setup.
n_tracks = 0;
titlename = ''; % (file)name to display in the window title
metadata = struct([]); % meta data for the file. Can contain fields: pixelsize (nm), framerate (1/s), filename
% TNTdata = struct([]); % Keep TNTdata
    
parse_inputs_and_setup(nargin);


if importMode
    set(h_all.panel_tabs,'Visible','off');
    set(h_all.panel_tabgroup,'Title','General options');
    set(h_all.button_preview,'Callback',@generatePreview);
    set(h_all.button_selectFile,'Callback',@callback_selectFile);
    set(h_all.button_SendToTNT,'Callback',@callback_sendToTNT);
    
    set(h_all.edit_firstFrame,'Callback',{@callback_IntEdit,1,inf});
    set(h_all.edit_lastFrame,'Callback',{@callback_IntEdit,1,inf});
    set(h_all.edit_binFrame,'Callback',{@callback_IntEdit,1,inf});
    set(h_all.cbx_useTimegate, 'Callback', @callback_toggleTimegate);
    set(h_all.edit_tgStart, 'Callback', {@callback_IntEdit,0,inf});
    set(h_all.edit_tgEnd, 'Callback', {@callback_IntEdit,0,inf});
    set(h_all.button_showTCSPC, 'Callback', @callback_showTCSPC);
    
else
    createTabs(h_all.panel_tabgroup,h_all.panel_tabs); % Create tabbed UI
    h_all = guihandles(h_main); % Get handles of all GUI objects. Refresh.
end
set(h_main,'Position',[1 1 0 1].*get(h_main,'Position')+[0 0 1 0].*(get(h_all.panel_tabgroup,'Position') *[2; 0; 1; 0])); % Adjust with
movegui(h_main,'center');

tracksVisibleInFrame = {}; % tracksVisibleInFrame{frame} is list of all track numbers visible in frame 'frame'
nr_VisibleTracksInFrame = []; % Number of visible tracks in each frame. 
maxNr_visibleTracksInFrame = 0; % Maximum number of tracks concurrently visible in one frame
compute_tracksVisibleInFrame(); % Sets tracksVisibleInFrame, nr_VisibleTracksInFrame and maxNr_visibleTracksInFrame.

% Set up vars for filters
candidateDataUnfiltered = candidateData;
refinementDataUnfiltered = refinementData;
trackingDataUnfiltered = trackingData;
postprocDataUnfiltered = postprocData;
candidateDataFilter = cellfun(@(d)true(size(d,1),1),candidateDataUnfiltered,'UniformOutput',false);
refinementDataFilter = cellfun(@(d)true(size(d,1),1),refinementDataUnfiltered,'UniformOutput',false);
trackingDataFilter = true(size(trackingDataUnfiltered,1),1);
postprocDataFilter = true(size(postprocDataUnfiltered,1),1);
% Drift correction
driftcalc  = zeros(size(movie,3),2); % Might be empty if movie was empty
driftcand  = driftcalc;
driftref   = driftcalc;
drifttrack = driftcalc;
driftpost  = driftcalc;
% ---------------------------

% Show the GUI!
% (We start invisible in case a movie is loaded during parse_inputs_and_setup)
set(h_main,'visible','on'); 
set(h_main,'CurrentAxes',h_all.axes);
% axes(h_all.axes); % Select axis for drawing plots. In some versions this also forces visablity.

%  -- Setup UI elements --
set(h_all.toptext,'ButtonDownFcn',@callback_Goto); % For text elements this unfortunately only fires for right clicks.
% Buttons
set(h_all.button_movieMode,'Callback', {@callback_changeMode,'movie'});
set(h_all.button_candidateMode,'Callback', {@callback_changeMode,'candidate'});
set(h_all.button_refinementMode,'Callback', {@callback_changeMode,'refinement'});
set(h_all.button_trackingMode,'Callback', {@callback_changeMode,'tracking'});
set(h_all.button_postprocMode,'Callback', {@callback_changeMode,'postproc'});

set(h_all.but_play,'Callback',@playCallback);
set(h_all.but_contrast,'Callback',@contrastCallback);
set(h_all.but_autocontrast,'Callback',@autocontrastCallback);
set(h_all.but_autocontrast,'TooltipString',sprintf('Set contrast automatically.\n The used algorithm can be selected in the popup menu to the right.\n Press  and hold "Shift" key during playback to adjust contrast automatically for each frame.'));
set(h_all.popup_autocontrast, 'TooltipString', sprintf('Algorithm used for autocontrast.\n   Spots: Emphasize highest 25%% intensity values.\n   Min/Max: Spans all values.\n   98%% range: Cuts the lower and upper 1%% of intensities. '));
set(h_all.popup_export, 'TooltipString', sprintf('These options export an RGB image of either the current frame or the whole stack with the current display options (colormap, contrast/limits, zoom).\n The images contain the markers for the localisations/tracks and for single frames also any set datatips.\n - Tiff (Stack)\n - Gif (Stack); Note: the GIF''s playback speed is set to the current FPS.\n - Tiff (Frame)\n - Png (Frame)\n\nThese options export the raw data (intensity / fast lifetime) with the native resolution and scale and without any overlay as single channel tiff.\n - Tiff (Stack, int)\n - Tiff (Stack, tau)\n'));
set(h_all.but_distribution,'Callback',@distributionCallback);
set(h_all.but_2d_distribution,'Callback',@distribution2DCallback);
set(h_all.but_weighted_distribution,'Callback',@distributionWeightedCallback);

% Slider
% Enableing/disabling is done in loadMovie()
hLstn = addlistener(h_all.slider,'ContinuousValueChange',@updateSlider); %#ok<NASGU> % Add event listener for continous update of the shown slider value

% Edit fields
set(h_all.edit_FPS,'String',sprintf('%i',FPS), 'Callback', @fpsCallback);
set(h_all.edit_distributionBins, 'Callback', {@callback_IntEdit,1,inf});
setNum(h_all.edit_distributionBins, 50, true);
set(h_all.edit_distributionRange,'Callback',{@callback_FloatEdit,0,100});
setNum(h_all.edit_distributionRange,100);
set(h_all.cb_distribution_allFrames,'Value',true);
set(h_all.edit_gamma,'Callback',@gammaCallback);
setNum(h_all.edit_gamma,img_gamma);

% Checkbox
% set(h_all.cb_bw, 'Value', use_bw, 'Callback',@bwCallback);
set(h_all.cb_invert, 'Value', cm_invert, 'Callback',@setColormap);
set(h_all.cb_flim, 'Value', use_flim, 'Callback',@flimCallback);

% Popupmenu
set(h_all.button_export, 'Callback', @exportCallback);
set(h_all.popup_distribution, 'String', 'No data');
set(h_all.popup_distribution_second, 'String', 'No data');
set(h_all.popup_filterParam, 'String', 'No data');
set(h_all.popup_reconstruct_mean, 'String', 'No data');

createColormapPopup(h_all.popup_colormap, cm_maps, cm_names);
set(h_all.popup_colormap, 'Callback', @setColormap);
set(h_all.popup_colormap, 'Value', cm_ind);

% Filtering panel
set(h_all.but_applyFilter, 'Callback', @applyFilter);
set(h_all.but_filterInsert, 'Callback', @callback_filter_insertParam);
set(h_all.but_showList, 'Callback', @showList);
set(h_all.but_exportList, 'Callback', @exportList);
set(h_all.but_resetFilter, 'Callback', @resetFilter);

% Reconstruction panel
set(h_all.but_reconstruct,'Callback',@callback_reconstruct);
set(h_all.but_reconstruct_mean,'Callback',{@callback_reconstruct,true});
set(h_all.edit_reconstruct_res,'Callback',{@callback_FloatEdit,1,1e3});
set(h_all.edit_reconstruct_pixelsize,'Callback',{@callback_FloatEdit,1,1e4});
set(h_all.edit_reconstruct_locprec,'Callback',{@callback_FloatEdit,0,1e4,'includeNaN'}); % Allow nan
if isfield(metadata,'pixelsize')&&~isempty(metadata.pixelsize)&&~isnan(metadata.pixelsize)
    if metadata.pixelsize < 1
        % Guess that it's in um
        set(h_all.edit_reconstruct_pixelsize,'Value',metadata.pixelsize*1e3);
    else
        % Guess that it's in nm
        set(h_all.edit_reconstruct_pixelsize,'Value',metadata.pixelsize);
    end
end

% Drift panel
set(h_all.but_drift_calc,'Callback',@calcDrift);
set(h_all.but_drift_apply,'Callback',{@callback_applyDrift,true});
set(h_all.but_drift_reset,'Callback',{@callback_applyDrift,false});
set(h_all.but_drift_show,'Callback',@callback_showDrift);
set(h_all.edit_drift_res,'Callback',{@callback_FloatEdit,1,inf});
set(h_all.edit_drift_rmax,'Callback',{@callback_FloatEdit,0.01,100});
set(h_all.edit_drift_seg,'Callback',{@callback_IntEdit,1,100000});
set(h_all.but_drift_export,'Callback',@callback_exportDrift);
set(h_all.but_drift_import,'Callback',@callback_importDrift);

%  -- Candidate UI elements --

%  -- Refinement UI elements --

%  -- Tracking UI elements --
set(h_all.edit_lifetime,'String',sprintf('%i',traj_lifetime), 'Callback', @callback_TrajLifetime);
set(h_all.edit_colors,'String',sprintf('%i',nr_track_colors), 'Callback', @callback_trackColors);
set(h_all.edit_trajDisplayLength,'String',sprintf('%i', traj_displayLength), 'Callback', @callback_dispLength);

% -- Timer -> this controls playing the movie --
h_all.timer = timer(...
    'ExecutionMode', 'fixedDelay', ...    % Run timer repeatedly
    'Period', round(1/FPS*1000)/1000, ... % Initial period is 1 sec. Limited to millisecond precision
    'TimerFcn', @onTimerUpdate, ...
    'StartFcn', @onTimerStart, ...
    'StopFcn',  @onTimerStop); % Specify callback

% Store handles to the plot objects (which is faster)
% Handles are set on first use (mostly in plotFrame)
linehandles = -1*ones(maxNr_visibleTracksInFrame,1); % Note: Size can increase dynamically if more lines are needed.
linehandleNr_to_TrackNr = -1*ones(maxNr_visibleTracksInFrame,1); % Stores which track the linehandle is assigned to
dothandle_fit = -1;
dothandle_cand = -1;
imagehandle = -1;

% Draw the marker color depending on background color
track_colors = [];
marker_color = [];
marker_fill_color = [];
drawColors(nr_track_colors);

% -- Variables for playback --
timePerFrame = round(1/FPS*1000)/1000; % limit to millisecond precision
elapsed_time = 0;
frame = 1;

%% initPlot
% Plot first frame to get limits right
% Set x,y,color limits
xl = [];
yl = [];
zl = [];
zl_alpha = [];

initPlot();

%%

% -- Change into the right mode (candidate/refinement/tracking) --
callback_changeMode();

% h_all.axes.HandleVisibility = 'callback'; % Prevents overplotting from outside

% For is_blocking==true we stop scripts/functions
% calling the GUI until the figure is closed
if(is_blocking)
    uiwait(h_main);
    drawnow; % makes figure disappear instantly (otherwise it looks like it is existing until script finishes)
end


% --- Nested Functions ---


% Change visualizer into the chosen mode 'movie' 'candidate','refinement','tracking'
% and display its relevant content.
% If "modus" input is given, the mode is set to "modus". Its
% implemented this way to use one callback for all buttons selecting the modes.
    function callback_changeMode(~,~,modus)
        % Note: nargin>2 is true if callback was invoked by a button
        % (and not from a direct call in this file)
        if nargin>2
            % Do nothing if the current modes button is pressed again
            if strcmp(modus,mode)
                return
            end
            mode = modus;
        end
        DEFAULT_COLOR = [0.941,0.941,0.941]; % Default color of buttons.
        SELECTED_COLOR = [0.65, 0.9, 0]; % Color of selected button.
        
        %Reset button colors
        set(h_all.button_movieMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_candidateMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_refinementMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_trackingMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_postprocMode,'BackgroundColor', DEFAULT_COLOR);
        
        
        % Set all mode specific panels invisible
%         set(h_all.panel_tracking,'Visible','off');
%         set(h_all.panel_histogram,'Visible','off');
        set(h_all.panel_tabgroup,'Visible','off');
             
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        % Delete all graphics objects, except the movie frame and invalidate all handles       
        resetGraphics();
        
        % Delete active datatip
        % WARNING this also deletes other hggroup objects associated with the figure which are also invisible and draggable.
        delete(findall(h_main,'Type','hggroup','HandleVisibility','off','Draggable','on'));
        
        % Mode specific changes (setting datatip function, highlight buttons etc.)
        switch mode
            case 'movie'
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_movieMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', 'No data');
                set(h_all.popup_distribution_second, 'String', 'No data');
                set(h_all.popup_filterParam, 'String', 'No data');
                set(h_all.popup_reconstruct_mean, 'String', 'No data');
            case 'candidate'
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_candidateMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', candidateParams);
                set(h_all.popup_distribution_second, 'String', candidateParams);
                set(h_all.popup_reconstruct_mean, 'String', candidateParams);
                set(h_all.popup_filterParam, 'String', matlab.lang.makeValidName(candidateParams));
%                 set(h_all.panel_histogram,'Visible','on');
                set(h_all.panel_tabgroup,'Visible','on');
            case 'refinement'
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_refinementMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', refinementParams);
                set(h_all.popup_distribution_second, 'String', refinementParams);
                set(h_all.popup_reconstruct_mean, 'String', refinementParams);
                set(h_all.popup_filterParam, 'String', matlab.lang.makeValidName(refinementParams));
%                 set(h_all.panel_histogram,'Visible','on');
                set(h_all.panel_tabgroup,'Visible','on');
            case 'tracking'
                initializeLinehandles(); %Initializes all needed line handles, if this is uncommented, linehandles are created on the fly
                
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_trackingMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', trackingParams);
                set(h_all.popup_distribution_second, 'String', trackingParams);
                set(h_all.popup_reconstruct_mean, 'String', trackingParams);
                set(h_all.popup_filterParam, 'String', matlab.lang.makeValidName(trackingParams));
%                 set(h_all.panel_histogram,'Visible','on');
%                 set(h_all.panel_tracking,'Visible','on');  
                set(h_all.panel_tabgroup,'Visible','on'); 
            case 'postproc'
                
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_postprocMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', postprocParams);
                set(h_all.popup_distribution_second, 'String', postprocParams);
                set(h_all.popup_reconstruct_mean, 'String', postprocParams);
                set(h_all.popup_filterParam, 'String', matlab.lang.makeValidName(postprocParams));
%                 set(h_all.panel_histogram,'Visible','on');
                set(h_all.panel_tabgroup,'Visible','on');
            otherwise
                error('Unkown mode ''%s''!', mode);
        end
        if importMode
            set(h_all.panel_tabgroup,'Visible','on');
        end
        % Select first histogram entry (There should always be a first entry!)
        % If not set, this runs into problems if selecting entry X in some
        % mode and there are less then X parameters when the mode is switched.
        set(h_all.popup_distribution, 'Value', 1);        
        set(h_all.popup_distribution_second, 'Value', 1);    
        set(h_all.popup_filterParam, 'Value', 1);        
        set(h_all.popup_reconstruct_mean, 'Value', max([1 find(strcmp(get(h_all.popup_reconstruct_mean, 'String'),'lt-tau'),1)]));
        % Resize GUI
        resizeGUIforMode();
        
        % Replot
        updateFrameDisplay();
        
        if isTimerOn
            start(h_all.timer);
        end
    end

    % Delete all graphics objects, except the movie frame and invalidate all handles       
    function resetGraphics()
        allHandles = get(h_all.axes,'Children');
        if imagehandle ~= -1 % Remove image handle from list
            allHandles(allHandles==imagehandle) = [];
        end
        delete(allHandles);
        dothandle_fit = -1;
        dothandle_cand = -1;
        linehandles = -1*ones(maxNr_visibleTracksInFrame,1);
    end

% Update size of GUI, show mode specific panels
    function resizeGUIforMode()
        units = 'characters';
        BOTTOM_SPACING = 0.5;
        
        % Get position of last element in GUI (mode dependent)
        if importMode
            if strcmpi(get(h_all.panel_importOptions,'Visible'),'on')
                lowestUI = h_all.panel_importOptions;
            else
                lowestUI = h_all.panel_tabgroup;
            end
        else
            switch mode
                case 'movie'
                    lowestUI = h_all.panel_player;
                case {'candidate','refinement','tracking','postproc'}
                    lowestUI = h_all.panel_tabgroup;
            end
        end
        set(lowestUI,'Units',units);
        pos = get(lowestUI,'Position');
        
        diff_height = pos(2)-BOTTOM_SPACING;
        
        if abs(diff_height)>0.1
            % Only set positions when they change.
            
            % Remove colorbar if exist, since the automatic resizeing
            % is disabled when the axis position is changed.
            hasCB = isappdata(h_all.axes,'LayoutPeers') && ~isempty(findobj(h_all.axes.Parent,'Tag','Colorbar'));
            if hasCB
                colorbar(h_all.axes,'off');
                drawnow; % To update the axes position
            end           
            
            set(h_main,'Units',units);
            win_pos = get(h_main,'Position');
            win_top_pos = win_pos(2)+win_pos(4); % Save top position
            % To resize the figure properly, we first need to move all objects
            % inside.. (Matlab ..)
            % Except the legend which is a child of the figure (not the axes)
            % but is nevertheless positioned relative to the axes.
            all_uiObjects = findobj(get(h_main,'Children'),'flat','-not',{'Tag','legend','-or','Tag','Colorbar'});
            
            for iObj = 1:numel(all_uiObjects)
                try %#ok<TRYNC>
                    set(all_uiObjects(iObj),'Units',units);
                    pos = get(all_uiObjects(iObj),'Position');
                    pos(2) = pos(2) - diff_height;
                    set(all_uiObjects(iObj),'Position',pos);
                end
            end
            
            win_pos(4) = win_pos(4)-diff_height; % Set new window height
            win_pos(2) = win_top_pos-win_pos(4); % Keeps the top position constant
            set(h_main,'Position', win_pos);
            
            % Reset units back to normalized, so figure resizes "properly" (cough..)
            set(h_main,'Units','normalized');
            for iObj = 1:numel(all_uiObjects)
                try %#ok<TRYNC>
                    set(all_uiObjects(iObj),'Units','normalized');
                end
            end
            % Restore colorbar
            if hasCB
                colorbar(h_all.axes);
            end
        else
            set(lowestUI,'Units','normalized');            
        end
    end

    % Custom function for datacursor which shows data relevant to the
    % current mode when clicking the currently plotted data.
    function txt = modeSpecificDatatipFunction(~,event_obj)
       
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        graphObjHandle = get(event_obj,'Target'); % The target object (line/image) of the cursor
        
        if(isgraphics(graphObjHandle,'image')) % Image is selected
            txt = {['X: ',num2str(pos(1))],...
                   ['Y: ',num2str(pos(2))],...
                   ['Intensity: ', num2str(movie(pos(2),pos(1),frame))]};
            if ~isempty(movieLT)
                txt = [txt,...
                    {['\tau_{fast}: ', num2str(movieLT(pos(2),pos(1),frame),3) ' ns']}];
            end
        else % Plotted position is selected
            I = get(event_obj, 'DataIndex');
            txt = {};
            switch mode
                case 'movie' % -> Image is selected handled above
                case 'candidate'
                    % Plot all parameters available for selected datapoint in the datacursor window
                    for iPar=1:numel(candidateParams)
                        txt = [txt, {[candidateParams{iPar},': ', num2str(candidateData{frame}(I,iPar))]}]; %#ok<AGROW>
                    end
                case 'refinement'
                    % Plot all parameters available for selected datapoint in the datacursor window
                    for iPar=1:numel(refinementParams)
                        txt = [txt, {[refinementParams{iPar},': ', num2str(refinementData{frame}(I,iPar))]}]; %#ok<AGROW>
                    end
                case 'tracking'
                    handleNr = linehandles==graphObjHandle; % Find lineobject for the selected point
                    TrackNr = linehandleNr_to_TrackNr(handleNr);
                    PointData = cell_traj{TrackNr}(I,:); % Data of the selected point
                    TrackID = sprintf('%i',id_tracks(TrackNr)); % Get track ID from its index (in case TracIDs go from 1 to N without missing numbers, TrackNr==TrackID)
                    % Plot all parameters available for selected datapoint in the datacursor window
                    txt = [txt, {['TrackID: ', TrackID]}];
                    for iPar=2:numel(trackingParams)
                        txt = [txt, {[trackingParams{iPar},': ', num2str(PointData(iPar-1))]}]; %#ok<AGROW>
                    end
                case 'postproc'
                    ind_frame = find(postprocData(:,2) == frame);
                    PointData = postprocData(ind_frame(I),:); % Data of the selected point
                    % Plot all parameters available for selected datapoint in the datacursor window
                    for iPar=1:numel(postprocParams)
                        txt = [txt, {[postprocParams{iPar},': ', num2str(PointData(iPar))]}]; %#ok<AGROW>
                    end
                otherwise
                    error('Unsupported mode ''%s'' for datatip function.',mode)
            end
        end
    end

% Get numeric value of edit field
    function value = getNum(hObj)
        value = str2double(get(hObj,'String'));
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

% % Callback for edit fields containing floats. Checks if a correct
% % number was entered and restricts it to the given bounds.
%     function callback_FloatEdit(hObj,~, minVal, maxVal)
%         if nargin<3 || isempty(minVal);
%             minVal=-inf;
%         end
%         if nargin<4 || isempty(maxVal);
%             maxVal=inf;
%         end
%         
%         % Check if a valid number was entered
%         value = str2double(get(hObj, 'String'));
%         if isempty(value)
%             set(hObj,'ForegroundColor','r');
%             set(hObj,'String','INVALID');
%             uicontrol(hObj);
%         else
%             value = max(minVal,value);
%             value = min(maxVal,value);
%             set(hObj,'ForegroundColor','k');
%             set(hObj,'String',sprintf('%.2f',value));
%         end
%     end
% 
% % Callback for edit fields containing integer values. Checks if a correct
% % number was entered and restricts it to the given bounds.
%     function callback_intEdit(hObj,~, minVal,maxVal)
%         if nargin<3 || isempty(minVal);
%             minVal=0;
%         end
%         if nargin<4 || isempty(maxVal);
%             maxVal=inf;
%         end
%         
%         value = round(str2double(get(hObj,'String')));
%         if isempty(value)
%             set(hObj,'ForegroundColor','r');
%             set(hObj,'String','INVALID');
%             uicontrol(hObj);
%         else
%             value = max(minVal,value);
%             value = min(maxVal,value);
%             set(hObj,'ForegroundColor','k');
%             set(hObj,'String',sprintf('%i',value));
%         end
%     end

% The main function of the application. This plays the movie if the timer is running
    function onTimerUpdate(~, ~)
        % Progress frame counter, clip at length of movie and stop at last frame.
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

    % Resets the measured elapsed time since last drawn frame when the timer starts
    function onTimerStart(~, ~)
        set(h_all.but_play,'String','Pause');
        elapsed_time = 0;
    end

    function onTimerStop(~, ~)
        set(h_all.but_play,'String','Play');
        updateFrameDisplay();
    end

% Sets the limits and plots the first frame
    function initPlot()
        frame = 1;
        elapsed_time = 0;
        
        if isempty(movie)
            set(h_all.text_noPreview,'Visible','on');
        else
            set(h_all.text_noPreview,'Visible','off');
            xl = [0.5,size(movie,2)+0.5];
            yl = [0.5,size(movie,1)+0.5];
            firstImg = movie(:,:,1).^img_gamma;
            zl = [min(firstImg(:)), max(firstImg(:))];
            if zl(1)==zl(2)
                zl = zl + [0 1];
            end
            if ~isempty(movieLT)
                firstImg = movieLT(:,:,1);
                zl_alpha = [min(firstImg(:)), max(firstImg(:))];% remeber zl when switching to FLIM
            end
            
            plotFrame(frame);
            
            xlim(h_all.axes,xl);
            ylim(h_all.axes,yl);
            caxis(h_all.axes,zl);
        end
        
        updateTopText();
        
        % Update scalebar unit
        if ~isempty(sb_obj) && ~isempty(metadata) && isstruct(metadata) 
            if isfield(metadata,'pixelsize') && ~isempty(metadata.pixelsize)
                sb_obj.Pixelsize = metadata.pixelsize;
            else
                sb_obj.Pixelsize = 1;
            end
            if isfield(metadata,'pixelsize_unit') && ~isempty(metadata.pixelsize_unit)
                sb_obj.Unit = [' ',metadata.pixelsize_unit];
            else
                sb_obj.Unit = '';
            end
        end
        % Prepare FLIM
        try
            set(h_all.axes,'Alphamap',linspace(0,1,256));% Increase resolution of alpha to 8 bit
        catch
            set(get(h_all.axes,'ColorSpace'),'Alphamap',linspace(0,1,256));
        end
    end

% Call to display the current frame as selected by the 'frame' variable
% Also this sets and saves the axis states (e.g. for zooming);
    function updateFrameDisplay()
        % Needed to minimize interference with other figures the user
        % brings into focus. It can be that the images are not plotted to
        % the GUI then but to the selected figure window
%         set(0,'CurrentFigure',h_main);
        
        % Delete active datatip
        % WARNING this also deletes other hggroup objects associated with the figure which are also invisible and draggable.
        delete(findall(h_main,'Type','hggroup','HandleVisibility','off','Draggable','on'));

        plotFrame(frame);
        
        % Adjust contrast continously if shift key is pressed
        modifiers = get(h_main,'currentModifier');
        shiftIsPressed = ismember('shift',modifiers);
        if(shiftIsPressed)
            autocontrastCallback([],[]);
        end
        
        % Important! Or Matlab will skip drawing entirely for high FPS
        if(MATLAB_2015b_or_newer)
            drawnow nocallbacks; 
        else
            drawnow expose update;
        end
    end

    % Plots contents of the frame with the input index iF
    function plotFrame(iF)
        % Draw frame iF of the movie 
        if imagehandle == -1
            imagehandle = imagesc(h_all.axes,movie(:,:,iF).^img_gamma); 
            axis(h_all.axes,'image');
            setColormap();
            return % setColormap() internally calls plotFrame(frame)
        end
        if use_flim
            set(imagehandle,'CData',movieLT(:,:,iF));
            set(imagehandle,'AlphaData',movie(:,:,iF).^img_gamma);
        else
            set(imagehandle,'CData',movie(:,:,iF).^img_gamma);
        end
        
        % Draw mode dependent data
        switch mode
            case 'movie'
                % Nothing additional to draw.
            case 'candidate'
                if(isempty(candidateData{iF}));
                    if dothandle_cand ~= -1 % Skip uninitialized handles (must be drawn once)
                        set(dothandle_cand,'xdata',[],'ydata',[]);
                    end
                    return % Jump empty frames
                end 
                
                % Plot markers of candidates
                hold(h_all.axes,'on');
                if dothandle_cand == -1 % Draw unitialized handles
                    dothandle_cand = plot(h_all.axes,candidateData{iF}(:,1), candidateData{iF}(:,2), 's','Color',marker_color,'MarkerSize',5,'MarkerFaceColor', marker_fill_color','DisplayName','Localization');
                else % For initialized handles set their data (MUCH faster than plot)
                    set(dothandle_cand,'xdata',candidateData{iF}(:,1),'ydata',candidateData{iF}(:,2));
                end
                hold(h_all.axes,'off');
                
                
            case 'refinement'
                if(isempty(refinementData{iF}));
                    if dothandle_fit ~= -1 % Skip uninitialized handles (must be drawn once)
                        set(dothandle_fit,'xdata',[],'ydata',[]);
                    end
                    return % Jump empty frames
                end 
                
                % Plot markers of fitted positions
                hold(h_all.axes,'on');
                if dothandle_fit == -1 % Draw unitialized handles
                    dothandle_fit = plot(h_all.axes,refinementData{iF}(:,1), refinementData{iF}(:,2), 'o','Color',marker_color,'MarkerSize',5,'MarkerFaceColor', marker_fill_color ,'Linewidth',1,'DisplayName','Localization');
                else % For initialized handles set their data (MUCH faster than plot)
                    set(dothandle_fit,'xdata',refinementData{iF}(:,1),'ydata',refinementData{iF}(:,2));
                end
                hold(h_all.axes,'off');
                
                
            case 'tracking'
                % Draw the tracks of currently visible particles
                hold(h_all.axes,'on');
                handleNr = 1;
                for iTr = tracksVisibleInFrame{iF}                    
                    % Plot trajectories a) only the last traj_displayLength positions AND  b) up to the current frame
                    mask_toPlot = ((cell_traj{iTr}(:,1)>iF-traj_displayLength) & cell_traj{iTr}(:,1)<=iF);
                    
                    % We use the next free linehandle. If there is no
                    % free handle left, we create a new one on the fly.
                    if (handleNr>numel(linehandles) || linehandles(handleNr) == -1)
                        linehandles(handleNr) = plot(h_all.axes,cell_traj{iTr}(mask_toPlot, 2), cell_traj{iTr}(mask_toPlot, 3), '.-','Color',track_colors(iTr,:),'Linewidth',1,'DisplayName',sprintf('Track %.0f',handleNr));
                        linehandleNr_to_TrackNr(handleNr) = iTr;
                        handleNr = handleNr+1;
                    else
                        set(linehandles(handleNr),'xdata',cell_traj{iTr}(mask_toPlot, 2),'ydata', cell_traj{iTr}(mask_toPlot, 3),'Color',track_colors(iTr,:));
                        linehandleNr_to_TrackNr(handleNr) = iTr;
                        handleNr=handleNr+1;
                    end
                end
                % Empty data of unused linehandles
                while(handleNr<= numel(linehandles))
                    if(linehandles(handleNr) ~= -1)                        
                        set(linehandles(handleNr),'xdata',[],'ydata',[]);
                        linehandleNr_to_TrackNr(handleNr) = -1;
                    end
                    handleNr=handleNr+1;
                end
                hold(h_all.axes,'off');
                
                
            case 'postproc'
                ind_frame = postprocData(:,2) == iF;
                if sum(ind_frame)==0
                    if dothandle_fit ~= -1 % Skip uninitialized handles (must be drawn once)
                        set(dothandle_fit,'xdata',[],'ydata',[]);
                    end
                    return % Jump empty frames
                end 
                
                % Plot markers of fitted positions
                hold(h_all.axes,'on');
                if dothandle_fit == -1 % Draw unitialized handles
                    dothandle_fit = plot(h_all.axes,postprocData(ind_frame,3), postprocData(ind_frame,4), 'o','Color',marker_color,'MarkerSize',5,'MarkerFaceColor', marker_fill_color ,'Linewidth',1,'DisplayName','Localization');
                else % For initialized handles set their data (MUCH faster than plot)
                    set(dothandle_fit,'xdata',postprocData(ind_frame,3),'ydata',postprocData(ind_frame,4));
                end
                hold(h_all.axes,'off');
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end
    end

    % Switch play/pause by button
    function playCallback(~, ~)
        if frame == size(movie,3)
            frame = 1;
        end
        if strcmp(get(h_all.timer, 'Running'), 'off')
            start(h_all.timer);
            
            % Disable datatip during playback (otherwise MATLAB tries to
            % put a datatip on the continously updated movie frame, which
            % can lead to errors).
            set(dcm_obj,'Enable','off'); 
        else
            stop(h_all.timer);
            
            % Enable datatip
            set(dcm_obj,'Enable','on'); 
        end
    end

    % Stop playing, adjust contrast manually, continue
    function contrastCallback(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        modifiers = get(h_main,'currentModifier');
        ctrlIsPressed = ismember('control',modifiers);
        % We clip the visible display range before the contrast dialog
        % to prevent a warning dialog about display range outside data range
%         axes(h_all.axes);

        if use_flim && ctrlIsPressed
            currImg = movie(:,:,frame).^img_gamma;
            zImg = [min(currImg(:)), max(currImg(:))];
            zl_alpha = double([max(zImg(1),zl_alpha(1)), min(zImg(2),zl_alpha(2))]);
            alim(h_all.axes,zl_alpha); 
            
            % Show contrast dialog, update color axis
            him = imcontrast_alpha(h_all.axes);
            % Disable adjust data button
            set(findobj(him,'String','Adjust Data'),'Visible','off');
            uiwait(him);
            zl_alpha = alim(h_all.axes);
            
        else
            if use_flim
                currImg = movieLT(:,:,frame);
            else
                currImg = movie(:,:,frame).^img_gamma;
            end
            zImg = [min(currImg(:)), max(currImg(:))];
            zl = double([max(zImg(1),zl(1)), min(zImg(2),zl(2))]);
            caxis(h_all.axes,zl);
            
            % Show contrast dialog, update color axis
            him = imcontrast(h_all.axes);
            % Disable adjust data button
            set(findobj(him,'String','Adjust Data'),'Visible','off');
            uiwait(him);
            zl = caxis(h_all.axes);
        end
        
        if isTimerOn
            start(h_all.timer);
        end
    end

    % Stop playing, set contrast to match image min/max values, continue
    function autocontrastCallback(~, ~)
%         axes(h_all.axes);
        xl = xlim(h_all.axes); % update the axis limits in case the user zoomed
        yl = ylim(h_all.axes);
        
        %         currImg = movie(:,:,frame); % Take whole frame for autocontrast
        % Take visible image cutout for autocontrast
        visibleXRange = max(1,floor(xl(1))):min(size(movie,2),ceil(xl(2)));
        visibleYRange = max(1,floor(yl(1))):min(size(movie,1),ceil(yl(2)));
        if use_flim
            currImg = movieLT(visibleYRange,visibleXRange,frame);
        else
            currImg = movie(visibleYRange,visibleXRange,frame).^img_gamma;            
        end
                
        selected_method =  get(h_all.popup_autocontrast,'Value');        

        switch selected_method
            case 1 % Focus on spots (upper quartil / 25% of data)
                if use_flim
                    currImgLT = currImg;
                    currImg = movie(visibleYRange,visibleXRange,frame).^img_gamma;
                end
                sortedIntensities = sort(currImg(:));                 
                minIdx = round(0.75*numel(sortedIntensities));
                if minIdx == 0; minIdx = 1; end
                minval = sortedIntensities( minIdx );
                maxval = sortedIntensities( end );
                if use_flim
                    % take bounds centre of centre 75% lifetime values of
                    % the brightest 25% pixels.
                    sortedIntensitiesLTs = sort(currImgLT(currImg>=minval & currImg<=maxval));                
                    minIdx = max(1,round(0.125*numel(sortedIntensitiesLTs)));
                    maxIdx = min(numel(sortedIntensitiesLTs),round(0.875*numel(sortedIntensitiesLTs)));
                    minval = sortedIntensitiesLTs(minIdx);
                    maxval = sortedIntensitiesLTs(maxIdx);
                end
            case 2 % Adjust contrast to match min/max intensity
                minval = min(currImg(:));
                maxval = max(currImg(:));
            case 3 % 98% range of data (cut upper / lower tails)
                sortedIntensities = sort(currImg(:)); 
                minIdx = round(0.01*numel(sortedIntensities));
                if minIdx == 0; minIdx = 1; end
                minval = sortedIntensities( minIdx );
                maxval = sortedIntensities( round(0.99*numel(sortedIntensities))); 
            case 4 % 50% range of data (cut upper / lower tails)
                sortedIntensities = sort(currImg(:)); 
                minIdx = round(0.25*numel(sortedIntensities));
                if minIdx == 0; minIdx = 1; end
                minval = sortedIntensities( minIdx );
                maxval = sortedIntensities( round(0.75*numel(sortedIntensities)));                        
            otherwise
                error('Unknown autocontrast mode!')
        end
        if minval<=maxval % this also excudes NaNs
            zl = double([minval, maxval]);
            caxis(h_all.axes,zl);
        end
    end

    % Plots the distribution/histogram of the parameter selected by popup_distribution
    % Selectable parameters vary depending on data available to the current mode
    function distributionCallback(~, ~)
        selected_parameter = get(h_all.popup_distribution,'Value');
        nBins = getNum(h_all.edit_distributionBins);
        distributionHelper(selected_parameter,[],nBins);
    end
    % Plots 2D histogram
    function distribution2DCallback(~, ~)
        selected_parameter = get(h_all.popup_distribution,'Value');
        selected_parameter2 = get(h_all.popup_distribution_second,'Value');
        nBins = getNum(h_all.edit_distributionBins);
        distributionHelper(selected_parameter,selected_parameter2,[nBins nBins]);% Give the number of bins in both directions to make it 2D
    end
    % Plots weighted histogram
    function distributionWeightedCallback(~, ~)
        selected_parameter = get(h_all.popup_distribution,'Value');
        selected_parameter2 = get(h_all.popup_distribution_second,'Value');
        nBins = getNum(h_all.edit_distributionBins);
        distributionHelper(selected_parameter,selected_parameter2,nBins);% Give the number of bins in both directions to make it 2D
    end

    function distributionHelper(selected_parameter, selected_parameter2,nBins)
        
        is2D = numel(nBins)==2;
        
        allFrames = get(h_all.cb_distribution_allFrames,'Value');
        avgTracks  = get(h_all.cb_distribution_avg_track,'Value'); 
        dataRange = getNum(h_all.edit_distributionRange); % Range of data to histogram
        ylabelText='frequency';
        switch mode
            case 'candidate'
                if ~allFrames
                    currCandidateData = candidateData{frame};
                    hdata1 = currCandidateData(:,selected_parameter);
                    hdata2 = currCandidateData(:,selected_parameter2);
%                     histobj = rangedHist(currCandidateData(:,selected_parameter),currCandidateData(:,selected_parameter2), nBins,dataRange);
                else
                    % Save concatenated data of all frames (done only once)
                    if(isempty(allFramesCandidateData))
                        allFramesCandidateData = vertcat(candidateData{:});
                    end
                    hdata1 = allFramesCandidateData(:,selected_parameter);
                    hdata2 = allFramesCandidateData(:,selected_parameter2);
%                     histobj = rangedHist(allFramesCandidateData(:,selected_parameter),allFramesCandidateData(:,selected_parameter2), nBins,dataRange);
                end
                xlabelText = candidateParams{selected_parameter};
                if is2D
                    ylabelText = candidateParams{selected_parameter2};
                end
%                 xlabel(candidateParams{selected_parameter},'Interpreter','none');
%                 if is2D
%                     ylabel(candidateParams{selected_parameter2},'Interpreter','none');
%                 else
%                     ylabel('frequency');
%                 end
                avgTracks = false;
            case 'refinement'
                if ~allFrames
                    currRefinementData = refinementData{frame};
                    hdata1 = currRefinementData(:,selected_parameter);
                    hdata2 = currRefinementData(:,selected_parameter2);
%                     histobj = rangedHist(currRefinementData(:,selected_parameter),currRefinementData(:,selected_parameter2), nBins,dataRange);                    
                else
                    % Save concatenated data of all frames (done only once)
                    if(isempty(allFramesRefinementData))
                        allFramesRefinementData = vertcat(refinementData{:});
                    end
                    hdata1 = allFramesRefinementData(:,selected_parameter);
                    hdata2 = allFramesRefinementData(:,selected_parameter2);
%                     histobj = rangedHist(allFramesRefinementData(:,selected_parameter),allFramesRefinementData(:,selected_parameter2), nBins,dataRange);
                end
                xlabelText = refinementParams{selected_parameter};
                if is2D
                    ylabelText = refinementParams{selected_parameter2};
                end
%                 xlabel(refinementParams{selected_parameter},'Interpreter','none');
%                 if is2D
%                     ylabel(refinementParams{selected_parameter2},'Interpreter','none');
%                 else
%                     ylabel('frequency');
%                 end
                avgTracks = false;
            case 'tracking'
                if ~allFrames
                    ind = trackingData(:,2)==frame;
                    hdata1 = trackingData(ind,selected_parameter);
                    hdata2 = trackingData(ind,selected_parameter2);
%                     histobj = rangedHist(trackingData(ind,selected_parameter),trackingData(ind,selected_parameter2), nBins,dataRange);                    
                elseif avgTracks
                    [~,~,ind] = unique(trackingData(:,1)); % The unique ensures that ind starts at one and has no gaps
                    hdata1 = accumarray(ind,trackingData(:,selected_parameter),[],@mean);
                    if isempty(selected_parameter2)
                        hdata2 = [];
                    else
                        hdata2 = accumarray(ind,trackingData(:,selected_parameter2),[],@mean);
                    end
                else
                    hdata1 = trackingData(:,selected_parameter);
                    hdata2 = trackingData(:,selected_parameter2);
%                     histobj = rangedHist(trackingData(:,selected_parameter),trackingData(:,selected_parameter2), nBins,dataRange);
                end
                xlabelText = trackingParams{selected_parameter};
                if is2D
                    ylabelText = trackingParams{selected_parameter2};
                end
%                 xlabel(trackingParams{selected_parameter},'Interpreter','none');
%                 if is2D
%                     ylabel(trackingParams{selected_parameter2},'Interpreter','none');
%                 else
%                     ylabel('frequency');
%                 end
            case 'postproc'
                if ~allFrames
                    ind = postprocData(:,2)==frame;
                    hdata1 = postprocData(ind,selected_parameter);
                    hdata2 = postprocData(ind,selected_parameter2);
%                     histobj = rangedHist(postprocData(ind,selected_parameter),postprocData(ind,selected_parameter2), nBins,dataRange);                    
                elseif avgTracks
                    [~,~,ind] = unique(postprocData(:,1)); % The unique ensures that ind starts at one and has no gaps
                    hdata1 = accumarray(ind,postprocData(:,selected_parameter),[],@mean);
                    if isempty(selected_parameter2)
                        hdata2 = [];
                    else
                        hdata2 = accumarray(ind,postprocData(:,selected_parameter2),[],@mean);
                    end
                else
                    hdata1 = postprocData(:,selected_parameter);
                    hdata2 = postprocData(:,selected_parameter2);
%                     histobj = rangedHist(postprocData(:,selected_parameter),postprocData(:,selected_parameter2), nBins,dataRange);
                end
                xlabelText = postprocParams{selected_parameter};
                if is2D
                    ylabelText = postprocParams{selected_parameter2};
                end
%                 xlabel(postprocParams{selected_parameter},'Interpreter','none');
%                 if is2D
%                     ylabel(postprocParams{selected_parameter2},'Interpreter','none');
%                 else
%                     ylabel('frequency');
%                 end
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end
        
        figure;
        histobj = rangedHist(hdata1,hdata2, nBins,dataRange);
        xlabel(xlabelText);
        ylabel(ylabelText);
        
        % Fit if requested. (2D histograms cannot be fitted.)
        if get(h_all.cb_distribution_fit,'Value') && isgraphics(histobj) && ~is2D
            fitHist(histobj);
        end        
    end

    % Switch black-white and hot display mode
    function bwCallback(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        use_bw = ~use_bw;
        
        drawColors(nr_track_colors); % Recompute colors
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay()
        end
        
    end

    % Switch flim and intensity display mode
    function flimCallback(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        use_flim = ~use_flim;
              
        if use_flim
            % Set up the image and axis for FLIM plotting
            set(imagehandle,'AlphaDataMapping','scaled');
            ax = get(imagehandle,'Parent'); 
            set(ax,'Color',[0 0 0]);
            colormap(ax,cmap_isoluminant65());
        else
            % Remove alpha channel
            set(imagehandle,'AlphaData',1);
            set(imagehandle,'AlphaDataMapping','none');
        end
        
        %swap colormap
        cm_ind_temp = cm_ind;
        cm_ind = cm_ind_inactive;
        cm_ind_inactive = cm_ind_temp;
        set(h_all.popup_colormap,'Value',cm_ind);
        setColormap();
        
        %swap zl and zl_inactive
        zl_temp = zl;
        zl = zl_alpha;
        zl_alpha = zl_temp;
        
        caxis(h_all.axes,zl);
        if use_flim
            alim(h_all.axes,zl_alpha);
        end
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay()
        end
        
    end

    function setColormap(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        cm_invert = logical(get(h_all.cb_invert,'Value'));
        cm_ind = get(h_all.popup_colormap,'Value');
        
        if imagehandle ~= -1 
            if isprop(get(imagehandle,'Parent'),'Colormap')
                set(get(imagehandle,'Parent'),'Colormap',getColormap(cm_maps{cm_ind},cm_invert));
            else % In older version the colormap property of the axis is not directly accessible. Can fail with empty (to be loaded) images
                try
                    colormap(get(imagehandle,'Parent'),getColormap(cm_maps{cm_ind}))
                catch
                end
            end
        end
                
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay()
        end
    end

    % Recompute the colors based on the current background
    function drawColors(num_colors)
%         if use_bw
%             bg = {'k'}; % background color
%         else
            bg = {'r'};
%         end
        
        colors = distinguishable_colors(2, bg);
        marker_color = colors(2,:);
        marker_fill_color = colors(1,:);
        if dothandle_cand ~= -1
            set(dothandle_cand,'Color',marker_color, 'MarkerFaceColor',marker_fill_color);
        end
        if dothandle_fit ~= -1
            set(dothandle_fit,'Color',marker_color, 'MarkerFaceColor',marker_fill_color);            
        end
        
        % Draw num_colors colors. If num_colors is less then the number of
        % tracks, we periodically repeat the colors.
        track_colors = repmat( distinguishable_colors(num_colors, bg), ceil(n_tracks/num_colors) ,1);
        track_colors = track_colors(1:n_tracks,:);
    end

    % Update image gamma. Gamma always acts only the intensity component,
    % not on the lifetime.
    function gammaCallback(hObj, ~)
        callback_FloatEdit(hObj,[], 1e-3, 1e3);
        value = str2double(get(hObj,'String'));
        if ~isnan(value)
            old_gamma = img_gamma;
            img_gamma = value;
            if use_flim
                zl_alpha = zl_alpha.^(img_gamma/old_gamma);                
                alim(h_all.axes,zl_alpha);
            else
                zl = zl.^(img_gamma/old_gamma);
                caxis(h_all.axes,zl);
            end
            
            updateFrameDisplay();
        end
    end

    % Update the movie FPS
    function fpsCallback(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        FPS = str2double(get(h_all.edit_FPS, 'String'));
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
    function sliderCallback(~, ~)
        updateFrameDisplay();
        elapsed_time = 0;
    end

    % This is called continously when dragging the slider
    function updateSlider(~,~)
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

    % Sets text on top of the visualizer according to the current frame
    function updateTopText()
        if firstFrame == 1
            set(h_all.toptext,'String',sprintf('Fr. %i/%i',frame,size(movie,3)));
        else
            set(h_all.toptext,'String',sprintf('Fr. %i/%i (%i)',frame,size(movie,3), frame+firstFrame-1));
        end
    end

    % Prompt to jump to a given frame number
    function callback_Goto(~,~)
        uians = inputdlg({'Frame number'},'Go to',[1 40],{num2str(frame)});
        if ~isempty(uians)
            tempframe = str2double(uians{1});
            if ~isnan(tempframe)
                frame = max(1,min(size(movie,3),tempframe));
                updateFrameDisplay();
                updateTopText();
            end
        end
    end
%% Data Filtering
    % Apply filter to the data of the current mode
    function applyFilter(~,~)
        % this modifies candidateData/refinementData/trackingData/postprocData.
        % Make sure that a copy exists for resetting the filter.
        
        % getname = @(x) inputname(1);
        getFilterFun = @(paramsNames)str2func(['@(posData)true(size(posData,1),1)&' regexprep(vectorize(['( ' h_all.edit_filter.String ' )']),...
            strcat('(?<!\w)',matlab.lang.makeValidName(paramsNames),'(?!\w)'),... % Replace spaces with underscore and makes sure the match is not within a name.
            cellfun(@(n)sprintf('posData(:,%i)',n),num2cell(1:numel(paramsNames)),'UniformOutput',false),...
            'ignorecase')]);
        try
            switch mode
                case 'movie'
                    % Do nothing. UI should be hiden in this case.
                case 'candidate'
                    filter_ind = cellfun(getFilterFun(candidateParams),candidateData,'UniformOutput',false); % get logical vector of applied filter
                    candidateData = cellfun(@(d,ind)d(ind,:),candidateData,filter_ind,'UniformOutput',false);% apply to data
                    candidateDataFilter = cellfun(@(f_all,f_new)accumarray(find(f_all),f_new,size(f_all)),candidateDataFilter,filter_ind,'UniformOutput',false); % merge new filter with previous filters to keep track of relation to unfiltered data
                    allFramesCandidateData = [];
                case 'refinement'
                    filter_ind = cellfun(getFilterFun(refinementParams),refinementData,'UniformOutput',false);
                    refinementData = cellfun(@(d,ind)d(ind,:),refinementData,filter_ind,'UniformOutput',false);
                    refinementDataFilter = cellfun(@(f_all,f_new)accumarray(find(f_all),f_new,size(f_all)),refinementDataFilter,filter_ind,'UniformOutput',false);
                    allFramesRefinementData = [];
                case 'tracking'
                    filter_ind = feval(getFilterFun(trackingParams),trackingData);
                    trackingData=trackingData(filter_ind,:);
                    trackingDataFilter(trackingDataFilter) = filter_ind;
                    updateTrackingCache();
                case 'postproc'
                    filter_ind = feval(getFilterFun(postprocParams),postprocData);
                    postprocData=postprocData(filter_ind,:);
                    postprocDataFilter(postprocDataFilter) = filter_ind;
            end
            
            h_all.edit_filter.ForegroundColor = [0 0.5 0];
            set(h_all.edit_filter,'KeyPressFcn',@callback_resetFilterColor);
            
            plotFrame(frame);
        catch
            h_all.edit_filter.ForegroundColor = [1 0 0];
            set(h_all.edit_filter,'KeyPressFcn',@callback_resetFilterColor);
        end
    end    

    function resetFilter(~,~)
        switch mode
            case 'movie'
                % Do nothing. UI should be hiden in this case.
            case 'candidate'
                candidateData = candidateDataUnfiltered;                
                candidateDataFilter = cellfun(@(d)true(size(d,1),1),candidateDataUnfiltered,'UniformOutput',false);
                allFramesCandidateData = [];
            case 'refinement'
                refinementData = refinementDataUnfiltered;
                refinementDataFilter = cellfun(@(d)true(size(d,1),1),refinementDataUnfiltered,'UniformOutput',false);
                allFramesRefinementData = [];
            case 'tracking'
                trackingData = trackingDataUnfiltered;
                trackingDataFilter = true(size(trackingDataUnfiltered,1),1);
                updateTrackingCache();
            case 'postproc'
                postprocData = postprocDataUnfiltered;
                postprocDataFilter = true(size(postprocDataUnfiltered,1),1);
        end
        h_all.edit_filter.ForegroundColor = [0 0 0];
        set(h_all.edit_filter,'KeyPressFcn',[]);
        
        % Reapply current drift correction (also updates frame display)
        applyDrift(true,false);
    end

    function updateTrackingCache()
        % Prepare tracking data for use with the visualizer
        % Put data for each track into its own cell. This is much more efficient comfortable for plotting.
         % id_tracks and n_tracks can change by filtering
        id_tracks = unique(trackingData(:,1)); % Note: (in case track IDs go from 1 to N without missing numbers, the track index is identical to the tracks ID.
        n_tracks = numel(id_tracks);
        
        cell_traj = cell(n_tracks ,1);
        cnt = 1;
        for iTrack = 1:n_tracks
            cell_traj{iTrack} = trackingData( trackingData(:,1)== id_tracks(cnt) , 2:end);
            cnt = cnt+1;
        end
        compute_tracksVisibleInFrame();
    end

    % Insert the variable name selected in the popup in the filter edit
    % field at the end.
    function callback_filter_insertParam(~,~)
        newParam = h_all.popup_filterParam.String{h_all.popup_filterParam.Value};
        if isempty(h_all.edit_filter.String)
            h_all.edit_filter.String = newParam;
        else
            h_all.edit_filter.String = [h_all.edit_filter.String newParam];
        end
        uicontrol(h_all.edit_filter);
    end
    
    % Helper function to turn text color back to black when modifying a
    % filter.
    function callback_resetFilterColor(obj,evt)
        set(obj,'ForegroundColor',[0 0 0]);
        set(obj,'KeyPressFcn',[]);
    end

    function [data,filter,header] = generateList(mode,showFiltered)
        data = [];
        filter = [];
        header = {};
        switch mode
            case 'movie'
                % Do nothing. UI should be hidden/disabled in this case.
            case 'candidate'
                if showFiltered
                    data = cellfun(@(d,n)[n.*ones(size(d,1),1) d],candidateDataUnfiltered,num2cell((1:numel(candidateDataUnfiltered))'),'UniformOutput',false);
                    data = vertcat(data{:});
                    filter = vertcat(candidateDataFilter{:});
                else
                    data = cellfun(@(d,n)[n.*ones(size(d,1),1) d],candidateData,num2cell((1:numel(candidateData))'),'UniformOutput',false);
                    data = vertcat(data{:});
                end
                header = [{'Frame'};candidateParams];
            case 'refinement'
                if showFiltered
                    data = cellfun(@(d,n)[n.*ones(size(d,1),1) d],refinementDataUnfiltered,num2cell((1:numel(refinementDataUnfiltered))'),'UniformOutput',false);
                    data = vertcat(data{:});
                    filter = vertcat(refinementDataFilter{:});
                else
                    data = cellfun(@(d,n)[n.*ones(size(d,1),1) d],refinementData,num2cell((1:numel(refinementData))'),'UniformOutput',false);
                    data = vertcat(data{:});
                end
                header = [{'Frame'};refinementParams];
            case 'tracking'
                if showFiltered
                    data = trackingDataUnfiltered;
                    filter = trackingDataFilter;
                else
                    data = trackingData;
                end
                header = trackingParams;
            case 'postproc'
                if showFiltered
                    data = postprocDataUnfiltered;
                    filter = postprocDataFilter;
                else
                    data = postprocData;
                end
                header = postprocParams;
        end
    end
    function showList(obj,evt)   
        modifiers = get(h_main,'currentModifier');
        ctrlIsPressed = ismember('control',modifiers);
        showFiltered = ctrlIsPressed;
        if showFiltered
            figTitle = sprintf('Localisation Data (%s, unfiltered)',mode);
        else
            figTitle = sprintf('Localisation Data (%s, filtered)',mode);
        end
        
        hdl_table = uitable('Parent',figure('Name',figTitle,'NumberTitle','off','ToolBar','none','MenuBar','none'),'Units','normalized','Position',[0 0 1 1],'RowName',[]);
        
        [data,filter,header] = generateList(mode,showFiltered);
        
        hdl_table.Data = data;
        if ~isempty(filter)
            hdl_table.BackgroundColor = padarray([1 1 1; 0.94 0.94 0.94],size(data,1)-2,'circular','post')-0.3+0.3*vertcat(filter);
        end
        hdl_table.ColumnName = header;
    end
    function exportList(~,~)
        modifiers = get(h_main,'currentModifier');
        ctrlIsPressed = ismember('control',modifiers);
        showFiltered = ctrlIsPressed;        
        
        [exportFile,exportPath] = uiputfile({'*.csv';'*.tsv';'*.xlsx'}, 'Export data...');
        
        if ~isnumeric(exportFile)
            exportFile = fullfile(exportPath,exportFile);
            [~,~,ext] = fileparts(exportFile);
            
            [data,filter,header] = generateList(mode,showFiltered);
            
            
            % writetable is working but very slow and requires sanitisation
            % of paramNames
            data = array2table(data,'VariableNames',matlab.lang.makeValidName(header));
            switch ext
                case 'csv'
                    writetable(data,exportFile,'Delimiter',',','FileType','text');                    
                case 'tsv'
                    writetable(data,exportFile,'Delimiter','tab','FileType','text');
                otherwise
                    writetable(data,exportFile);
            end
        end
        
    end

%% Drift correction
    function calcDrift(~,~)
        seglen = str2double(get(h_all.edit_drift_seg,'String'));
        mag = str2double(get(h_all.edit_drift_res,'String'));
        rmax = str2double(get(h_all.edit_drift_rmax,'String'));
        driftMode = get(h_all.popup_drift,'Value'); 
        
        pixelSize = 1;
        binSize = pixelSize/mag;
        
        movieSize = size(movie);
        movieSize = max(movieSize(1:2));
        
        switch driftMode
            case 1 % RCC
                driftfun = @(coords)RCC(coords,seglen,movieSize,pixelSize,binSize,rmax);
            case 2 % DCC
                driftfun = @(coords)DCC(coords,seglen,movieSize,pixelSize,binSize);
            case 3 % MCC
                driftfun = @(coords)MCC(coords,seglen,movieSize,pixelSize,binSize);
            case 4 % BaSDI
                driftfun = @(coords)calcDriftBaSDI(coords,binSize);
        end
        
        enabledUI = findall(h_main,'Enable','on','Type','UIControl','-not','Style','text');
        set(enabledUI,'Enable','off');
        drawnow();
        posdata = [];
        switch mode
            case 'movie'
                % no localsiations to process
                return
            case 'candidate'
                %TODO
            case 'refinement'
                %TODO
            case 'tracking'
                posdata = [trackingData(:,3:4) + drifttrack(trackingData(:,2),:),trackingData(:,2)]; % Remove any old drift correction
            case 'postproc'
                posdata = [postprocData(:,3:4) + driftpost(postprocData(:,2),:),postprocData(:,2)]; % Remove any old drift correction
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end
        if ~isempty(posdata)
            % Ensure that the drift corrrection can determine the true
            % number of frames
            if posdata(end,:)<size(movie,3)
                posdata(end+1,:) = [nan nan size(movie,3)];
            end
            try % Avoid blocking of the UI on error
                [~,driftcalc] = driftfun(posdata);
            catch err
                disp( getReport( err, 'extended', 'hyperlinks', 'on' ) )
            end
        end
        set(enabledUI,'Enable','on');
    end

    function [img,drift] = calcDriftBaSDI(coords,binSize)
        % Genrates the same output format as for RCC, etc
        BaSDI_struct = BaSDI(coords,binSize);
        drift = processing_result(BaSDI_struct.g).*binSize;
        img = BaSDI_struct.theta;
    end
    function callback_showDrift(~,~)
        modifiers = get(h_main,'currentModifier');
        ctrlIsPressed = ismember('control',modifiers);
        if isempty(driftcalc)
            driftcalc = zeros(size(movie,3),2);
        end
        fig = figure;
        ax = axes(fig);
        if ~ctrlIsPressed
            %%
            hdl = plot(ax,(1:size(driftcalc,1))',driftcalc(:,1),(1:size(driftcalc,1))',driftcalc(:,2));
            legend({'x','y'});
            xlabel('Frame');
            ylabel('Drift (px)');
            title('Drift with respect to first frame');
            %%        
            driftcurrent = false; 
            switch mode
                case 'movie'
                    % no localsiations to process
                    return
                case 'candidate'
                    if ~isempty(driftcand)
                        driftcurrent = driftcand;
                    end
                case 'refinement'
                    if ~isempty(driftref)
                        driftcurrent = driftref;
                    end
                case 'tracking'
                    if ~isempty(drifttrack)
                        driftcurrent = drifttrack;
                    end
                case 'postproc'
                    if ~isempty(driftpost)
                        driftcurrent = driftpost;
                    end
                otherwise
                    error('Unkown display mode ''%s''!',mode);
            end
            if numel(driftcurrent)==numel(driftcalc) && any(driftcurrent(:)~=driftcalc(:))
                hold all;
                ax.ColorOrderIndex = 1;
                plot(ax,1:size(driftcurrent,1),driftcurrent(:,1),1:size(driftcurrent,1),driftcurrent(:,2),'LineStyle','--');
                legend({'x (calculated)','y (calculated)','x (current)','y (current)'});
                hold(h_all.axes,'off');
            end
            
        else
            %%
            hdl = scatter(ax,driftcalc(:,2), driftcalc(:,1),1,1:size(driftcalc,1));
            title('Drift with respect to first frame');
            xlabel('x (px)');
            ylabel('y (px)');
            axis image;
            axis ij; % flip y-axis upside down
            colorbar();
        end
    end

    function callback_importDrift(~,~)
        [importFile,importPath] = uigetfile({'*.csv';'*.tsv';'*.xlsx'}, 'Import drift...');
        
        if ~isnumeric(importFile)
            try
                wstate = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
                datatab = readtable(fullfile(importPath,importFile),'ReadVariableNames',true);
                warning(wstate);
            catch
                warning('Cannot read csv file.');
                return
            end
            
            containsIsolated = @(str,x)~cellfun(@isempty,regexpi(str, ['(?<![a-zA-Z0-9])' x '(?![a-zA-Z0-9])'],'forcecelloutput'));
            
            xpos = find(containsIsolated(datatab.Properties.VariableNames,'x'));
            ypos = find(containsIsolated(datatab.Properties.VariableNames,'y'));
            fpos = find(containsIsolated(datatab.Properties.VariableNames,'frame'));
            
            if ~isempty(xpos)&&~isempty(ypos)&&~isempty(fpos)
                data = datatab.Variables;
                driftcalc = interp1(data(:,fpos),data(:,[xpos ypos]),(1:size(driftcalc,1)),'spline','extrap');
            end
        end
    end
    
    function callback_exportDrift(~,~)
        % ToDo
        [exportFile,exportPath] = uiputfile({'*.csv';'*.tsv';'*.xlsx'}, 'Export drift...');
        
        if ~isnumeric(exportFile)
            exportFile = fullfile(exportPath,exportFile);
            [~,~,ext] = fileparts(exportFile);
            
            data = [(1:size(driftcalc,1))',driftcalc];            
            
            % writetable is working but very slow and requires sanitisation
            % of paramNames
            data = array2table(data,'VariableNames',{'frame','x','y'});
            switch ext
                case 'csv'
                    writetable(data,exportFile,'Delimiter',',','FileType','text');                    
                case 'tsv'
                    writetable(data,exportFile,'Delimiter','tab','FileType','text');
                otherwise
                    writetable(data,exportFile);
            end
        end
    end
    
    function callback_applyDrift(~,~,varargin)
        applyDrift(varargin{:})
    end

    function applyDrift(driftnew,driftold)
        % applyDrift(true)          -> changes drift to driftcalc
        % applyDrift(false)         -> resets drift to 0
        % applyDrift(true,false)    -> applies driftcalc
        if numel(driftnew)==1 % Exspanned singletons (e.g. 0)
            if islogical(driftnew) && driftnew
                driftnew = driftcalc;
            else
                driftnew = driftnew .* ones(size(movie,3),2);
            end
        end
        if ~exist('driftold','var') || isempty(driftold)
            driftold = true;
        elseif ~all(size(driftold)==size(driftnew)) % Exspanned singletons (e.g. 0)
            driftold = driftold .* ones(size(driftnew));
        end
        switch mode
            case 'movie'
                % no localsiations to process
                return
            case 'candidate'
                %TODO
            case 'refinement'
                %TODO
            case 'tracking'
                if islogical(driftold) && driftold
                    driftold = drifttrack;
                end
                trackingData(:,3:4) = trackingData(:,3:4) + driftold(trackingData(:,2),:) - driftnew(trackingData(:,2),:);
                drifttrack = driftnew;
                updateFrameDisplay();
            case 'postproc'
                if islogical(driftold) && driftold
                    driftold = driftpost;
                end
                postprocData(:,3:4) = postprocData(:,3:4) + driftold(postprocData(:,2),:) - driftnew(postprocData(:,2),:);
                driftpost = driftnew;
                updateFrameDisplay();
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end

    end
%% Reconstruction
    function callback_reconstruct(~,~,sthist_zParam)
        recMode = get(h_all.popup_reconstruct_mode,'Value');
        pixelSize = str2double(get(h_all.edit_reconstruct_pixelsize, 'String')); %nm
        sthist_superRes =  str2double(get(h_all.edit_reconstruct_res,'String'));
        sthist_locprec =  str2double(get(h_all.edit_reconstruct_locprec,'String'));
        
        avgTracks  = get(h_all.cb_reconstruct_avgTrack,'Value');
        if isnan(sthist_locprec)
            sthist_guessLocPrecision = true;
            sthist_locprec = [];
        else
            sthist_guessLocPrecision = false;
            sthist_locprec = sthist_locprec/pixelSize; % convert to px
        end
        if nargin>2 && numel(sthist_zParam)==1 && sthist_zParam
            sthist_zParam = get(h_all.popup_reconstruct_mean,'Value');
        else
            sthist_zParam = false;
        end
        sthist_z = [];
        
        % Generate locData = [x,y,amp,bg,sigma]
        sigma_astigmatic = false;
        switch mode
            case 'movie'
                % no localsiations to process
                return
            case 'candidate'
                % Save concatenated data of all frames (done only once)
                if(isempty(allFramesCandidateData))
                    allFramesCandidateData = vertcat(candidateData{:});
                end
                locData = allFramesCandidateData(:,1:2);
                avgTracks = false;
                
                recMode = 1; % We do not have amplitudes
                sthist_guessLocPrecision = false;
                if sthist_zParam
                    sthist_z = allFramesCandidateData(:,sthist_zParam);
                end
            case 'refinement'
                % Save concatenated data of all frames (done only once)
                if(isempty(allFramesRefinementData))
                    allFramesRefinementData = vertcat(refinementData{:});
                end
                if strcmpi('sigma_x',refinementParams{6}) && size(allFramesRefinementData,2)>6
                    sigma_astigmatic = true;
                    locData = allFramesRefinementData(:,[1 2 4 5 6 7]);
                else
                    locData = allFramesRefinementData(:,[1 2 4 5 6]);
                end
                if sthist_zParam
                    sthist_z = allFramesRefinementData(:,sthist_zParam);
                end
                avgTracks = false;
            case 'tracking'
                if strcmpi('sigma_x',trackingParams{8}) && size(trackingData,2)>8
                    sigma_astigmatic = true;
                end
                if avgTracks
                    [~,~,ind] = unique(trackingData(:,1)); % The unique ensures that ind starts at one and has no gaps
                    locData = zeros(max(ind),5);
                    locData(:,1) = accumarray(ind,trackingData(:,3),[],@mean);
                    locData(:,2) = accumarray(ind,trackingData(:,4),[],@mean);
                    locData(:,3) = accumarray(ind,trackingData(:,6),[],@sum);
                    locData(:,4) = accumarray(ind,trackingData(:,7),[],@sum);
                    locData(:,5) = accumarray(ind,trackingData(:,8),[],@mean);
                    if sigma_astigmatic
                        locData(:,6) = accumarray(ind,trackingData(:,9),[],@mean);
                    end
                    if sthist_zParam
                        sthist_z = accumarray(ind,trackingData(:,sthist_zParam),[],@mean);
                    end
                else
                    if sigma_astigmatic
                        locData = trackingData(:,[3 4 6 7 8 9]);
                    else
                        locData = trackingData(:,[3 4 6 7 8]);
                    end
                    if sthist_zParam
                        sthist_z = trackingData(:,sthist_zParam);
                    end
                end
            case 'postproc'
                if strcmpi('sigma_x',postprocParams{8}) && size(postprocData,2)>8
                    sigma_astigmatic = true;
                end
                if avgTracks
                    [~,~,ind] = unique(postprocData(:,1)); % The unique ensures that ind starts at one and has no gaps
                    locData = zeros(max(ind),5);
                    locData(:,1) = accumarray(ind,postprocData(:,3),[],@mean);
                    locData(:,2) = accumarray(ind,postprocData(:,4),[],@mean);
                    locData(:,3) = accumarray(ind,postprocData(:,6),[],@sum);
                    locData(:,4) = accumarray(ind,postprocData(:,7),[],@sum);
                    locData(:,5) = accumarray(ind,postprocData(:,8),[],@mean);
                    if sigma_astigmatic
                        locData(:,6) = accumarray(ind,postprocData(:,9),[],@mean);
                    end
                    if sthist_zParam
                        sthist_z = accumarray(ind,postprocData(:,sthist_zParam),[],@mean);
                    end
                else
                    if sigma_astigmatic
                        locData = postprocData(:,[3 4 6 7 8]);
                    else
                        locData = postprocData(:,[3 4 6 7 8 9]);
                    end
                    if sthist_zParam
                        sthist_z = postprocData(:,sthist_zParam);
                    end
                end
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end
        if sthist_guessLocPrecision && size(locData,2)>4
            if sigma_astigmatic
                sigma_sq = mean(locData(:,5)+locData(:,6),2).^2;
            else
                sigma_sq = locData(:,5).^2;
            end
            N = locData(:,3).*2*pi.*sigma_sq; %total number of photons
            tau = 2*pi*locData(:,4).*(sigma_sq+1/12)./N;
            sthist_locprec = sqrt((sigma_sq+1/12)./N.*(1+4*tau+sqrt(2*tau./(1+4*tau)))); %Rieger et al, DOI 10.1002/cphc.201300711
        end
        if sthist_zParam
            locData = [locData(:,1:2) sthist_z];
            sthist_map =  cell(2,1);
        else
            locData = locData(:,1:2);
            sthist_map =  cell(1,1);
        end
        switch recMode
            case 1 % Gaussian
                recMode = 'gaussian';
                sthist_weight = [];
            case 2 % Weighted
                recMode = 'point';
                sthist_weight = sthist_locprec;
                sthist_locprec = [];                
            case 3 % Jitter
                recMode = 'jitter';
                sthist_weight = [];                
            otherwise % Points
                recMode = 'point';
                sthist_weight = [];
        end
        [sthist_map{:}] = reconstructSMLM(locData,sthist_locprec,sthist_weight,sthist_superRes,size(movie),recMode);
                
        TNTvisualizer(sthist_map,struct('title','Reconstruction','metadata',struct('pixelsize',pixelSize/sthist_superRes*1e-3,'pixelsize_unit',[char(181) 'm'])));

    end

%% I/O
    % Parse input variables and set them to their default values if they
    % are not given. This is also responsible for loading the movie or a
    % TNT file if these are given as input. Additionally it prepares the
    % trackingData for plotting and disables buttons of modes where the
    % corresponding data is missing.
    function parse_inputs_and_setup(num_argin)
        try  % convertStringsToChars is not availabe in older matlab version. 
            movieOrTNTfile = convertStringsToChars(movieOrTNTfile);
            if num_argin>1
                candidateDataOrTNTfile = convertStringsToChars(candidateDataOrTNTfile);
            end
        catch 
        end
        
        fullPathToThisFile = mfilename('fullpath');
        [TNTpath,~,~] = fileparts(fullPathToThisFile);
        
        % Set unspecified input to empty
        if num_argin<2 || isempty(candidateData)
            candidateData = {};
        end
        if num_argin<3
            candidateParams = {};
        end
        if num_argin<4 || isempty(refinementData)
            refinementData = {};
        end
        if num_argin<5
            refinementParams = {};
        end
        if num_argin<6
            trackingData = [];
        end
        if num_argin<7
            trackingParams ={};
        end
        if num_argin<6
            postprocData = [];
        end
        if num_argin<7
            postprocParams ={};
        end
        if num_argin<10
            FPS = [];
        end
        if num_argin<11
            is_blocking = [];
        end
        
        if num_argin<12
            firstFrame = [];
        end
        
        loadMovie(movieOrTNTfile,num_argin);
        
        if iscell(movieOrTNTfile)||isnumeric(movieOrTNTfile)
            % Clear movie to save memory
            movieOrTNTfile = [];
        end

        % Check if candidateData or TNT file was given as second argument
        if num_argin>1
            if ischar(candidateDataOrTNTfile)||isstruct(candidateDataOrTNTfile) % TNT file/struct was given             
                if isstruct(candidateDataOrTNTfile)
                    TNTdata = candidateDataOrTNTfile;
                    if isfield(TNTdata,'filename_movie')
                        [~,titlename] = fileparts(TNTdata.filename_movie);
                        titlename = ['Analysis of ', titlename];
                    end
                else
                    fprintf('Loading TNT file..\n')
                    % Before loading, The plugin subfolder path is removed, because loading an XXXoptions
                    % struct with a (main/init/post) function handle where the plugin file exists,
                    % but the subfunction does not (e.g. because it was renamed) throws an
                    % error on loading the file
                    TNTpluginpath = genpath([TNTpath,filesep,'plugins']);
                    s = warning('off','all');
                    rmpath(TNTpluginpath);
                    TNTdata = load(candidateDataOrTNTfile);
                    warning(s);
                    % Add plugin path again (otherwise we loading a TNT
                    % file gives warning that the function to the function
                    % handles could not be found.
                    addpath(TNTpluginpath);
                    
                    [~,titlename] = fileparts(fileORmovie);
                end
                
                if isfield(TNTdata,'movieSize') && ~isequal(size(movie), TNTdata.movieSize)
                   error('Input movie and given TNT file do not fit together! Movie size [%i,%i,%i], TNT file [%i,%i,%i].\n',size(movie),TNTdata.movieSize) 
                end              
                transferDatafromTNTdata(TNTdata);
            elseif iscell(candidateDataOrTNTfile)
                candidateData = candidateDataOrTNTfile;
            end            
        end
        if ~isempty(titlename)
            h_main.Name = ['TNT Visulizer: ' titlename];
        end        
        
        % -- Initialize variables left empty --
        % Is candidate data available?
        if isempty(candidateData)
            set(h_all.button_candidateMode,'Enable','off');
        else
            candidateParams = checkParameterDescription(candidateData, candidateParams);
            mode = 'candidate';
        end
        
        if isempty(candidateParams)
            candidateParams = {};
        end
        
        % Is refinement data available?
        if isempty(refinementData)
            set(h_all.button_refinementMode,'Enable','off');
        else
            refinementParams = checkParameterDescription(refinementData, refinementParams);
            mode = 'refinement';
        end
        
        if isempty(refinementParams)
            refinementParams = {};
        end
        
        % Is tracking data available?
        if isempty(trackingData)
            set(h_all.button_trackingMode,'Enable','off');
            
            id_tracks = []; % Note: (in case track IDs go from 1 to N without missing numbers, the track index is identical to the tracks ID.
            n_tracks = 0;
        else
            trackingParams = checkParameterDescription(trackingData, trackingParams);
            mode = 'tracking';
            
            id_tracks = unique(trackingData(:,1)); % Note: (in case track IDs go from 1 to N without missing numbers, the track index is identical to the tracks ID.
            n_tracks = numel(id_tracks);
        end
        % Prepare tracking data for use with the visualizer
        % Put data for each track into its own cell. This is much more efficient comfortable for plotting.
%         cell_traj = cell(n_tracks ,1);
%         cnt = 1;
%         for iTrack = 1:n_tracks
%             cell_traj{iTrack} = trackingData( trackingData(:,1)== id_tracks(cnt) , 2:end);
%             cnt = cnt+1;
%         end
        if ~isempty(trackingData)
            trackingData_sorted = sortrows(trackingData,[1 2]);
            cell_traj = mat2cell(trackingData_sorted(:,2:end),accumarray(trackingData_sorted(:,1),1),size(trackingData_sorted,2)-1);
        else
            cell_traj = cell(0 ,1);
        end
        
        if isempty(trackingParams)
            trackingParams ={};
        end        
        % Is postproc data available?
        if isempty(postprocData)
            set(h_all.button_postprocMode,'Enable','off');
        else
            postprocParams = checkParameterDescription(postprocData, postprocParams);
            mode = 'postproc';
        end
        
        if isempty(postprocParams)
            postprocParams ={};
        end
        
        if isempty(FPS)
            FPS = 30;
        end
        
        if isempty(is_blocking)
            is_blocking = false;
        end
        
        if isempty(firstFrame)
            firstFrame = 1;
        end
    end

    % Transfer data from TNTfile into our local GUI variables
    function transferDatafromTNTdata(TNTdata)
        if(isfield(TNTdata,'candidateData'))
            candidateData = TNTdata.candidateData;
            candidateParams = TNTdata.candidateOptions.outParamDescription;
        else
            candidateData = {};
            candidateParams = {};
        end
        if(isfield(TNTdata,'refinementData'))
            refinementData = TNTdata.refinementData;
            refinementParams = TNTdata.refinementOptions.outParamDescription;
        else
            refinementData = {};
            refinementParams = {};
        end
        if(isfield(TNTdata,'trackingData'))
            trackingData = TNTdata.trackingData;
            trackingParams = TNTdata.trackingOptions.outParamDescription;
        else
            trackingData = [];
            trackingParams = {};
        end
        if(isfield(TNTdata,'postprocData'))
            postprocData = TNTdata.postprocData;
            postprocParams = TNTdata.postprocOptions.outParamDescription;
        else
            postprocData = [];
            postprocParams = {};
        end
        if(isfield(TNTdata,'firstFrame_lastFrame'))
            firstFrame = TNTdata.firstFrame_lastFrame(1);
            if isfield(TNTdata,'globalOptions')&&isfield(TNTdata,'binFrame')
                firstFrame = round(firstFrame-1 / TNTdata.globalOptions.binFrame)+1;% Can be of by half a frame. Not critically. Only for display.
            end
        end
        if(isfield(TNTdata,'FPS'))
            FPS = TNTdata.FPS;
        end
        if(isfield(TNTdata,'is_blocking'))
            is_blocking = TNTdata.is_blocking;
        end
        if(isfield(TNTdata,'metadata'))
            metadata = TNTdata.metadata;
        end
        if isfield(TNTdata,'title')
            titlename = TNTdata.title;
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
        
        % If there is no data, do nothing
        if nrParams == 0
            return
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
    
    function loadMovie(fileORmovie,num_argin)
        if isempty(fileORmovie)
            set(h_all.slider,'Enable','off');
            set(h_all.but_play,'Enable','off');
            return;
        end
        if importMode && ischar(fileORmovie) && ~endsWith(fileORmovie,'.mat')
            [importPlugins,importFormats] = loadImportPlugins();
            
            if exist(fileORmovie,'file')
                importPlugin = selectImportPlugin(fileORmovie,importPlugins);
                if importPlugin>0
                    importPlugin = importPlugins(importPlugin);
                else
                    error('Unsupported file type.');
                end
            end
            [~,filename,ext] = fileparts(fileORmovie);
            filename = [filename,ext];
            set(h_all.edit_filename,'String',filename);
            set(h_all.edit_filename,'Tooltip',fileORmovie);
            
            if size(importPlugin.param_specification,1)>0 %if import plugin has options
                set(h_all.panel_importOptions,'Visible','on');
                importPlugin.createOptionsPanel(h_all.panel_importOptions);                
            else
                set(h_all.panel_importOptions,'Visible','off');
            end
            if isfield(importPlugin.info,'hasTCSPC') && importPlugin.info.hasTCSPC
                set(h_all.cbx_useTimegate, 'Enable', 'on');
                if isfield(importPlugin.info,'getTCSPC') && isa(importPlugin.info.getTCSPC,'function_handle')
                    set(h_all.button_showTCSPC, 'Enable', 'on');
                else
                    set(h_all.button_showTCSPC, 'Enable', 'off');
                end
            else
                set(h_all.cbx_useTimegate, 'Enable', 'off');
                set(h_all.cbx_useTimegate, 'Value', false);
                set(h_all.button_showTCSPC, 'Enable', 'off');
            end
        % Check if movie, path to movie or path to TNT file was given as first argument
        elseif(isnumeric(fileORmovie)||iscell(fileORmovie)) % movie OR FLIM movie was given
            movie = fileORmovie;
        elseif (ischar(fileORmovie)||isstruct(fileORmovie)) % path was given
            if ischar(fileORmovie) && ~endsWith(fileORmovie,'.mat','IgnoreCase',true)
                % Not a TNT file. Try importMode for import with Plugin
                importMode = true;
                loadMovie(fileORmovie)
                return;
            end
            
            % Load movie and data
            fprintf('Loading TNT file..\n')
            fullPathToThisFile = mfilename('fullpath');
            [TNTpath,~,~] = fileparts(fullPathToThisFile);

            if isstruct(fileORmovie)
                TNTdata = fileORmovie;
                [~,titlename] = fileparts(TNTdata.filename_movie);
                titlename = ['Analysis of ', titlename];
            else
                % Before loading, The plugin subfolder path is removed, because loading an XXXoptions
                % struct with a (main/init/post) function handle where the plugin file exists,
                % but the subfunction does not (e.g. because it was renamed) throws an
                % error on loading the file
                TNTpluginpath = genpath([TNTpath,filesep,'plugins']);
                s = warning('off','all');
                rmpath(TNTpluginpath);
                TNTdata = load(fileORmovie);
                warning(s);
                % Add plugin path again (otherwise we loading a TNT
                % file gives warning that the function to the function
                % handles could not be found.
                addpath(TNTpluginpath);
                
                [~,titlename] = fileparts(fileORmovie);
            end
            
            fprintf('Loading movie specified in TNT file..\n')
            
            loadMovieWithPlugin(TNTdata,fileORmovie);
            
            transferDatafromTNTdata(TNTdata);
            
            % If TNT file is given as first AND second argument, throw error
            if nargin>1 && num_argin>1
                if(ischar(candidateDataOrTNTfile))
                    error('Invalid input in TNTvisualizer (Two TNT/MAT files).')
                end
            end
        else
            error('Invalid movie input.');
        end
        if ~isempty(movie)
            if iscell(movie) % FLIM movie
                if numel(movie)>1
                    movieLT = movie{2};
                end
                movie = movie{1};
                if ~isempty(movieLT)
                    set(h_all.cb_flim, 'Enable', 'on');
                end
            end
            if(size(movie,3)>1)
                set(h_all.slider,'Value',1, 'Min',1,'Max',size(movie,3),'SliderStep',[1/size(movie,3) 1/size(movie,3)],'Callback', @sliderCallback);
                set(h_all.slider,'Enable','on');
                set(h_all.but_play,'Enable','on');
            else % For single images we disable slider and play button
                set(h_all.slider,'Enable','off');
                set(h_all.but_play,'Enable','off');
            end
        end
    end

    function loadMovieWithPlugin(TNTdata,TNTfile)
        [~,fname,ext] = fileparts(TNTdata.filename_movie);
        if ~exist(TNTdata.filename_movie,'file')
            [filename, path] = uigetfile({['*' ext],'Movie files';'*.*','All files'},['File ' fname ext ' not found.'],TNTdata.filename_movie);
            if isfloat(filename)
                movie = zeros(TNTdata.movieSize);
                warning('File not found: %s',TNTdata.filename_movie);
                return;                
            else
                TNTdata.filename_movie = [path,filename];
            end
        end
        try
            if strcmpi('UNKNOWN Function',func2str(TNTdata.importOptions.mainFunc))
                % This can happen if the path to the plugin folder has
                % changed. In this case load importOptions again. This time
                % with plugins folder in the path.
                temp = load(TNTfile,'importOptions');
                TNTdata.importOptions = temp.importOptions;
            end
            movieArgs = {TNTdata.filename_movie, [TNTdata.globalOptions.firstFrame, TNTdata.globalOptions.lastFrame], TNTdata.globalOptions.binFrame,true};
            if TNTdata.globalOptions.useTimegate
                movieArgs = [movieArgs,{[TNTdata.globalOptions.tgStart TNTdata.globalOptions.tgEnd]}];
            end
            [movie, metadata] = TNTdata.importOptions.mainFunc(TNTdata.importOptions,movieArgs{:});
        catch err
            if any(strcmpi(ext,{'.tif','.tiff'})) % legacy for old TNT files without importOptions (tif only, no binning)
                if(isfield(TNTdata,'firstFrame_lastFrame'))
                    firstFrame = TNTdata.firstFrame_lastFrame(1);
                    movie = read_tiff(TNTdata.filename_movie,false, TNTdata.firstFrame_lastFrame);
                else
                    movie = read_tiff(TNTdata.filename_movie,false);
                end
            else
                movie = zeros(TNTdata.movieSize);
                warning('Cannot load movie from file %s',TNTdata.filename_movie);
                warning(getReport(err,'extended','hyperlinks','on'));
            end
        end
    end

    function generatePreview(~,~)
        importOptions = importPlugin.getOptions();
        
        globalOptions.firstFrame = getNum(h_all.edit_firstFrame);
        globalOptions.lastFrame =  getNum(h_all.edit_lastFrame);
        globalOptions.binFrame =  getNum(h_all.edit_binFrame);
        globalOptions.useTimegate = logical(get(h_all.cbx_useTimegate,'Value'));
        globalOptions.tgStart = getNum(h_all.edit_tgStart);
        globalOptions.tgEnd = getNum(h_all.edit_tgEnd);
        
        
        movieArgs = {movieOrTNTfile, [globalOptions.firstFrame, globalOptions.lastFrame], globalOptions.binFrame,true};
        if globalOptions.useTimegate
            movieArgs = [movieArgs,{[globalOptions.tgStart globalOptions.tgEnd]}];
        end
        enabledUI = findall(h_main,'Enable','on','Type','UIControl','-not','Style','text');
        set(enabledUI,'Enable','off');
        try
            [tempmovie,metadata] = importOptions.mainFunc(importOptions,movieArgs{:});
        catch err
            warning(getReport(err,'extended','hyperlinks','on'));
            tempmovie = [];
        end
        set(enabledUI,'Enable','on');
        loadMovie(tempmovie);
        initPlot();
    end

    function callback_selectFile(~,~)
        modifiers = get(h_main,'currentModifier');
        ctrlIsPressed = ismember('control',modifiers);
        
        if isempty(importFormats)
            [~,importFormats] = loadImportPlugins();
            importFormats{1,1} = [importFormats{1,1},';*.mat'];
            importFormats = [importFormats;{'*.mat','Visualizer files';'*.*','All files'}];
        end
        newpath = pwd;
        if ~isempty(movieOrTNTfile)
            [newpath,~,~] = fileparts(movieOrTNTfile);
        end
        [newMovie, newpath] = uigetfile(importFormats,'Select files to load',[newpath,filesep]);
        if ~isfloat(newMovie)
            if endsWith(newMovie,'.mat') %% TNT file
                if ctrlIsPressed
                    TNTvisualizer([newpath,newMovie]);
                else
                    enabledUI = findall(h_main,'Enable','on','Type','UIControl','-not','Style','text');
                    set(enabledUI,'Enable','off'); % Disable the old GUI since it will be deleted anyway
                    TNTvisualizer([newpath,newMovie]);
                    delete(h_main);
                end
            elseif ctrlIsPressed %% import plugin
                TNTvisualizer([newpath,newMovie]);
            else
                movieOrTNTfile = [newpath,newMovie];
                path = newpath;
                loadMovie(movieOrTNTfile);
                resizeGUIforMode();
            end
        end
    end
    
    % Cleanup function. This is neccessary to delete the timer!
    function onAppClose(~, ~)
        if isfield(h_all, 'timer') % Timer is not created if a error during loading occured
            if strcmp(get(h_all.timer, 'Running'), 'on')
                stop(h_all.timer);
            end
            delete(h_all.timer);
        end
        delete(h_main);
    end

%% Export images
% 
    function exportCallback(~,~)
        modifiers = get(h_main,'currentModifier');
        ctrlIsPressed = ismember('control',modifiers);
        exportMode = get(h_all.popup_export,'Value');
        
        switch exportMode
            case 1 % Tiff (Stack)
                exportFormat = '.tif';
                exportFrames = 1:size(movie,3);
                exportFunction = @(file,img) export_tiff(file,permute(img,[1 2 4 3]),1);

            case 2 % Gif (Stack)
                exportFormat = '.gif';
                exportFrames = 1:size(movie,3);
                exportFunction = @export_gif;
                
            case 3 % Tiff (Frame)
                exportFormat = '.tif';
                exportFrames = frame;
                exportFunction = @(file,img) export_tiff(file,permute(img,[1 2 4 3]),1);
                
            case 4 % Png (Frame)
                exportFormat = '.png';
                exportFrames = frame;
                exportFunction = @(file,img) imwrite(img,file);
                
            case 5 % Tiff (int)
                exportFormat = '.tif';
                exportFrames = NaN;
                exportData = movie;
                exportFunction = @(file,img) export_tiff(file,img,1);
                
            case 6 % Tiff (tau)
                exportFormat = '.tif';
                exportFrames = NaN;
                if isempty(movieLT)
                    warning('No lifetime data available.');
                    return
                end
                exportData = movieLT;
                exportFunction = @(file,img) export_tiff(file,img,1);
                
            case 7 % Copy to new figure window
                nfig = figure('Color',[1 1 1],'Units','pixels','Visible','off');
                temp_unit = h_main.Units;
                h_main.Units = 'pixels';
                nfig.InnerPosition(3:4) = h_main.InnerPosition(3);
                h_main.Units = temp_unit;
                movegui(nfig,'onscreen');
                temp_unit = h_all.axes.Units;
                h_all.axes.Units = 'pixels';
                nax = copyobj([imagehandle.Parent imagehandle.Parent.Legend],nfig);
                nax(1).Position(3:4) = max(h_all.axes.Position(3:4));
                nax(1).Position(1:2) = (nax(1).Parent.InnerPosition(4)-nax(1).Position(4))/2;
                h_all.axes.Units = temp_unit;
                % Recreate colorbar if visible to keep it interactive
                if isgraphics(imagehandle.Parent.Colorbar,'Colorbar')
                    colorbar(nax(1),'Location',imagehandle.Parent.Colorbar.Location);
                end
                % Recreate scalebar if visible to keep it interactive
                delete(findall(nax,'Tag','TNTscalebar_Line','-or','Tag','TNTscalebar_Text'));
                if ~isempty(sb_obj) && strcmpi(sb_obj.Visible,'on')
                    TNTscalebar(findobj(nax,'Type','axes'),'Pixelsize',sb_obj.Pixelsize,'Unit',sb_obj.Unit);
                end
                nfig.Visible = 'on';
                return
            otherwise
                error('Unknown export mode.');
        end
        if ctrlIsPressed
            if numel(exportFrames)>1
                uians = inputdlg({'First frame','Last frame','Step'},'Export range',[1 40],{num2str(exportFrames(1)),num2str(exportFrames(end)),num2str(1)});
                if isempty(uians)
                    return
                else
                    exportFrames = max(1,str2double(uians{1})):max(1,str2double(uians{3})):min(size(movie,3),str2double(uians{2}));
                end
            elseif isnan(exportFrames)
                uians = inputdlg({'First frame','Last frame','Step'},'Export range',[1 40],{num2str(1),num2str(size(exportData,3)),num2str(1)});
                if isempty(uians)
                    return
                else
                    exportData = exportData(:,:,max(1,str2double(uians{1})):max(1,str2double(uians{3})):min(size(movie,3),str2double(uians{2})));
                end                
            end
        end
        
        [exportFile,exportPath] = uiputfile(exportFormat, 'Export image...');
        
        
        if ~isnumeric(exportFile)
            exportFile = fullfile(exportPath,exportFile);
            if isnan(exportFrames)
                exportFunction(exportFile,exportData);
            else
                enabledUI = findall(h_main,'Enable','on','Type','UIControl','-not','Style','text');
                try
                    % try export, to ensure that the UI is reenabled if it fails
                    set(enabledUI,'Enable','off');
                    tempframe = frame;
                    for iframe = 1:numel(exportFrames)
                        if ~isnan(exportFrames(iframe))
                            frame = exportFrames(iframe);
                            updateFrameDisplay();
                            updateTopText();
                        end
                        img(:,:,:,iframe) = export_fig('-m2','-a1','-nocrop',imagehandle.Parent);% This gets the content of the axes with double the screen resolution. Nocrop to ensure equal size of all frames.
                    end
                    exportFunction(exportFile,crop_borders(img,[],-1));%Crop one pixel more to remove border from interpolating
                    
                    frame = tempframe;
                catch err
                    warning(getReport(err,'extended','hyperlinks','on'));
                end
                updateFrameDisplay();
                updateTopText();
                
                set(enabledUI,'Enable','on');
            end
        end
    end

    function export_gif(file,img)
        for idx = 1:size(img,4)
            [A,map] = rgb2ind(img(:,:,:,idx),256);
            if idx == 1
                imwrite(A,map,file,'gif','LoopCount',Inf,'DelayTime',1/FPS);
            else
                imwrite(A,map,file,'gif','WriteMode','append','DelayTime',1/FPS);
            end
        end
    end
%%  Tracking mode related functions

    % Compute which tracks are visible in each frame
    % This is done to iterate only over visible tracks when plotting. Otherwise
    % a large number of tracks in the movie slows us down, even if we skip them
    % on the fly (because we have to check each track for every frame)
    %
    % Sets:
    %  tracksVisibleInFrame
    %  maxNr_visibleTracksInFrame
    function compute_tracksVisibleInFrame()
        nr_frames = size(movie,3);
        tracksVisibleInFrame = cell(nr_frames,1);
        tracksVisibleInFrame_bool = false(n_tracks,nr_frames);
        
        %Note: We need the bool to compute this fast (at the cost of memory)
        % but convert the data to a cell array for using it later
        for iTrack = 1:n_tracks
            track_firstFrame = cell_traj{iTrack}(1,1);
            track_lastFrame = cell_traj{iTrack}(end,1);

            tracksVisibleInFrame_bool(iTrack, track_firstFrame:min(track_lastFrame+traj_lifetime, nr_frames)) = true;
        end
        
        % Convert bool matrix to indices
        for iFrame = 1:nr_frames
            tracksVisibleInFrame{iFrame} = find(tracksVisibleInFrame_bool(:,iFrame)).';
        end
        % Alternative:
        % tracksVisibleInFrame = accumarray(trackingData(:,1),trackingData(:,2),[],@(x){min(x):min(max(x)+traj_lifetime,nr_frames)});

        [~,nr_VisibleTracksInFrame] = cellfun(@size,tracksVisibleInFrame);
        maxNr_visibleTracksInFrame = max(nr_VisibleTracksInFrame);
    end

    % Simply plots the frame with the maximum number of visible tracks,
    % which initializes all linehandles needed for displaying the movie. If
    % this is not called, linehandles are simply initialized on the fly
    % (which might cause jerky playback on playing the movie)
    %
    % Make sure updateFrameDisplay() is called afterwards, so that the
    % current frame is displayed again!
    function initializeLinehandles()
       [~, maxFrame] = max(nr_VisibleTracksInFrame);
       plotFrame(maxFrame);
    end

     % Sets the traj_lifetime variable if it changes in the GUI
    function callback_TrajLifetime(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        traj_lifetime = round(str2double(get(h_all.edit_lifetime,'String')));
        if traj_lifetime<=0 || isempty(traj_lifetime)
            traj_lifetime = 0;
        end
        set(h_all.edit_lifetime,'String',sprintf('%i%',traj_lifetime));
        
        % We have to update which tracks are visible in which frames
        % This alters the number of linehandles which are needed for
        % display (the maximum number of concurrently visible tracks)
        compute_tracksVisibleInFrame();
        resetGraphics();
        initializeLinehandles();
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end

    % Sets the traj_displayLength variable if it changes in the GUI
    function callback_dispLength(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        traj_displayLength = round(str2double(get(h_all.edit_trajDisplayLength,'String')));
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
    function callback_trackColors(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        nr_track_colors = round(str2double(get(h_all.edit_colors,'String')));
        if nr_track_colors<=0 || isempty(nr_track_colors)
            nr_track_colors = 1;
        end
        set(h_all.edit_colors,'String',sprintf('%i%',nr_track_colors));
        
        drawColors(nr_track_colors);
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end
    %% Callbacks for import mode
    % Show the time gates and the TCSPC histogram
    function callback_showTCSPC(~, ~)
        importOptions = importPlugin.getOptions();
        if isfield(importOptions.info,'getTCSPC') && isa(importOptions.info.getTCSPC,'function_handle')
            figure;
            showTimegate(movieOrTNTfile,importOptions,[str2double(get(h_all.edit_tgStart,'String')),str2double(get(h_all.edit_tgEnd,'String'))]);
        else
            fprintf('TNT: The import plugin %s does not support preview of the TCSPC.\n', importOptions.plugin_name)
        end
    end

    function callback_toggleTimegate(~,~)
        if get(h_all.cbx_useTimegate, 'Value')
            set(h_all.text_tgStart, 'Enable', 'on');
            set(h_all.text_tgEnd, 'Enable', 'on');
            set(h_all.edit_tgStart, 'Enable', 'on');
            set(h_all.edit_tgEnd, 'Enable', 'on');
            set(h_all.button_showTCSPC, 'Enable', 'on');
        else
            set(h_all.text_tgStart, 'Enable', 'off');
            set(h_all.text_tgEnd, 'Enable', 'off');
            set(h_all.edit_tgStart, 'Enable', 'off');
            set(h_all.edit_tgEnd, 'Enable', 'off');
            set(h_all.button_showTCSPC, 'Enable', 'off');
        end
    end

    function callback_sendToTNT(~,~)
        modifiers = get(h_main,'currentModifier');
        ctrlIsPressed = ismember('control',modifiers);
        
        if ~isempty(importPlugin) && ~isempty(movieOrTNTfile)
            
            importOptions = importPlugin.getOptions();
            
            globalOptions.firstFrame = getNum(h_all.edit_firstFrame);
            globalOptions.lastFrame =  getNum(h_all.edit_lastFrame);
            globalOptions.binFrame =  getNum(h_all.edit_binFrame);
            globalOptions.useTimegate = logical(get(h_all.cbx_useTimegate,'Value'));
            globalOptions.tgStart = getNum(h_all.edit_tgStart);
            globalOptions.tgEnd = getNum(h_all.edit_tgEnd);
            
            if ~ctrlIsPressed
                onAppClose();
            end
            RunTrackNTrace(movieOrTNTfile,'globalOptions',globalOptions,'importOptions',importOptions);
        end
    end
        
end


%% --- General functions ---

% Adds pathes needed for the visualizer.
function addPathsVisualizer()
    fullPathToThisFile = mfilename('fullpath');
    [path,~,~] = fileparts(fullPathToThisFile);
    addpath(genpath([path,filesep,'subfun']));
    addpath(genpath([path,filesep,'external']));
    addpath(genpath([path,filesep,'helper']));
end

% Function that cuts data from upper and lower tails of the distribution
% keeping at least 'percentOfData' percent of all values.
function histobj = rangedHist(data1,data2, nbins, percentOfData)
is2D = numel(nbins)==2 && ~isempty(data2);
if(percentOfData<100)
    % Limits for the cumulative density function (which goes from 0 to 1)
    lower_limit = (1-percentOfData/100)/2; % below this we throw away
    upper_limit = 1-(1-percentOfData/100)/2; % above this we throw away
    
    [data1,ord] = sort(data1);
    minIdx = round(lower_limit*numel(data1));
    if minIdx < 1
       minIdx = 1; 
    end
    maxIdx = round(upper_limit*numel(data1));
    if maxIdx < 1
       maxIdx = 1; 
    end
    data1 = data1(minIdx:maxIdx);
    if ~isempty(data2)
        data2 = data2(ord);
        data2 = data2(minIdx:maxIdx);
    end
    
    if is2D    % Limits for the cumulative density function (which goes from 0 to 1)
        [data2,ord] = sort(data2);
        minIdx = round(lower_limit*numel(data2));
        if minIdx < 1
            minIdx = 1;
        end
        maxIdx = round(upper_limit*numel(data2));
        if maxIdx < 1
            maxIdx = 1;
        end
        data2 = data2(minIdx:maxIdx);
        data1 = data1(ord);
        data1 = data1(minIdx:maxIdx);
    end
end

if is2D
    histobj = histogram2(data1,data2,nbins,'DisplayStyle','tile');
else
    histobj = histogramw(data1,data2,nbins,'Normalization','count');
end

end

function gaussresmono = fitHist(h)
    bincenter = movmean(h.BinEdges,2,'Endpoints','discard'); % Bin Centres
    hsum = sum(h.Values);
    hmean = sum(h.Values.*bincenter)./hsum;
    hstd = sqrt(abs(sum(h.Values.*bincenter.^2)-(sum(h.Values.*bincenter)).^2/hsum)/hsum)+eps; % +eps ensures sigma>0
    fmono = fittype(@(mu,sigma,A,x)A*exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi)));
    gaussresmono = fit(bincenter',h.Values',fmono,'Start',[hmean hstd hsum]);
    
    
    hold(h.Parent,'on');
    plot(linspace(h.BinLimits(1),h.BinLimits(2)),gaussresmono(linspace(h.BinLimits(1),h.BinLimits(2))),'LineWidth',1.5)
    text(feval(@(x)x(2)-0.05*diff(x),xlim(h.Parent)),feval(@(x)x(2)-0.05*diff(x),ylim(h.Parent)),{'','','',sprintf('Histogram of %g values',sum(h.Values)),sprintf('Data: %.2f \\pm %.2f',hmean, hstd),sprintf('Fit: %.2f \\pm %.2f',gaussresmono.mu,abs(gaussresmono.sigma))},'HorizontalAlignment','right');
    
    hold(h.Parent,'off');
end

%% Special UI elements

% Callback for edit fields containing floats. Checks if a correct
% number was entered and restricts it to the given bounds.
    function callback_FloatEdit(hObj,~, minVal, maxVal, flag)
        if nargin<3 || isempty(minVal)
            minVal=-inf;
        end
        if nargin<4 || isempty(maxVal)
            maxVal=inf;
        end
        
        value = str2num(get(hObj, 'String')); %#ok<ST2NM> str2num distinguishes NaN and invalid input, str2double does not.
        if isempty(value)
            set(hObj,'ForegroundColor','r');
            set(hObj,'String','INVALID');
            uicontrol(hObj);
        elseif nargin>4 && isnan(value) && ischar(flag) && strcmpi(flag,'includeNaN')
            set(hObj,'ForegroundColor','k');
            set(hObj,'String','NaN');
        else
            value = max(minVal,value);
            value = min(maxVal,value);
            set(hObj,'ForegroundColor','k');
            set(hObj,'String',sprintf('%.2f',value));
        end
    end

% Callback for edit fields containing integer values. Checks if a correct
% number was entered and restricts it to the given bounds.
    function callback_IntEdit(hObj,~, minVal,maxVal)
        if nargin<3 || isempty(minVal)
            minVal=0;
        end
        if nargin<4 || isempty(maxVal)
            maxVal=inf;
        end
        
        % Accept end as inf, as MATLAB users are used to end as the last element
        if(strcmp(get(hObj,'String'), 'end'))
            value = inf;
        else
            value = round(str2num(get(hObj,'String'))); %#ok<ST2NM> str2num distinguishes NaN and invalid input, str2double does not.
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

function popup_handle = createColormapPopup(popup_handle,maps,names,cm_width)
    % Creates one entry in the popup for each colormap.
    %
    narginchk(2,4);
    if nargin<3 || isempty(names)
        names = {};
    end
    if nargin<4 || isempty(cm_width)
        cm_width = 1.1;                                                 % the 1.2 is deterimned empirical. It is roughly the ratio between the width of a space and an x.
    end
    if isinteger(cm_width)
        cm_width = double(cm_width);
    else
        temp_units = popup_handle.Units;
        popup_handle.Units = 'characters';
        cm_width = round(cm_width*popup_handle.Position(3));
        popup_handle.Units = temp_units;
    end
    popup_values = {};
    for imap = 1:numel(maps)
        icol = getColormap(maps{imap});

        short_map = interp1(1:size(icol,1),icol,linspace(1,size(icol,1),cm_width));
        short_map_hex = cellfun(@(x)dec2hex(min(255,max(0,(round(x*255)))),2)',num2cell(short_map,2),'UniformOutput',false);

        iHTML = cellfun(@(hex)['<FONT bgcolor="#' hex(:)' '">&nbsp;</FONT>'],short_map_hex,'UniformOutput',false);
        if ~isempty(names)
            iHTML = ['<HTML><nobr>',iHTML{:},'&nbsp;' names{imap} '</nobr></HTML>'];        
        else
            iHTML = ['<HTML><nobr>',iHTML{:},'</nobr></HTML>'];
        end

        popup_values = [popup_values; iHTML];
    end
    popup_handle.String = popup_values;
end
function outmap = getColormap(cmap,flipFlag)
% This function always returns a colormap as [n 3] matrix. The input can be
% a martix, colormap name or function handle. [n 1] maps are expanded.
    if isnumeric(cmap) && ismatrix(cmap)                   % is matrix of colors
        outmap = cmap;
    elseif isa(cmap,'function_handle')                     % is map stored in a function
        outmap = cmap();
    else                                                   % is the name of a built-in map
        outmap = feval(cmap);
    end
    if size(outmap,2) == 1
        outmap = outmap .* ones(1,3);
    end
    if nargin>1 && flipFlag
        outmap = flip(outmap,1);
    end
end

function createTabs(tabplaceholder,tabcontainer)
    % Creates a tabbed UI at the position of tabplaceholder. All panels
    % contained in tabcontainer are moved as single tabs to the new
    % tabgroup. The order of the tabs is controlled by the tab (key) order.
    temp_units = get(tabplaceholder,'Units');
    set(tabplaceholder,'Units','normalized'); % uitabgroup only supports normalized
    tabgroup = uitabgroup(get(tabplaceholder,'Parent'),'Position',get(tabplaceholder,'Position'),'Units','normalized','Tag',get(tabplaceholder,'Tag'));
    delete(tabplaceholder);
    set(tabgroup,'Units',temp_units);
    
    panelToProcess = flip(findobj(get(tabcontainer,'Children')','-depth',0,'Type','UIPanel'));
    for pidx = 1:numel(panelToProcess)
        currenttab = uitab(tabgroup,'Title',get(panelToProcess(pidx),'Title'),'Tag',get(panelToProcess(pidx),'Tag'));
        set(get(panelToProcess(pidx),'Children'),'Parent',currenttab);
    end
    delete(tabcontainer);
end
