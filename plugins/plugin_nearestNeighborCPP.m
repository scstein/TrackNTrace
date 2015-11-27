function [plugin_name, plugin_type] = plugin_nearestNeighborCPP(h_panel, inputOptions)
if nargin < 2
    inputOptions = [];
end

% Name of the component these options are for
plugin_name = 'NearestNeighbor C++';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
plugin_type = 3;

% Enter names of the parameters
% These translate to the names of variables inside options struct this plugin
% outputs by removing all white spaces.
par_name  = {'minTrajLength','maxTrackRadius','maxFrameGap','minSegLength','maxGapRadius','verbose'};

% Enter type of the parameters
% possible: 'float', 'int', 'bool','list'
par_type  = {'int','float','float','int','int','bool'};

% Default value for parameters
% Should be a number for 'float'/'int', true/false for 'bool'
% or a cell array string list of possible choices for 'list' (first entry is default)
par_defaultValue = {2,6,0,2,6,false};

% Tooltip for the parameters
par_tooltip = {'Minimum length of trajectories AFTER gap closing in [frames].',...
    'Maximum allowed linking distance between two spots in [pixels].',...
    'Maximum allowed time gap between two segments used for gap closing in [frames]. 0 = no gap closing.',...
    'Minimum length of trajectory segments BEFORE gap closing in [frames].',...
    'Maximum allowed linking distance between two segments uses for gap closing in [pixels].',...
    'Switch on to see tracking progress in command window.'};

% Calling the plugin function without arguments just returns its name and type
if (nargin == 0); return; end

% Create the panel for this plugin
createOptionsPanel(h_panel, plugin_name, par_name, par_type, par_defaultValue, par_tooltip,inputOptions);

% Save handle of the plugins function
options = getappdata(h_panel,'options');
options.functionHandle = @trackParticles_nearestNeighborCPP;
setappdata(h_panel,'options',options);
end


%FUNCTION CODE STARTS HERE
function [trajData] = trackParticles_nearestNeighborCPP(fitData,options)
% Wrapper function for nearest neighbor C++ tracker (see below). Refer to
% tooltips above and to nn_tracker_cpp help to obtain information on input
% and output variables.
%
% INPUT:
%     fitData: Cell array of localizations created by locateParticles.m
%     Refer to that function or to TrackNTrace manual for more information.
%
%     options: Struct of input parameters provided by GUI.
%
% OUTPUT:
%     trajData: 2D double array of trajectories in format
%     [id,frame,xpos,ypos,amp]. Refer to trackParticles.m or to TrackNTrace
%     manual for more information.


% convert fitData cell array to adjust to tracker function input
nrFrames = size(fitData,1);
nrPos = zeros(nrFrames,1);
for iFrame=1:nrFrames
    if(isempty(fitData{iFrame})); continue; end; % Jump empty frames
    nrPos(iFrame) = sum((fitData{iFrame}(:,6)==1)); %error flag is 1?
end
pos = zeros(6,sum(nrPos));

nrPos_cs = [0;cumsum(nrPos)];

for iFrame = 1:nrFrames
    if(isempty(fitData{iFrame})); continue; end; % Jump empty frames
    pos_frame_now = fitData{iFrame}.';
    valid_pos = pos_frame_now(6,:)==1; %error flag is 1?
    pos_frame_now(pos_frame_now==0) = 1e-6; %this is a dirty hack for particles which run out of the frame
    pos(:,nrPos_cs(iFrame)+1:nrPos_cs(iFrame+1)) = [repmat(iFrame,1,nrPos(iFrame));pos_frame_now(1:5,valid_pos)];
end


% call main function
trajData = nn_tracker_cpp(pos,options.minSegLength,options.maxTrackRadius,options.maxGapRadius,options.maxFrameGap,options.minTrajLength,options.verbose).';
end


function [ tracks ] = nn_tracker_cpp(  Localizations, varargin)
% Performs gap-closing nearest neighbor tracking of 2d point data recorded over time.
% Gap closing connects endpoints of tracks to the next startpoint of a
% track popping up within the specified distance.

% SYNTAX [ tracks ] = nn_tracker_cpp(  Localizations, min_tak-Length,max_linking_distance,max_distance_gap,max_frame_gap,min_track_length_afterGapClosing,verbose )
%
% Input:
%     Localizations - (3+n)xN matrix of points, where N is the number of points
%                     and the rows correspond to (frame,x,y, + n additional data).
%                     Data must be sorted with increasing frame number.
%     min_track_length - Tracks with less then min_track_length connected
%                        positions are removed before gap closing. | default: 0
%     max_linking_distance - Points farther away then this distance will
%                            not be linked. | default: inf
%     max_distance_gap - Max distance to close gaps over | default: 0
%     max_frame_gap - Max distance in time to connect gaps | default: 0
%     min_track_length_afterGapClosing - Tracks shorter than this length
%                                        are rejected after gap closing | default: 0
%     verbose: If true, information is printed to the console | default: false
%
%  Parameters can be left empty [] to use their default values.
%
% Output:
%     tracks - (4+n)xP matrix of points, where P is the number of points in
%              the tracks and the rows correspond to (trackID, frame,x,y, + n
%              additional data from input).
%
%  Example: to get the data for track 1 and plot its y-x trajectory:
%     track1data = tracks(tracks(1,:)==1, :);
%     plot(track1data(3,:), track1data(4,:));
%  To plot the y-t over time movement
%     plot(track1data(2,:), track1data(4,:));

% Make sure logicals are passed as correct datatype
if numel(varargin) == 6
    varargin{6} = logical(varargin{6});
end

tracks = mx_nn_tracker(Localizations, varargin{:});
end

