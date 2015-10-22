function [ tracks ] = nn_tracker_cpp(  Localizations, varargin)
% Performs gap-closing nearest neighbor tracking of 2d point data recorded over time.
% Gap closing connects endpoints of tracks to the next startpoint of a
% track popping up within the specified distance.

% SYNTAX [ drift ] = drift_calculation_cpp(  Localizations, min_track_length, max_linking_distance)
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

