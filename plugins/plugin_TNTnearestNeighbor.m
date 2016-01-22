function [plugin] = plugin_TNTnearestNeighbor()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'TNT NearestNeighbor';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 3;

% The functions this plugin implements
mainFunc =  @trackParticles_nearestNeighborCPP;

% Description of output parameters
outParamDescription = {'Placeholder'}; % is set in mainFunc adapting to input data

% Create the plugin
plugin = TNTplugin(name,type, mainFunc, outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Simple and fast nearest neighbor tracking implemented in C++\n\n Implementation uses the nanoflann library written by Marius Muja, David G. Lowe and Jose Luis Blanco. Nanoflann is under the BSD license.';

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('minTrajLength',...
    'int',...
    {2, 0, inf},...
    'Minimum length of trajectories AFTER gap closing in [frames].');
plugin.add_param('maxTrackRadius',...
    'float',...
    {6, 0, inf},...
    'Maximum allowed linking distance between two spots in [pixels].');
plugin.add_param('maxFrameGap',...
    'float',...
    {0, 0, inf},...
    'Maximum allowed time gap between two segments used for gap closing in [frames]. 0 = no gap closing.');
plugin.add_param('minSegLength',...
    'int',...
    {2, 0, inf},...
    'Minimum length of trajectory segments BEFORE gap closing in [frames].');
plugin.add_param('maxGapRadius',...
    'int',...
    {6, 0, inf},...
    'Maximum allowed linking distance between two segments uses for gap closing in [pixels].');
plugin.add_param('verbose',...
    'bool',...
    false,...
    'Switch on to see tracking progress in command window.');
end


%   -------------- User functions --------------

function [trackingData, options] = trackParticles_nearestNeighborCPP(fitData,options)
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
%     [id,frame,xpos,ypos,zpos,amp]. Refer to trackParticles.m or to TrackNTrace
%     manual for more information.


% convert fitData cell array to tracker input format
nrFrames = size(fitData,1);
nrPosInFrame = cellfun('size',fitData,1);

% builds a vector with the frame number for every position
% e.g. [1,1,1,2,2,3,3,3,3] for 3 particles in frame 1, 2 p. in fr. 2, 4 p. in fr. 3 etc.
frameVec = arrayfun( @(val,nr) repmat(val,nr,1), [1:nrFrames].',nrPosInFrame,'UniformOutput',false);
frameVec = vertcat(frameVec{:});

% build a matrix with the frame number attached to each localization
pos = [frameVec, vertcat(fitData{:})].';

% call main function
trackingData = nn_tracker_cpp(pos,options.minSegLength,options.maxTrackRadius,options.maxGapRadius,options.maxFrameGap,options.minTrajLength,options.verbose).';

% Set name of output variables
global fittingOptions % Needed as we drag along names specified here
options.outParamDescription = ['Track-ID'; 'Frame'; fittingOptions.outParamDescription(:)];

end


function [ tracks ] = nn_tracker_cpp(  Localizations, varargin)
% Performs gap-closing nearest neighbor tracking of 2d point data recorded over time.
% Gap closing connects endpoints of tracks to the next startpoint of a
% track popping up within the specified distance.

% SYNTAX [ tracks ] = nn_tracker_cpp( Localizations, min_track_length,max_linking_distance,max_distance_gap,max_frame_gap,min_track_length_afterGapClosing,verbose )
%
% Input:
%     Localizations - (4+n)xN matrix of points, where N is the number of points 
%                     and the rows correspond to (frame,x,y,z + n additional data). 
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
%     tracks - (5+n)xP matrix of points, where P is the number of points in
%              the tracks and the rows correspond to (trackID, frame,x,y,z + n 
%              additional data from input). 
%
%  Example: to get the data for track 1 and plot its y-x trajectory:
%     track1data = tracks(tracks(1,:)==1, :);
%     plot(track1data(3,:), track1data(4,:));
%  To plot the y-t over time movement
%     plot(track1data(2,:), track1data(4,:));

% Make sure logicals are passed as correct datatype
%
% C++ implementation is written using the nanoflann library. Licensing see below.

if numel(varargin) == 6
    varargin{6} = logical(varargin{6});
end

tracks = mx_nn_tracker(Localizations, varargin{:});
end

% /***********************************************************************
 % * Software License Agreement (BSD License) for the nanoflann library
 % *
 % * Copyright 2008-2009  Marius Muja (mariusm@cs.ubc.ca). All rights reserved.
 % * Copyright 2008-2009  David G. Lowe (lowe@cs.ubc.ca). All rights reserved.
 % * Copyright 2011-2013  Jose Luis Blanco (joseluisblancoc@gmail.com).
 % *   All rights reserved.
 % *
 % * THE BSD LICENSE
 % *
 % * Redistribution and use in source and binary forms, with or without
 % * modification, are permitted provided that the following conditions
 % * are met:
 % *
 % * 1. Redistributions of source code must retain the above copyright
 % *    notice, this list of conditions and the following disclaimer.
 % * 2. Redistributions in binary form must reproduce the above copyright
 % *    notice, this list of conditions and the following disclaimer in the
 % *    documentation and/or other materials provided with the distribution.
 % *
 % * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 % * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 % * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 % * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 % * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 % * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 % * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 % * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 % * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 % * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 % *************************************************************************/

