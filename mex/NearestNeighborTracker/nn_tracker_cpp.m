function [ tracks ] = nn_tracker_cpp(  Localizations, varargin)
% Performs gap-closing nearest neighbor tracking of 3d point data recorded over time.
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
% C++ implementation is unter the FreeBSD license and uses the
% nanoflann library, which is licensed under the BSD license (see below).
% 

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