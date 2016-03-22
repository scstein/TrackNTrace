% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
%     Copyright (C) 2016  Simon Christoph Stein, scstein@phys.uni-goettingen.de
%     Copyright (C) 2016  Jan Thiart, jthiart@phys.uni-goettingen.de
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
function [trackingData, trackingOptions] = trackParticles(fittingData,trackingOptions)
% trajData = trackParticles(fittingData,trackingOptions)
% This function handles particle positions within movies obtained by the
% TrackNTrace routine "locateParticles" and passes them to a particle
% tracking plugin to obtain trajectories.
%
% For details on input parameters or tracking algorithm, refer to the
% respective tracking plugin.
%
% INPUT:
%     fittingData: Cell array of particle positions to track. Refer to
%     locateParticles or TrackNTrace manual for details.
%
%     trackingOptions: Struct of input paramter options for tracking
%     plugins.
%
% OUTPUT:
%     trackingData: 2D double array, list of trajectories with columns 
%     [id,frame,xpos,ypos,zpos] + additional columns. Every trajectory is
%     given an id, starting at 1, after which the list is sorted. 
%     Frame number starts with 1.
% 
%     trackingOptions: see above

fprintf('Tracking particles using %s .. \n',trackingOptions.plugin_name);

% Execute the initializing function
if ~isempty(trackingOptions.initFunc)
    trackingOptions = trackingOptions.initFunc(trackingOptions);
end

% Execute the tracking function
switch nargout(trackingOptions.mainFunc)
    case 1 % Only tracking data assigned during call
        trackingData = trackingOptions.mainFunc(fittingData,trackingOptions);
    case 2 % Tracking data assigned + options changed
        [trackingData, trackingOptions] = trackingOptions.mainFunc(fittingData,trackingOptions);
    otherwise
        warning('Too many output arguments in mainFunc of tracking plugin ''%s''. Ignoring additional outputs', trackingOptions.plugin_name);
end

% Execute the post-processing function
if ~isempty(trackingOptions.postFunc)
    [trackingData,trackingOptions] = trackingOptions.postFunc(trackingData,trackingOptions);
end

fprintf('\b done\n');

% Verify the outParamDescription, make it fit to the data if neccessary
trackingOptions = verifyOutParamDescription(trackingData, trackingOptions);

end
