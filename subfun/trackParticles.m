function [trajectoryData, trackingOptions] = trackParticles(fitData,trackingOptions)
% trajData = trackParticles(fitData,trackingOptions)
% This function handles particle positions within movies obtained by the
% TrackNTrace routine "locateParticles" and passes them to a particle
% tracking plugin to obtain trajectories.
%
% For details on input parameters or tracking algorithm, refer to the
% respective tracking plugin.
%
% INPUT:
%     fitData: Cell array of particle positions to track. Refer to
%     locateParticles or TrackNTrace manual for details.
%
%     trackingOptions: Struct of input paramter options for tracking
%     plugins.
%
% OUTPUT:
%     trajectoryData: 2D double array, list of trajectories in a suitable
%     format [id,frame,xpos,ypos,amp]. Every trajectory is given an id,
%     starting at 1, after which the list is sorted. Frame number starts
%     with 1, positions are given in pixels and the amplitude is the
%     gaussian amplitude value A in A*exp(...)+B given by locateParticles.
% 
%     trackingOptions: see above

fprintf('Tracking particles using %s .. \n',trackingOptions.plugin_name);

if ~isempty(trackingOptions.initFunc)
    trackingOptions = trackingOptions.initFunc(trackingOptions);
end

trajectoryData = trackingOptions.mainFunc(fitData,trackingOptions);

if ~isempty(trackingOptions.postFunc)
    trajectoryData = trackingOptions.postFunc(trajectoryData,trackingOptions);
end

fprintf('\b done\n');

end
