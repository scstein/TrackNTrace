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
