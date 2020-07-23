function [postprocData, postprocOptions] = postprocTracks(trackingData,globalOptions,importOptions,postprocOptions)
% [postprocData, postprocOptions] = postprocTracks(trackingData,globalOptions,postprocOptions)
% This function handles tracks within movies obtained by the
% TrackNTrace routine "trackParticles" and passes them to a postprocessing
% plugin (e.g. for lifetime fitting).
%
% For details on input parameters or fitting algorithm, refer to the
% respective postprocessing plugin.
%
% INPUT:
%     trackingData: Matrix with trackingData. Refer to
%       trackParticles or TrackNTrace manual for details. Uses
%       id,frame,xpos,ypos and sigma. (Not all refinement plugins)
%       [id,frame,xpos,ypos,zpos,Amp,BG,sigma]
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
%
% Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de, 2019

fprintf('TNT: Postprocessing tracks using ''%s''. \n',postprocOptions.plugin_name);

totalTime_start = tic;

% Execute the initializing function
tic;
if ~isempty(postprocOptions.initFunc)
    postprocOptions = postprocOptions.initFunc(postprocOptions,globalOptions,importOptions);
end
initTime = toc;

% Execute the tracking function
tic;
switch nargout(postprocOptions.mainFunc)
    case 1 % Only tracking data assigned during call
        postprocData = postprocOptions.mainFunc(trackingData,postprocOptions);
    case 2 % Tracking data assigned + options changed
        [postprocData, postprocOptions] = postprocOptions.mainFunc(trackingData,postprocOptions);
    otherwise
        warning('Too many output arguments in mainFunc of postprocessing plugin ''%s''. Ignoring additional outputs', trackingOptions.plugin_name);
end
mainTime = toc;

% Execute the post-processing function
tic;
if ~isempty(postprocOptions.postFunc)
    [postprocData, postprocOptions] = postprocOptions.postFunc(postprocData,postprocOptions);
end
postTime = toc;

totalTime = toc(totalTime_start);
fprintf('TNT: Postprocessing took %im %is (init: %im %is, main: %im %is, post: %im %is).\n', floor(totalTime/60), floor(mod(totalTime,60)),floor(initTime/60), floor(mod(initTime,60)) ,floor(mainTime/60), floor(mod(mainTime,60)) ,floor(postTime/60), floor(mod(postTime,60)));

% Verify the outParamDescription, make it fit to the data if neccessary
postprocOptions = verifyOutParamDescription(postprocData,postprocOptions);

end
