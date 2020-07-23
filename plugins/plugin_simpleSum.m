function [plugin] = plugin_simpleSum()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Simple Sum';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
type = 2;

% The functions this plugin implements
mainFunc =  @sumup;

% Description of output parameters
outParamDescription = {'x','y','z','sum'}; % set depending on plugin options in init function

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Sum up intensity in a window around the candidate positions.'];

% Whether to use parallel processing or not
plugin.useParallelProcessing = true;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('halfwX',...
    'int',...
    {2, 0,inf},...
    'Half window size in X-direction. Window edge size will be 2*halfw+1.');

plugin.add_param('halfwY',...
    'int',...
    {2, 0,inf},...
    'Half window size in Y-direction. Window edge size will be 2*halfw+1.');
end


%   -------------- User functions --------------

function [refinementData] = sumup(img,candidatePos,options,currentFrame)
% Wrapper function for psfFit_Image function (see below). Refer to
% tooltips above and to psfFit_Image help to obtain information on input
% and output variables.
%
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidatePos: 2D double row array of localization candidates created
%     by locateParticles.m. Refer to that function or to TrackNTrace manual
%     for more information.
%
%     options: Struct of input parameters provided by GUI.
%
% OUTPUT:
%     refinementData: 1x1 cell of 2D double array of fitted parameters
%     [x,y,z,A,B,[other parameters]]. Other parameters can be q1, q2, q3
%     (refer to locateParticles.m or to TrackNTrace manual for more
%     information). q_i will be calculated back to sigma_x,sigma_y,
%     rotation angle and possibly z in post-processing function (see
%     below).

refinementData = zeros(size(candidatePos,1),4);
for iCand = 1:size(candidatePos,1)
    % Get position of candidate
    x = candidatePos(iCand,1);
    y = candidatePos(iCand,2);
    
    % Integrate intensity
    yWin = max(1,(y-options.halfwY)):min((y+options.halfwY),size(img,1));
    xWin = max(1,(x-options.halfwX)):min((x+options.halfwX),size(img,2));
    intensity = sum(sum(img(yWin,xWin)));
    
    % Save data
    refinementData(iCand,1:2) = candidatePos(iCand,1:2);    
    refinementData(iCand,4) = intensity;
end

end
