function [plugin] = plugin_useCandidateData()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Use candidate data';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 2;

% The functions this plugin implements
mainFunc =  @convert_candidateData;

% Description of output parameters
outParamDescription = {' '}; % set in init function

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Function to execute before frame-by-frame processing starts
plugin.initFunc = @updateOutParamDescription;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Copies the positions acquired by the candidate detection plugin.\n\n'];

end


%   -------------- User functions --------------

function fittingOptions = updateOutParamDescription(fittingOptions)
global candidateOptions % Needed as we drag along names specified here
if(numel(candidateOptions.outParamDescription)>2)
    fittingOptions.outParamDescription = [candidateOptions.outParamDescription(1:2); 'z'; candidateOptions.outParamDescription(3:end)];
else
    fittingOptions.outParamDescription = [candidateOptions.outParamDescription(1:2); 'z'];
end
end


function [fitData] = convert_candidateData(img,candidatePos,options,currentFrame)
% We extend the candidate data with the z position column to adhere to TNT specifications of fitting data.
% All additional columns are simply copied.
if(size(candidatePos,2)>2)
    fitData = [candidatePos(:,1:2),zeros(size(candidatePos,1),1),candidatePos(:,3:end)]; %adding z = 0
else
    fitData = [candidatePos(:,1:2),zeros(size(candidatePos,1),1)]; %adding z = 0
end

end

