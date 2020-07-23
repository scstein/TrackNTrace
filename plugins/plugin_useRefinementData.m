function [plugin] = plugin_useRefinementData()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Use refinement data';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
type = 3;

% The functions this plugin implements
mainFunc =  @convert_refinementData;

% Description of output parameters
outParamDescription = {'Track-ID';'Frame';'x';'y';'z';'Amplitude'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Function to execute before frame-by-frame processing starts
plugin.initFunc = @updateOutParamDescription;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Simply copies the positions acquired by the refinement plugin and makes each one track.\n'...
               'If the data contains the parameter track or ID it is used. (Relvant for imported tracks)'];

end


%   -------------- User functions --------------

function trackingOptions = updateOutParamDescription(trackingOptions)
    global refinementOptions % Needed as we drag along names specified here
    
    containsIsolated = @(str,x)~cellfun(@isempty,regexpi(str, ['(?<![a-zA-Z0-9])' x '(?![a-zA-Z0-9])'],'forcecelloutput'));
    trackingOptions.trackID_pos = containsIsolated(refinementOptions.outParamDescription,'track')|containsIsolated(refinementOptions.outParamDescription,'id');
    if any(trackingOptions.trackID_pos)
        trackingOptions.outParamDescription = ['Track-ID';'Frame';refinementOptions.outParamDescription(~trackingOptions.trackID_pos)];
    else
        trackingOptions.outParamDescription = ['Track-ID';'Frame';refinementOptions.outParamDescription];
    end
end

function [trackingData,trackingOptions] = convert_refinementData(refinementData,trackingOptions)
    trackingData = [];
    for iFrame=1:size(refinementData,1)
        if ~isempty(refinementData{iFrame})
            trackingData = [trackingData; [iFrame*ones(size(refinementData{iFrame},1),2), refinementData{iFrame}]];
        end
    end
    if any(trackingOptions.trackID_pos)
        pos = find(trackingOptions.trackID_pos)+2;
        trackingData = [trackingData(:,pos), trackingData(:,[(2:(pos-1)), ((pos+1):end)])];
        trackingData = sortrows(trackingData,1);
    else
        trackingData(:,1) = 1:size(trackingData,1);
    end
    trackingOptions = rmfield(trackingOptions,'trackID_pos');
end

