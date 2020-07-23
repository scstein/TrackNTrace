function [plugin] = plugin_importCSV_can()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Import CSV';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
type = 1;

% The functions this plugin implements
mainFunc =  @convertData;

% Description of output parameters
outParamDescription = {'x';'y'};

% Create the plugin
plugin = TNTplugin(name,type, mainFunc,outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Import candidateData from a CSV file (eg. as exported by ThunderSTORM). \n\n',...
               'Make sure that the unit of x, y, and sigma is px and that the first line contains lables including x and y.'];
 
plugin.initFunc = @importData;
plugin.postFunc = @cleanupOptions;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('CSVfile',...
    'filechooser',...
    {'','csv'},...
    'CSV file to import.');
plugin.add_param('sigma',...
    'bool',...
    true,...
    'If enabled sigma is imported.');
plugin.add_param('otherValues',...
    'bool',...
    true,...
    ['If enabled all any aditional values are imported. Can interfer with the refinement plugin. \\n'...
     'If you want to import tracks enable ''sigma'' and ''otherValues'' and use ''Use candidate data'' and ''Use refinement data''.\\n'...
     'Setting the first frame to something other than one shifts the localisations.']);

end


function candidateData = convertData(~,options,currentFrame)
% 
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidateOptions: Struct of input parameters provided by GUI.
%
% OUTPUT:
%     candidateData
    if isempty(options.csvData)
        candidateData = [];
    else
        ind = options.csvData(:,1)==currentFrame;
        candidateData = options.csvData(ind,2:end);
    end
end


function [candidateOptions] = importData(candidateOptions)

    try
        wstate = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
        datatab = readtable(candidateOptions.CSVfile,'ReadVariableNames',true);
        warning(wstate);
    catch
        warning('TNT: %s: Cannot read csv file.',candidateOptions.plugin_name);
        
        candidateOptions.csvData = [];
        return
    end
    
    containsIsolated = @(str,x)~cellfun(@isempty,regexpi(str, ['(?<![a-zA-Z0-9])' x '(?![a-zA-Z0-9])'],'forcecelloutput'));

    xpos = find(containsIsolated(datatab.Properties.VariableNames,'x'));
    ypos = find(containsIsolated(datatab.Properties.VariableNames,'y'));
    fpos = find(containsIsolated(datatab.Properties.VariableNames,'frame'));
    spos = find(containsIsolated(datatab.Properties.VariableNames,'sigma'));
    otherpos = setdiff(1:numel(datatab.Properties.VariableNames),[xpos, ypos, fpos, spos]);
    
    data = datatab.Variables;
    
    
    if isempty(xpos) || isempty(ypos)
        warning('TNT: %s: Cannot identify x and y colum. Make sure the header exists and contains x and y.',candidateOptions.plugin_name);
        
        candidateOptions.csvData = [];
        return
    end
    outData = [data(:,xpos),data(:,ypos)];
    
    candidateOptions.outParamDescription = {datatab.Properties.VariableNames{xpos};datatab.Properties.VariableNames{ypos}};
    
    if candidateOptions.sigma
        if isempty(spos)
            fprintf('TNT: %s: Cannot identify sigma colum. Set to zero.',candidateOptions.plugin_name);
            outdata = [outData, zeros(size(outData,1),1)];
            candidateOptions.outParamDescription = [candidateOptions.outParamDescription; {'sigma'}];
        else
            outData = [outData,data(:,spos)];
            candidateOptions.outParamDescription = [candidateOptions.outParamDescription; datatab.Properties.VariableNames{spos}];
        end
    end
    
    if candidateOptions.otherValues && ~isempty(otherpos)
        outData = [outData,data(:,otherpos)];
        candidateOptions.outParamDescription = [candidateOptions.outParamDescription; datatab.Properties.VariableNames(otherpos)'];
    end
    
    try % This is just cosmetics and requires R2016b. 
        candidateOptions.outParamDescription = strip(candidateOptions.outParamDescription,'_');
    catch
        
    end
    % the frame is removed by convertData
    if isempty(fpos)
        fprintf('TNT: %s: Cannot identify frame colum. All localisations are assigned to the first frame.',candidateOptions.plugin_name);
        outData = [ones(size(outData,1),1), outData];
    else
        outData = [data(:,fpos), outData];
    end
    candidateOptions.csvData = outData;
end

% Cleanup function. Removes the csvData from the options.
function [candidateData,candidateOptions] = cleanupOptions(candidateData,candidateOptions)
    candidateOptions = rmfield(candidateOptions,'csvData');
end
