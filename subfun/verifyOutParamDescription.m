function options = verifyOutParamDescription(outputData, options)
% Verifies there is a description for every column of the output data.
% Columns without a descriptor get the description '<Unknown>'. If
% there are too many descriptors, they are simply cut down to the
% number of columns present in the output.

outParamDescription = options.outParamDescription;
unknown_name = '<Unknown>';

% Find number of parameters in candidateData
if iscell(outputData) % if outputData is in cell format
    nrOutParams = 0;
    for iFrame = 1:size(outputData,1)
        if ~isempty(outputData{iFrame})
            nrOutParams = size(outputData{iFrame},2);
            break;
        end
    end
else % if outputData is in matrix format
    nrOutParams = size(outputData,2);
end

% Is there a descriptor for every parameter?
if numel(outParamDescription) == nrOutParams
    return;
end

warning off backtrace;
warning(['The number of descriptors specified in outParamDescription by Plugin ''%s'' does not match the number ', ...
    'of columns in its output. (specified: %i, output: %i). Please update the outParamDescription in the ', ...
    'plugin file by either setting plugin.outParamDescription or setting options.outParamDescription in the ',...
    'initFunc/outFunc if the output depends on the selected options.'],options.plugin_name, numel(outParamDescription), nrOutParams);
warning on backtrace;

% Make outParamDescription the same size as the number of params
% If there are more columns than specified in the plugin, they are
% named <Unknown>. If there are less columns, the description is
% cut down to the number of columns present in the output data
if isempty(outParamDescription)
    outParamDescription = repmat({unknown_name}, nrOutParams,1);
else
    if numel(outParamDescription) > nrOutParams
        outParamDescription = outParamDescription(1:nrOutParams);
    elseif numel(outParamDescription) < nrOutParams
        tmp = outParamDescription;
        outParamDescription = cell(nrOutParams,1);
        outParamDescription(:) = {unknown_name};
        outParamDescription(1:numel(tmp)) = tmp(1:numel(tmp));
    end
end

% Assign to Options
options.outParamDescription = outParamDescription;
end
