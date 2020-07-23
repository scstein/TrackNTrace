% Returns the index of the first import plugin which can handle the
% current filetype.
% Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de, 2020
function pluginIndex = selectImportPlugin(filename,plugins)
    for iPlug = 1:numel(plugins)
        % Take the first plugin which accepts the extension.
        if any(cellfun(@(m_end)~isempty(m_end)&&m_end==numel(filename),cellfun(@(ext)regexpi(filename, regexptranslate('wildcard',ext),'end'),strsplit(plugins(iPlug).info.supportedFormats{:,1},';'),'UniformOutput',false)))
            if isfield(plugins(iPlug).info,'setFile') && isa(plugins(iPlug).info.setFile,'function_handle')
                % If the function setFile is defined call it with the file
                % name. If it has output arguments, a false indicates a
                % incompatible file.
                if nargout(plugins(iPlug).info.setFile)>0
                    validFile = plugins(iPlug).info.setFile(filename);
                    if ~isempty(validFile)&&isscalar(validFile)&&~validFile
                        continue;
                    end
                else
                    plugins(iPlug).info.setFile(filename);
                end
            end
            pluginIndex = iPlug;
            return;
        end
    end
    pluginIndex = 0;
    warning('No compatible import plugin found.');
end