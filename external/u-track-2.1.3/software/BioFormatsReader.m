classdef  BioFormatsReader < Reader
    % BioFormatsReader is a Reader subclass which reads metadata/pixels
    % from image files using the Bio-Formats library
%
% Copyright (C) 2014 LCCB 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
    
    properties (SetAccess=protected, Transient=true)
        id
        formatReader
        series
    end
    
    methods
        %% Constructor
        function obj = BioFormatsReader(varargin)
            
            % Input check
            ip = inputParser();
            ip.addRequired('id', @ischar);
            ip.addOptional('series', 0, @(x) isscalar(x) && isnumeric(x));
            ip.addParamValue('reader', [], @(x) isa(x, 'loci.formats.IFormatReader'));
            ip.parse(varargin{:});
            
            % Initialize Bio-Formats and log4j
            bfInitLogging();
            
            obj.id = ip.Results.id;
            if ~isempty(ip.Results.reader),
                obj.formatReader = ip.Results.reader;
            else
                obj.formatReader = bfGetReader(obj.id, false);
            end
            obj.series = ip.Results.series;
        end
        
        function metadataStore = getMetadataStore(obj)
            metadataStore = obj.getReader().getMetadataStore();
        end
        
        function r = getReader(obj)
            r = obj.formatReader;
            if isempty(r.getCurrentFile())
                r.setId(obj.id);
            end
            r.setSeries(obj.getSeries());
        end
        
        function series = getSeries(obj)
            series = obj.series;
        end
        
        function sizeX = getSizeX(obj)
            sizeX = obj.getMetadataStore().getPixelsSizeX(obj.getSeries()).getValue();
        end
        
        function sizeY = getSizeY(obj)
            sizeY = obj.getMetadataStore().getPixelsSizeY(obj.getSeries()).getValue();
        end
        
        function sizeZ = getSizeZ(obj)
            sizeZ = obj.getMetadataStore().getPixelsSizeZ(obj.getSeries()).getValue();
        end
        
        function sizeT = getSizeT(obj)
            sizeT = obj.getMetadataStore().getPixelsSizeT(obj.getSeries()).getValue();
        end
        
        function sizeC = getSizeC(obj)
            sizeC = obj.getMetadataStore().getPixelsSizeC(obj.getSeries()).getValue();
        end
        
        function bitDepth = getBitDepth(obj)
            pixelType = obj.getReader().getPixelType();
            bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
            bitDepth = 8 * bpp;
        end
        
        function fileNames = getImageFileNames(obj, iChan, iFrame, varargin)
            % Generate image file names
            usedFiles = obj.getReader().getUsedFiles();
            [~, fileName] = fileparts(char(usedFiles(1)));
            if obj.getReader().getSeriesCount() > 1
                basename = sprintf('%s_s%g_c%d_t',fileName, obj.getSeries()+1, iChan);
            else
                basename = sprintf('%s_c%d_t',fileName, iChan);
            end
            
            fileNames = arrayfun(@(t) [basename num2str(t, ['%0' num2str(floor(log10(obj.getSizeT))+1) '.f']) '.tif'],...
                1:obj.getSizeT,'Unif',false);
            if(nargin > 2)
                fileNames = fileNames(iFrame);
            end
        end
        
        function channelNames = getChannelNames(obj, iChan)
            usedFiles = obj.getReader().getUsedFiles();
            [~, fileName, fileExt] = fileparts(char(usedFiles(1)));
            base = [fileName fileExt];
            if obj.getReader().getSeriesCount() > 1
                base = [base ' Series ' num2str(obj.getSeries()+1)];
            end
            base = [base ' Channel '];
            
            channelNames = arrayfun(@(x) [base num2str(x)], iChan, 'Unif',false);
        end
        
        function I = loadImage(obj, c, t, varargin)
            % Retrieve single plane specified by its (c, t, z) coordinates
            
            ip = inputParser;
            ip.addRequired('c', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC()));
            ip.addRequired('t', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT()));
            ip.addOptional('z', 1, @(x) isscalar(x) && ismember(x, 1 : obj.getSizeZ()));
            ip.parse(c, t, varargin{:});
            
            % Using bioformat tools, get the reader and retrieve dimension order
            I = loadImage_(obj, c, t, ip.Results.z);
        end
        
        function delete(obj)
            obj.formatReader.close()
        end
    end
    methods ( Access = protected )
        function I = loadImage_(obj, c, t, z, varargin)
            javaIndex =  obj.getReader().getIndex(z - 1, c - 1, t - 1);
            I = bfGetPlane(obj.getReader(), javaIndex + 1);
        end
    end
end
