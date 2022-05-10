function [plugin] = plugin_loadPTU()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Import TIF';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
% 4: Postprocessing
% 5: Import
% 6: Export
type = 5;

% The functions this plugin implements
mainFunc =  @read_tiff_helper;

% Description of output parameters
outParamDescription = {'movie'}; % set in init function

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Function to execute before frame-by-frame processing starts
% plugin.initFunc = @updateOutParamDescription;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = struct(...
    'Description','Load a TIF file.',...
    'supportedFormats',{{'*.tif;*.tiff','TIF'}},...
    'hasFLIM',false,...
    'hasTCSPC',false);
end


%   -------------- User functions --------------

function [movie,metadata] = read_tiff_helper(~,filename_movie, frame_range, frame_binning,varargin)
    metadata =  []; % Empty atm. Good place to add resolution information when available.
    movie = read_tiff(filename_movie, false, frame_range);
    % Bin movie if requested.
    if nargin>=4 && ~isempty(frame_binning) && frame_binning > 1
        frame_binning = min(frame_binning,size(movie,3));
        movie(:,:,(end+1):(frame_binning*ceil(end/frame_binning))) = 0; % Pad the 3rd dimensions with zeros to ensure a multiple of the binning
        movie = permute(sum(reshape(movie,size(movie,1),size(movie,2),frame_binning,size(movie,3)/frame_binning),3),[1 2 4 3]);
    end
    movie = {movie};
    if nargout>1
        metadata = struct(...
            'filename',filename_movie...
            );
        tiff_info = imfinfo(filename_movie);
        tiff_info = tiff_info(1); % Take only the first frame
        % It is assumed that XResolution is the true resolution of the
        % whole movie and that the unit is Centimeter (not inches).
        if isfield(tiff_info,'XResolution') && isfield(tiff_info,'ResolutionUnit') &&...
                ~isempty(tiff_info.XResolution) && strcmpi(tiff_info.ResolutionUnit,'Centimeter')
            metadata.pixelsize = 1e7/tiff_info.XResolution;%pixel per centimeter to pixel size in nm
            metadata.pixelsize_unit = 'nm';
        end
    end
end

