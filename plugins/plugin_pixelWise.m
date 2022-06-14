function [plugin] = plugin_pixelWise()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Pixel wise grid';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
type = 1;

% The functions this plugin implements
mainFunc =  @pxGrid;

% Description of output parameters
outParamDescription = {'x';'y';'z';'Amp';'sigma'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Places a localization at every pixel position. \n\nUseful to quickly test lifetime fitting methods, e.g. on solution measurements.\nUse the "use candidate data" refinement plugin to pass the localizations to tracking.';

% plugin.initFunc = @makeGrid;
% plugin.postFunc = @cleanupOptions;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
plugin.add_param('PSFsigma',...
    'float',...
    {0.5, 0, inf},...
    'Standard deviation of the PSF in pixels. \nsigma = FWHM/(2*sqrt(2*log(2))) ~ 0.21*lambda/NA where lambda is the emission wavelength in pixels and NA is the numerical aperture of the objective.');
plugin.add_param('Spacing X',...
    'float',...
    {1, 1, inf},...
    'Distance between localizations in x in pixels.');
plugin.add_param('Spacing Y',...
    'float',...
    {1, 1, inf},...
    'Distance between localizations in y in pixels.');
end


%   -------------- User functions --------------

function candidatePos = pxGrid(img, options, currentFrame)
%Wrapper function for cross-correlation candidate finding. Refer to
%matchPattern below or tooltips above to obtain information on input and
%output variables.
%
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidateOptions: Struct of input parameters provided by GUI.
% 
%     currentFrame: Integer, current movie frame in main loop.
%
% OUTPUT:
%     candidatePos - 2D array of Nx2 matrix of particle candidate positions
%     [column pixel, row pixel] without subpixel position. Middle of upper
%     left pixel would be [1,1].

sigma = options.PSFsigma;
deltaX = options.Spacing_X;
deltaY = options.Spacing_Y;
imgSZ = size(img);

posX = 1+ceil((deltaX-1)/2):deltaX:imgSZ(2);
posY = 1+ceil((deltaY-1)/2):deltaY:imgSZ(1);

[posX,posY] = meshgrid(posX,posY);

if ceil(sigma)>1
    halfw = ceil(3*sigma);
    X = -halfw:halfw;
    mask = exp(-(X(:).^2+X.^2)./2./sigma.^2);
    img = conv2(img,mask./sum(mask(:)),'same');
end

% ATTENTION! Here we switch the coordinates
% candidatePos(:,1) is columns (x-coordinate)
% candidatePos(:,2) is rows    (y-coordinate)
candidatePos = [posX(:),posY(:),zeros(numel(posX),1),img(sub2ind(imgSZ,posY(:),posX(:))),sigma*ones(numel(posX),1)];

end


