function [plugin] = plugin_FilterWithLocmax()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Image filtering';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 1;

% The functions this plugin implements
mainFunc =  @findCandidates_filter;

% Description of output parameters
outParamDescription = {'x';'y'};

% Create the plugin
plugin = TNTplugin(name,type, mainFunc,outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Filter image with appropriate kernel and locate features by local maximum suppression.';
plugin.initFunc = @findCandidates_prepareFilters;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('kernelName',...
    'list',...
    {'Average','Average Difference','Gauss','Gauss Zero','Gauss Difference','Median'},...
    'Choose filter kernel.');
plugin.add_param('PSFSigma',...
    'float',...
    {1.3, 0, inf},...
    'PSF standard deviation in [pixel]. FWHM = 2*sqrt(2*log(2))*sigma. Applies to Gaussian filters and overwrites kernelSize.');
plugin.add_param('kernelSize',...
    'int',...
    {5,3,inf},...
    'Window size of first kernel. Does not apply to Gaussian filters.');
plugin.add_param('filterBlowup',...
    'float',...
    {1.5,1.01,inf},...
    'Blow up PSFSigma or kernelSize by this factor when applying second filtering step. Only applies to Average Difference and Gauss Difference.');
plugin.add_param('detectionRadius',...
    'int',...
    {3,1,inf},...
    'Size of local maximum detection window. Should be similar to kernelSize.');
plugin.add_param('detectionThreshold',...
    'float',...
    {0.5,0,1},...
    'Detection threshold. Higher means less detected candidates.');
end


function [candidateData] = findCandidates_filter(img,options,currentFrame)
% Find candidates in 2D image by filtering, image dilation, and
% thresholding (in that order). Possible filters are Average, Difference of
% averages, Gauss, Difference of Gauss, Zeroed Gauss and Median.
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

% filter image
if iscell(options.kernel) || ~isnan(options.kernel)
    img_filtered = filterImageGeneral(img,options.kernel,options.isSeparable,'valid');
else
    img_filtered = medfilt2(img,[options.kernelSize,options.kernelSize]);
    img_filtered = img_filtered(ceil(options.kernelSize/2):end-floor(options.kernelSize/2),ceil(options.kernelSize/2):end-floor(options.kernelSize/2));
end
% normalize
img_filtered = (img_filtered-min(img_filtered(:)))/(max(img_filtered(:))-min(img_filtered(:)));

% dilate - uses Image Processing Toolbox.
dilated_mask = imdilate(img_filtered,ones(2*options.detectionRadius+1))==img_filtered;
dilated_mask = dilated_mask & img_filtered>options.detectionThreshold;
locmax_idx = find(dilated_mask);
[rowIdx,colIdx] = ind2sub(size(img_filtered),locmax_idx);

candidateData = [colIdx,rowIdx]+(ceil(options.kernelSize/2)-1);

end


function [candidateOptions] = findCandidates_prepareFilters(candidateOptions)

PSFSigma = candidateOptions.PSFSigma;
kernelSize = candidateOptions.kernelSize;
sizeFactor = candidateOptions.filterBlowup;
kernelSizeLarge = kernelSize;

switch candidateOptions.kernelName
    case 'Average'
        kernel1D = ones(kernelSize,1)/kernelSize;
%         kernel2D = kernel1D*kernel1D.';
%         [kernel,isSeparable] = testFilterPerformance(kernel1D,kernel2D);
        
    case 'Average Difference'
        kernelSizeLarge = round(kernelSize*sizeFactor)+1*(round(kernelSize*sizeFactor)==kernelSize);
        kernel1D = {ones(kernelSize,1)/kernelSize,ones(kernelSizeLarge,1)/kernelSizeLarge};
%         kernel2D = {ones(kernelSize)/kernelSize^2,ones(kernelSizeLarge)/kernelSizeLarge^2};
%         [kernel,isSeparable] = testFilterPerformance(kernel1D,kernel2D);
        
    case 'Gauss'
        kernelSize = 2*ceil(3*PSFSigma)+1;
        kernel1D = gauss_window(kernelSize,PSFSigma);
%         kernel2D = kernel1D*kernel1D.';
%         [kernel,isSeparable] = testFilterPerformance(kernel1D,kernel2D);
        
    case 'Gauss Zero'
        kernelSize = 2*ceil(3*PSFSigma)+1;
        kernel1D = {gauss_window(kernelSize,PSFSigma),ones(kernelSize,1)/kernelSize};
%         kernel2D = kernel1D{1}*kernel1D{1}.'; kernel2D = kernel2D-mean(kernel2D(:));
%         [kernel,isSeparable] = testFilterPerformance(kernel1D,kernel2D);
        
    case 'Gauss Difference'
        PSFSigmaLarge = sizeFactor*PSFSigma;
        kernelSize = 2*ceil(3*PSFSigma)+1;
        kernelSizeLarge = 2*ceil(3*PSFSigmaLarge)+1;
        kernel1D = {gauss_window(kernelSize,PSFSigma),gauss_window(kernelSizeLarge,PSFSigmaLarge)};
%         kernel2D = {kernel1D{1}*kernel1D{1}.',kernel1D{2}*kernel1D{2}.'};
%         [kernel,isSeparable] = testFilterPerformance(kernel1D,kernel2D);
        
    case 'Median'
        kernel1D = NaN;
        
    otherwise
        warning off backtrace
        warning('Unknown kernel. Falling back to moving average.');
        warning on backtrace
        kernel1D = ones(kernelSize,1)/kernelSize;
%         kernel2D = kernel1D*kernel1D.';
%         [kernel,isSeparable] = testFilterPerformance(kernel1D,kernel2D);

end
% testing on a variety of kernels and image sizes has determined that 2D
% convolution only makes sense for either a large image and a very small
% kernel or a very small image and a very large kernel, and even then the
% difference is neglible. 
kernel = kernel1D; isSeparable = true;


candidateOptions.kernel = kernel;
candidateOptions.isSeparable = isSeparable;
candidateOptions.kernelSize = max(kernelSizeLarge,kernelSize);

    function [gauss_norm] = gauss_window(w,s)
    gauss_norm = exp(-((1:w).'-(w+1)/2).^2/(2*s^2));
    gauss_norm = gauss_norm/sum(gauss_norm(:));
    end
end


% function [kernel,isSeparable] = testFilterPerformance(kernel1D,kernel2D)
% global movie
% 
% [nrRow,nrCol] = size(movie);
% img_test = rand(nrRow,nrCol);
% 
% if iscell(kernel2D)
%     kernelSize = size(kernel2D{2},1);
% else
%     kernelSize = size(kernel2D,1);
% end
% 
% if min(nrRow,nrCol) <= kernelSize
%     error('Kernel size cannot be larger than image size. Please adjust.');
% end
% 
% % Matlab does some weird optimizations sometimes, so call each function
% % once and "initialize" the memory.
% [img_out1] = filterImageGeneral(img_test,kernel1D,true,'valid');
% [img_out2] = filterImageGeneral(img_test,kernel2D,false,'valid');
% 
% tic;
% for i=1:5
%     [img_out1] = filterImageGeneral(img_test,kernel1D,true,'valid');
% end
% t1 = toc;
% 
% 
% tic;
% for i=1:5
%     [img_out2] = filterImageGeneral(img_test,kernel2D,false,'valid');
% end
% t2 = toc;
% 
% % if 2x1D convolution is faster, pick that
% if t1<t2
%     kernel = kernel1D;
%     isSeparable = true;
% else
%     kernel = kernel2D;
%     isSeparable = false;
% end
% 
% end

