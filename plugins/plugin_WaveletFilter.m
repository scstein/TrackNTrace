function [plugin] = plugin_WaveletFilter()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Wavelet filtering';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 1;

% The functions this plugin implements
mainFunc =  @findCandidates_wavelet;

% Description of output parameters
outParamDescription = {'x';'y'};

% Create the plugin
plugin = TNTplugin(name,type, mainFunc,outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'This filter implemets wavelet filtering with B-splines. See Izeddin et al, OptExp 20(3),2012 for details.';
plugin.initFunc = @findCandidates_prepareFilter;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('scalingFactor',...
    'float',...
    {2, 0, inf},...
    'Wavelet scaling factor. This determines how fast the filter falls of at the edges. Lower = greater falloff.');
plugin.add_param('splineOrder',...
    'int',...
    {3,1,inf},...
    'B-spline order. Determines how smooth the filter falls off towards the edge. Higher = smoother falloff.');
plugin.add_param('detectionRadius',...
    'int',...
    {3,1,inf},...
    'Size of local maximum detection window. Should be similar to kernelSize.');
plugin.add_param('detectionThreshold',...
    'float',...
    {1,0,inf},...
    'Detection threshold. Higher means less detected candidates.');
plugin.add_param('automaticThreshold',...
    'bool',...
    true,...
    'If set to true, the standard deviation of the wavelet filtered image at wavelet level 1 times the detectionThreshold is taken as the absolute peak intensity threshold. In this case, setting detectionThreshold from 0.5 to 2 is recommended.');
end


function [candidateData] = findCandidates_wavelet(img,options,currentFrame)
% Find candidates in 2D image by wavelet filtering, image dilation, and
% thresholding (in that order). Currently, only wavelet level 2 is
% supported.
% 
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidateOptions: Struct of input parameters provided by GUI.
%
% OUTPUT:
%     candidatePos - 2D array of Nx2 matrix of particle candidate positions
%     [column pixel, row pixel] without subpixel position. Middle of upper
%     left pixel would be [1,1].

% filter image up to level 2
img_filtered = filterImageGeneral(img,options.kernel{1},options.isSeparable,'same');

if options.automaticThreshold
    padding_filtered_img = numel(options.kernel{1})/2;
    threshold_image = img-img_filtered;
    threshold_image = threshold_image(ceil(padding_filtered_img):end-floor(padding_filtered_img),ceil(padding_filtered_img):end-floor(padding_filtered_img)); % can use this for thresholding.
    threshold = std(threshold_image(:));
else
    threshold = 1;
end

% filtered image at level 2 = filter_1(img)-filter_2(img), therefore,
% "convolve" with 1 which is, funnily, faster.
img_filtered = filterImageGeneral(img_filtered,{1,options.kernel{2}},options.isSeparable,'valid');

% normalize
% img_filtered = (img_filtered-min(img_filtered(:)))/(max(img_filtered(:))-min(img_filtered(:)));

% dilate - uses Image Processing Toolbox.
dilated_mask = imdilate(img_filtered,ones(2*options.detectionRadius+1))==img_filtered;
dilated_mask = dilated_mask & img_filtered>(options.detectionThreshold*threshold);
locmax_idx = find(dilated_mask);
[rowIdx,colIdx] = ind2sub(size(img_filtered),locmax_idx);

candidateData = [colIdx,rowIdx]+(ceil(options.kernelSize/2)-1);

end


function [candidateOptions] = findCandidates_prepareFilter(candidateOptions)
% Function calculates Wavelet filters up to wavelet level 2 using
% normalized B-splines according to: 
% Izeddin, I., Boulanger, J., Racine, V., Specht, C. G., Kechkar, A., Nair,
% D., Triller, A., Choquet, D., Dahan, M., and Sibarita, J. B. Wavelet
% analysis for single molecule localization microscopy. Optics Express
% 20(3), 2081–95 (2012).

wavelet_level = 2;

s = candidateOptions.scalingFactor;
q = candidateOptions.splineOrder;
kernelSize = 2*ceil(q*s/2)-1;
kernel = calculateBSpline(((1:kernelSize).'-(kernelSize+1)/2)/s+q/2,q); 
kernel = kernel/sum(kernel(:));

kernel_levels = cell(wavelet_level,1);
kernel_levels{1} = kernel;
for i=2:wavelet_level
    kernel_levels_current = zeros(2*numel(kernel_levels{i-1})-1,1);
    kernel_levels_current(1:2:end) = kernel_levels{i-1};
    kernel_levels{i} = kernel_levels_current;
end

candidateOptions.kernel = kernel_levels;
candidateOptions.isSeparable = true;
candidateOptions.kernelSize = max(cellfun(@(var) numel(var),kernel_levels));

    function b = calculateBSpline(t,q)
        if not(mod(q,1)==0 && q>0)
            error('q-index must be a natural number greater 0.');
        end
        
        b = zeros(numel(t),1);
        
        for iT = 1:numel(t)
            b(iT) = bspline(t(iT),q);
        end
        
        function b = bspline(t,q)
            if q==1
                if t>=0 && t<1
                    b = 1;
                else
                    b = 0;
                end
            else
                b = t/(q-1)*bspline(t,q-1)+(q-t)/(q-1)*bspline(t-1,q-1);
            end
        end
        
    end %calculateBSpline

end

