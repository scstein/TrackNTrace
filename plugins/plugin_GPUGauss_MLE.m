function [plugin] = plugin_GPUGauss_MLE()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'GPU-Gauss MLE';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 2;

% The functions this plugin implements
mainFunc = @refineParticles_gpugaussmle;

% Description of output parameters
outParamDescription = {'x';'y';'z';'Intens. (integ.)'; 'Background'; 'sigma_x'; 'sigma_y'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Add initial function
plugin.initFunc = @refineParticles_gpugaussinit;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Fit a Gaussian PSF via GPU-MLE fitting. Absolutely requires photon conversion.';

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('PSFSigma',...
    'float',...
    {0, 1.3, inf},...
    'PSF standard deviation in [pixel]. FWHM = 2*sqrt(2*log(2))*sigma.');
plugin.add_param('fitType',...
    'list',...
    {'[x,y,A,BG]', '[x,y,A,BG,s]','[x,y,A,BG,sx,sy]'},...
    'Fit positions (xy), amplitude & background (A,BG), and sigma (s, or sx & sy for elliptic Gaussian).');
plugin.add_param('Iterations',...
    'int',...
    {3, 10, inf},...
    'Number of fit iterations.');

end


function [fitData] = refineParticles_gpugaussmle(img,candidatePos,options,currentFrame)
% Wrapper function for gaussmlev2 function (see below). Refer to tooltips
% above and to gaussmlev2 help to obtain information on input and output
% variables. gaussmlev2.m was released as part of the following
% publication under the GNU public licence: 
% Smith et al, Fast, single-molecule localization that achieves
% theoretically minimum uncertainty, Nature Methods 7, 373-375 (2010),
% doi:10.1038/nmeth.1449
% 
%
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidatePos: 2D double row array of localization candidates created
%     by locateParticles.m. Refer to that function or to TrackNTrace manual
%     for more information.
%
%     options: Struct of input parameters provided by GUI.
%
% OUTPUT:
%     fitData: 1x1 cell of 2D double array of fitted parameters
%     [x,y,z,A,B,[other parameters]]. Other parameters can be PSF standard
%     deviation, both in x and y.

fitData = zeros(size(candidatePos,1),5+options.nrParam);

for iCand = 1:size(candidatePos,1)
    pos = candidatePos(iCand,1:2);
    idx_x = max(1,pos(1)-options.halfWindowSize):min(size(img,2),pos(1)+options.halfWindowSize);
    idx_y = max(1,pos(2)-options.halfWindowSize):min(size(img,1),pos(2)+options.halfWindowSize);
    [param,~,~]=gaussmlev2(single(img(idx_y,idx_x)),options.PSFSigma,options.Iterations,options.fitType,options.functionName); %[x,y,z,A,BG,[rest params]]
    fitData(iCand,:) = [param(1:2)+[idx_x(1),idx_y(1)],0,param(3:end)]; %correct position relative to fit window, middle of first pixel is [0,0]
end

fitData = fitData(fitData(:,1)>0,:);

end


function [fittingOptions] = refineParticles_gpugaussinit(fittingOptions)
global globalOptions

if ~globalOptions.usePhotonConversion
    warning off backtrace
    warning('GPU-Gauss MLE strictly requires photon conversion!');
    warning on backtrace
end

fittingOptions.halfWindowSize = min(10,ceil(3*fittingOptions.PSFSigma));
switch fittingOptions.fitType
    case '[x,y,A,BG]'
        fitType = 1;
        nrParam = 0;
    case '[x,y,A,BG,s]'
        fitType = 2;
        nrParam = 1;
    case '[x,y,A,BG,sx,sy]'
        fitType = 4;
        nrParam = 2;
end

fittingOptions.fitType = fitType;
fittingOptions.nrParam = nrParam;
fittingOptions.outParamDescription = fittingOptions.outParamDescription(1:5+nrParam);

functions = {'gaussmlev2_cuda50','gaussmlev2_cuda42','gaussmlev2_cuda40',...
    'gaussmlev2_c_thread','gaussmlev2_c','gaussmlev2_matlab'};

% Try CUDA availability
try
    evalc(listgpus());
catch
    warning off backtrace
    warning('Could not find CUDA GPU. Verify your CUDA installation if possible.');
    warning on backtrace
    functions = functions(4:end);
end

% Try rest of functions
[x,y] = meshgrid(-10:10,-10:10);
z = single(round(100*exp(-x.^2/2-y.^2/2) + 10*rand(21,21)));

for iFun = 1:numel(functions)
    functionName = functions{iFun};
    
    try
        feval(functionName,z,1.0,2,1);
        fittingOptions.functionName = functionName;
        return
    catch 
        error_message = 'No working functions found. Check your MATLAB installation.';
    end
    
end

% If nothing works, throw exception
error(error_message);

end


function [P,CRLB,LL]=gaussmlev2(data,PSFSigma,iterations,fittype,function_name,Ax,Ay,Bx,By,gamma,d)
%gaussmlev2  MLE of single molecule positions 
%
%   [P CRLB LL t]=gaussmlev2(data,PSFSigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d)
%
%   This code performs a maximum likelihood estimate of particle position,
%   emission rate (Photons/frame), and background rate (Photons/pixels/frame).
%   The found values are used to calculate the Cramer-Rao Lower Bound (CRLB)
%   for each parameter and the CRLBs are returned along with the estimated parameters.
%   The input is one or more identically sized images, each of which is
%   assumed to contain a single fluorophore.  Since the maximum likelihood
%   estimate assumes Poisson statistics, images must be converted from arbitrary
%   ADC units to photon counts.
%
%       INPUTS:
%   data:       PxPxN stack of images where N is number of PxP images
%   PSFSigma:   Microscope PSF sigma (sigma_0 for z fit, starting sigma for sigma fit)
%   iterations: Number of iterations (default=10)
%   fittype:    (default=1)
%    1: XY Position, Photons, Background
%    2: XY Position, Photons, Background, PSF Sigma
%    3: XYZ Position, Photons, Background 
%    4: XY Position, Photons, Background, PSF Sigma_X, PSF Sigma_Y
%   
%       OUTPUTS:
%   P:      Found parameter values by fittype:
%    1: [X Y Photons BG] 
%    2: [X Y Photons BG Sigma] 
%    3: [X Y Photons BG Z]  
%    4: [X Y Photons BG SigmaX SigmaY] 
%   CRLB:   Cramer-Rao lower bound calculated using P
%   LL:     Log-Likelihood calculated using Stirling's approximation      
%   t:      Execution time of fitting code. 
%
%   This codes attempts first the GPU implementation (Nvidia CUDA),
%   then the c implementation, and finally the matlab implementation
%
%   REFERENCE:
%   "Fast, single-molecule localization that
%   achieves theoretical optimal accuracy." Carlas S. Smith, Nikolai Joseph,
%   Bernd Rieger and Keith A. Lidke

% functions = {'gaussmlev2_cuda50','gaussmlev2_cuda42','gaussmlev2_cuda40',...
%     'gaussmlev2_c_thread','gaussmlev2_c','gaussmlev2_matlab'};
% 
% if nargin<3
%     iterations=10;
% end
% if nargin<4
%     fittype=1;
% end
% if nargin<2
%     error('Minimal usage: gaussmlev2(data,PSFSigma)');
% end
% 
% %convert dipimage data to single
% if isa(data,'dip_image')
%     data=permute(single(data),[2 1 3]);
% end
% 
% %convert other matlab datatyes to single
% if ~ (isa(data,'double'))
%     data=single(data);
% end
% 
% if ~(isa(data,'single'))
%     error('Input data must be type dipimage, single, or double')
% end

% for iFun=1:numel(functions)
%     function_name = functions{iFun};
%     try
        switch fittype
            case 1
                [P, CRLB, LL]=feval(function_name,data,PSFSigma,iterations,fittype);
            case 2
                [P, CRLB, LL]=feval(function_name,data,PSFSigma,iterations,fittype);
            case 3
                %not reachable here
                [P, CRLB, LL]=feval(function_name,data,PSFSigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d);
            case 4
                [P, CRLB, LL]=feval(function_name,data,PSFSigma,iterations,fittype);
        end
%         return
%     catch ME
%         % do nothing
%     end
% end
% 
% throw(ME);

end
