function [plugin] = plugin_GPUGauss_MLE()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'GPU-Gauss MLE';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
type = 2;

% The functions this plugin implements
mainFunc = @refinePositions_gpugaussmle;

% Description of output parameters
outParamDescription = {'x';'y';'z';'Intens. (integ.)'; 'Background'; 'sigma_x'; 'sigma_y'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Add initial function
plugin.initFunc = @refinePositions_gpugaussinit;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Refine candidate positions by assuming a Gaussian PSF and finding the center and other parameters via GPU-MLE fitting. Absolutely requires photon conversion. \n\nFunction published in ''Smith et al, NatMet 7, 373-375(2010), doi:10.1038/nmeth.1449''.';

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('PSFSigma',...
    'float',...
    {1.3, 0, inf},...
    'Standard deviation of the PSF in pixels. \nsigma = FWHM/(2*sqrt(2*log(2))) ~ 0.21*lambda/NA where lambda is the emission wavelength in pixels and NA is the numerical aperture of the objective.');
plugin.add_param('fitType',...
    'list',...
    {'[x,y,N,BG]', '[x,y,N,BG,s]','[x,y,N,BG,sx,sy]'},...
    'Fit positions (xy), amplitude & background (A,BG), and sigma (s, or sx & sy for elliptic Gaussian).');
plugin.add_param('Iterations',...
    'int',...
    {10, 3, inf},...
    'Number of fit iterations.');

end


function [refinementData] = refinePositions_gpugaussmle(img,candidatePos,options,currentFrame)
% Wrapper function for gaussmlev2 function (see below). Refer to tooltips
% above and to gaussmlev2 help to obtain information on input and output
% variables. gaussmlev2.m was released as part of the following
% publication under the GNU public licence: 
% Smith et al, Fast, single-molecule localization that achieves
% theoretically minimum uncertainty, Nature Methods 7, 373-375 (2010),
% doi:10.1038/nmeth.1449
% All files can be downloaded at http://omictools.com/gaussmlev2-tool 
% (put in external folder).
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
%     refinementData: 1x1 cell of 2D double array of fitted parameters
%     [x,y,z,A,B,[other parameters]]. Other parameters can be PSF standard
%     deviation, both in x and y.

refinementData = zeros(size(candidatePos,1),5+options.nrParam);

for iCand = 1:size(candidatePos,1)
    pos = candidatePos(iCand,1:2);
    idx_x = max(1,pos(1)-options.halfWindowSize):min(size(img,2),pos(1)+options.halfWindowSize);
    idx_y = max(1,pos(2)-options.halfWindowSize):min(size(img,1),pos(2)+options.halfWindowSize);
    if numel(idx_x)~=numel(idx_y)
        continue; %gaussmlev2 cannot deal with asymmetric fit windows
    end
    [param,~,~]=gaussmlev2(single(img(idx_y,idx_x)),options.PSFSigma,options.Iterations,options.fitType,options.functionName); %[x,y,z,A,BG,[rest params]]
    refinementData(iCand,:) = [param(1:2)+[idx_x(1),idx_y(1)],0,param(3:end)]; %correct position relative to fit window, middle of first pixel is [0,0]
end

refinementData = refinementData(refinementData(:,1)>0,:);

end


function [refinementOptions] = refinePositions_gpugaussinit(refinementOptions)
global globalOptions

if ~globalOptions.usePhotonConversion
    warning off backtrace
    warning('GPU-Gauss MLE strictly requires photon conversion!');
    warning on backtrace
end

refinementOptions.halfWindowSize = min(10,ceil(3*refinementOptions.PSFSigma));
switch refinementOptions.fitType
    case '[x,y,N,BG]'
        nrParam = 0;
    case '[x,y,N,BG,s]'
        nrParam = 1;
    case '[x,y,N,BG,sx,sy]'
        nrParam = 2;
end

refinementOptions.nrParam = nrParam;
refinementOptions.outParamDescription = refinementOptions.outParamDescription(1:5+nrParam);

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
        refinementOptions.functionName = functionName;
        return
    catch 
        error_message = 'No working functions found. Check your MATLAB installation.';
    end
    
end

% If nothing works, throw exception
error(error_message);

end


function [P,CRLB,LL]=gaussmlev2(data,PSFSigma,iterations,fittype,function_name,Ax,Ay,Bx,By,gamma,d)
% modified version of gaussmlev2 which was changed to work with TrackNTrace
% gaussmlev2.m was released as part of the following
% publication under the GNU public licence: 
% Smith et al, Fast, single-molecule localization that achieves
% theoretically minimum uncertainty, Nature Methods 7, 373-375 (2010),
% doi:10.1038/nmeth.1449
% All files can be downloaded at http://omictools.com/gaussmlev2-tool (put in external folder).
%
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
            case '[x,y,N,BG]'
                [P, CRLB, LL]=feval(function_name,data,PSFSigma,iterations,1);
            case '[x,y,N,BG,s]'
                [P, CRLB, LL]=feval(function_name,data,PSFSigma,iterations,2);
            case 3
                %not reachable here
                [P, CRLB, LL]=feval(function_name,data,PSFSigma,iterations,3,Ax,Ay,Bx,By,gamma,d);
            case '[x,y,N,BG,sx,sy]'
                [P, CRLB, LL]=feval(function_name,data,PSFSigma,iterations,4);
        end
%         return
%     catch ME
%         % do nothing
%     end
% end
% 
% throw(ME);

end
