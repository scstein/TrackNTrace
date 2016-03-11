function [plugin] = plugin_TNTfitter()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'TNT Fitter';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 2;

% The functions this plugin implements
mainFunc =  @fitPositions_psfFitCeres;

% Description of output parameters
outParamDescription = {' '}; % set depending on plugin options in init function

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Additional functions to call before and after main function
plugin.initFunc = @fitPositions_psfFitCeres_consolidateOptions;
plugin.postFunc = @fitPositions_psfFitCeres_calculateZ;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Fast Gaussian PSF fitting implemented in C++.\n\n', ...
    'The fitting code utilizes the ceres-solver library for optimization currently developed by Google (2015).'];

% Deactivate TNT's parallel processing, as fitter is parallelized in C++
plugin.useParallelProcessing = false;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('PSFsigma',...
    'float',...
    {1.3, 0,inf},...
    'Standard deviation of the PSF in [pixels]. sigma = FWHM/(2*sqrt(2*log(2)))');
plugin.add_param('fitType',...
    'list',...
    {'[x,y,A,BG]', '[x,y,A,BG,s]', '[x,y,A,BG,sx,sy]', '[x,y,A,BG,sx,sy,angle]'},...
    'Fit positions (xy), amplitude & background (A,BG), sigma (s, or sx & sy for elliptic Gaussian), and angle (for a rotated elliptic Gaussian).');
plugin.add_param('usePixelIntegratedFit',...
    'bool',...
    true,...
    'Use a pixel integrated Gaussian PSF for fitting (true, recommended due to higher accuracy) or not.');
plugin.add_param('useMLE',...
    'bool',...
    false,...
    'Use Maximum Likelihood Estimation in addition to Least-squares optimization (true) or not (false).');
plugin.add_param('astigmaticCalibrationFile',...
    'filechooser',...
    {'','mat'},...
    'For astigmatic imaging, please select a calibration file created with z-calibration plugin. Leave empty otherwise!');

end


%   -------------- User functions --------------

function [fittingData] = fitPositions_psfFitCeres(img,candidatePos,options,currentFrame)
% Wrapper function for psfFit_Image function (see below). Refer to
% tooltips above and to psfFit_Image help to obtain information on input
% and output variables.
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
%     fittingData: 1x1 cell of 2D double array of fitted parameters
%     [x,y,z,A,B,[other parameters]]. Other parameters can be q1, q2, q3
%     (refer to locateParticles.m or to TrackNTrace manual for more
%     information). q_i will be calculated back to sigma_x,sigma_y,
%     rotation angle and possibly z in post-processing function (see
%     below).


[params] = psfFit_Image( img, candidatePos.',options.varsToFit,options.usePixelIntegratedFit,options.useMLE,options.halfw,options.PSFsigma);
params = [params(1:2,:);zeros(1,size(params,2));params(3:end,:)]; %adding z = 0
fittingData = params(1:end-1,params(end,:)==1).';

end


function [ params ] = psfFit_Image( img, varargin )
% Short usage: params = psfFit_Image( img, param_init );
% Full usage: [ params, exitflag ] = psfFit_Image( img, param_init, param_optimizeMask, useIntegratedGauss, useMLErefine, hWinSize, global_init )
%  Fit multiple spot candidates in a given image with a gaussian point spread function model.
%
% Coordinate convention: Integer coordinates are centered in pixels. I.e.
% the position xpos=3, ypos=5 is exactly the center of pixel img(5,3). Thus
% pixel center coordinates range from 1 to size(img,2) for x and 1 to
% size(img,1) for y coordinates.
%
% Use empty matrix [] for parameters you don't want to specify.
%
% Gaussian fitting covers three basic cases:
%   - fitting isotropic gaussian  (sigma_y and angle initial values not specified and they should not be optimized)
%   - fitting anisotropic gaussian  (angle initial value not specified and should not be optimized)
%   - fitting anisotropic rotated gaussian
%
% The fitting case is selected based on the set of specified initial
% parameters together with the set of parameters that should be optimized.
% The angle input/output should be in degree.
%
% Input:
%   img        - Image to fit to. (internally converted to double)
%   param_init - Initial values PxN to fit N spot candidates in the given image with
%                initial conditions for P parameters. At least position [xpos; ypos]
%                must be specified (P>=2). You can specify up to [xpos;ypos;A;BG;sigma_x,sigma_y,angle].
%                If negative or no values are given, the fitter estimates a value for that parameter if neccessary, 
%                with the exception of the angle, which is estimated if angle==0 and if it should be optimized.
%                If parameters are not optimized (see next entry) their values are kept constant during optimization.
%   param_optimizeMask - Must be true(1)/false(0) for every parameter [xpos,ypos,A,BG,sigma_x,sigma_y,angle].
%                Parameters with value 'false' are not fitted.  | default: [1,1,1,1,1,0,0] -> 'optimize x,y,A,BG,sigma' (isoptric gaussian)
%   useIntegratedGauss - Wether to use pixel integrated gaussian or not 
%                  not supported for non-isotropic arbitrarily roated gaussians | default: false
%   useMLErefine - Use Poissonian noise based maximum likelihood estimation after
%                  least squares fit. Make sure to input image intensities in photons
%                  for this to make sense. | default: false
%   hWinSize   - Each candidates fit takes intensites in a window of (2*hWinsize+1)x(2*hWinsize+1) into account | default: 5
%   global_init - For convenience (up to) [sigma_x,sigma_y,angle] can also be given as an extra parameter.
%               This simply sets all candidates initial values, eventually overwriting their param_init values.
%
% Output
%   params     -  Fitted parameters 8xN. Columns are in order
%                 [xpos,ypos,A,BG,sigma_x,sigma_y,angle; exitflag].
%
%           The last row 'exitflag' returns the state of optimizer.
%           Positive = 'good'. Negative = 'bad'.
%             1 - CONVERGENCE
%            -1 - NO_CONVERGENCE
%            -2 - FAILURE
%
% Authors: Simon Christoph Stein and Jan Thiart
% Date:   March 2016
% E-Mail: scstein@phys.uni-goettingen.de

% Copyright (c) 2016, Simon Christoph Stein
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The views and conclusions contained in the software and documentation are those
% of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of the FreeBSD Project.
%
%
% -- Licensing --
%
% License ceres-solver:
%
% Copyright 2015 Google Inc. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%
%     1. Redistributions of source code must retain the above copyright notice, this list
% of conditions and the following disclaimer.
%     2. Redistributions in binary form must reproduce the above copyright notice, this list
% of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%     3. Neither the name of Google Inc., nor the names of its contributors may be used to
% endorse or promote products derived from this software without specific prior written permission.
%
% This software is provided by the copyright holders and contributors "AS IS" and any express or
% implied warranties, including, but not limited to, the implied warranties of merchantability and
% fitness for a particular purpose are disclaimed. In no event shall Google Inc. be liable for any
% direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited
% to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption)
% however caused and on any theory of liability, whether in contract, strict liability, or tort (including
% negligence or otherwise) arising in any way out of the use of this software, even if advised of the
% possibility of such damage.

% Make sure logicals are passed as correct datatype
if numel(varargin) >= 2;  varargin{2} = logical(varargin{2});  end
if numel(varargin) >= 3;  varargin{3} = logical(varargin{3});  end
if numel(varargin) >= 4;  varargin{4} = logical(varargin{4});  end

% Set default optimization x,y,A,BG,sigma (isotropic gaussian)
if( numel(varargin)<2 || isempty(varargin{2}) )
    varargin{2} = logical([1,1,1,1,1,0,0]);
end

% Convert img to double if neccessary
[ params ] = mx_psfFit_Image( double(img), varargin{:} );

end


function [fittingOptions] = fitPositions_psfFitCeres_consolidateOptions(fittingOptions)
%initializer function
% plugin.add_param('fitType',...
%     'list',...
%     {'[x,y,A,BG]', '[x,y,A,BG,s]','[x,y,A,BG,sx,sy], [x,y,A,BG,sx,sy,angle]'},...

if ~isempty(fittingOptions.astigmaticCalibrationFile)
    cal_struct = load(fittingOptions.astigmaticCalibrationFile,'calibrationData');
    if isfield(cal_struct,'calibrationData')
        fittingOptions.calibrationData = cal_struct.calibrationData;
        fittingOptions.fitType = '[x,y,A,BG,sx,sy]';
    else
        error('%s is not a valid calibration file. Aborting.',fittingOptions.astigmaticCalibrationFile);
    end
end

switch fittingOptions.fitType
    case '[x,y,A,BG]'
        varsToFit = [ones(4,1);zeros(3,1)];
    case '[x,y,A,BG,s]'
        varsToFit = [ones(5,1);zeros(2,1)];
    case '[x,y,A,BG,sx,sy]'
        varsToFit = [ones(6,1);0];
    case '[x,y,A,BG,sx,sy,angle]'
        varsToFit = ones(7,1);
    otherwise
        warning off backtrace
        warning('Unrecognized fit type. Switching to [x,y,A,BG].');
        warning on backtrace
        varsToFit = [ones(4,1);zeros(3,1)];
end

if varsToFit(7) == 1 %fit angle?
    if fittingOptions.usePixelIntegratedFit
        fittingOptions.usePixelIntegratedFit = false;
        warning off backtrace
        warning(sprintf('Fitting a rotation angle is not possible while using pixel-integrated Gaussian model. \nSwitching to sampled Gaussian.'));
        warning on backtrace
    end
end

%always fit x,y,A,BG, determine if one has to fit sigma_x,
%sigma_y,theta_rot
fittingOptions.halfw = round(3*fittingOptions.PSFsigma);
fittingOptions.varsToFit = varsToFit;

% updating parameter description
switch sum(varsToFit(5:end))
    case 0 % Sigma not fitted
        fittingOptions.outParamDescription = {'x';'y';'z';'Amp (Peak)'; 'Background'; 'sigma'};
    case 1
        fittingOptions.outParamDescription = {'x';'y';'z';'Amp (Peak)'; 'Background'; 'sigma'};
    case 2
        fittingOptions.outParamDescription = {'x';'y';'z';'Amp (Peak)'; 'Background'; 'sigma_x'; 'sigma_y'};
    case 3
        fittingOptions.outParamDescription = {'x';'y';'z';'Amp (Peak)'; 'Background'; 'sigma_x'; 'sigma_y'; 'angle [°]'};
    otherwise
        warning('TNT fitter init func: Unknown case');
end

end %consolidateOptions


function [fittingData,fittingOptions] = fitPositions_psfFitCeres_calculateZ(fittingData,fittingOptions)
% If astigmatic imaging was performed, the axial position is
% extrapolated from the calibration file. The axial position is given in
% pixels (z in nm = z*fittingOptions.calibrationData.zPixel).
%     
% The rotation angle is associated with the axis of the largest eigenvalue
% given by [-pi/4:pi/4]. Therefore, an elliptic Gaussian with its major
% axis closer to the y-axis will have a larger sigma_y and a smaller
% sigma_x value and the rotation angle shows how much a hypothetical
% ellipse whose major axis is perfectly aligned with the y-axis is rotated.
% Thus, rotation and sigma values are uniquely defined up to a symmetric
% 180° rotation.
% A positive angle corresponds to a counter-clockwise rotation
% (mathematically positive).

calibrationFileExists = isfield(fittingOptions,'calibrationData');
emptyFrames = cellfun('isempty',fittingData);

if ~calibrationFileExists && sum(fittingOptions.varsToFit(end-1:end))==0
    fittingData(~emptyFrames) = cellfun(@(var) var(:,1:6),fittingData(~emptyFrames),'UniformOutput',false);
    return
end


% Do we have to calculate z or correct the angle?
for iFrame = 1:numel(fittingData)
    if emptyFrames(iFrame)
        continue
    end
    
    fittingData_frame = fittingData{iFrame};
    if fittingOptions.varsToFit(end)
        sigma_x = fittingData_frame(:,6);
        sigma_y = fittingData_frame(:,7);
        angle = fittingData_frame(:,8);
        
        % if the angle is 0 or very close to 0, calculating sigma will fail
        idxFaultyValues = (angle == 0 | sum(isnan([sigma_x,sigma_y]),2)>0);
        
        fittingData_frame(idxFaultyValues,:) = [];
    end
    
    if calibrationFileExists
        fittingData_frame(:,3) = interp1(fittingOptions.calibrationData.aspectRatioSigmaSmooth(:,2),...
            fittingOptions.calibrationData.aspectRatioSigmaSmooth(:,1),...
            fittingData_frame(:,6)./fittingData_frame(:,7),'linear','extrap')-repmat(fittingOptions.calibrationData.zMidpoint,size(fittingData_frame,1),1);
        %we have sigma_x/sigma_y(z) curve -> interpolate from inverse function
        %alternative: evaluate polynom and find minimum of curve difference
        fittingData_frame = fittingData_frame(:,1:7); %angle not fitted when determining z
    end
    
    fittingData{iFrame} = fittingData_frame;
end


end %calculateZ
