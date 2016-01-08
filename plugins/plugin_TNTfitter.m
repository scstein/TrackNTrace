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

% Create the plugin
plugin = TNTplugin(name, type, mainFunc);

% Additional functions to call before and after main function
plugin.initFunc = @consolidateOptions;
plugin.postFunc = @calculateZRot;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Fast Gaussian PSF fitting implemented in C++.\n\n', ...
    'The fitting code utilizes the ceres-solver library for optimization currently developed by Google (2015).'];

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
plugin.add_param('PSFsigma',...
    'float',...
    {1.2, 0,inf},...
    'Standard deviation of the PSF in [pixels]. sigma = FWHM/(2*sqrt(2*log(2)))');
plugin.add_param('fitPSFsigma',...
    'bool',...
    false,...
    'Controls if the PSF standard deviation is optimized by the fitting routine (true) or kept fixed (false).');
plugin.add_param('usePixelIntegratedFit',...
    'bool',...
    true,...
    'Use a pixel integrated Gaussian PSF for fitting (true, recommended due to higher accuracy) or not.');
plugin.add_param('useMLE',...
    'bool',...
    false,...
    'Use Maximum Likelihood Estimation in addition to Least-squares optimization (true) or not (false).');
plugin.add_param('fitAngle',...
    'bool',...
    false,...
    'Fit rotation angle of anisotropic Gaussian function.');
plugin.add_param('astigmaticCalibrationFile',...
    'filechooser',...
    {'','mat'},...
    'For astigmtic imaging, please select a calibration file created with z-calibration plugin. Leave empty otherwise!');

end


%   -------------- User functions --------------

function [fitData] = fitPositions_psfFitCeres(img,candidatePos,options,currentFrame)
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
%     fitData: 1x1 cell of 2D double array of fitted parameters
%     [x,y,z,A,B,[other parameters],flag]. Other parameters can be q1, q2,
%     q3 (refer to locateParticles.m or to TrackNTrace manual for more
%     information). q_i will be calculated back to sigma_x,sigma_y and
%     rotation angle later.


% varsToFit = [ones(6,1);0]; %fit everything except angle
% halfw = round(4*options.PSFsigma);

[params] = psfFit_Image( img, candidatePos.',options.varsToFit,options.usePixelIntegratedFit,options.useMLE,options.halfw,options.PSFsigma);
params = [params(1:2,:);zeros(1,size(params,2));params(3:end,:)]; %adding z = 0
fitData = params(:,params(end,:)==1).';

end

function [ params ] = psfFit_Image( img, varargin)
% Short usage: params = psfFit_Image( img, param_init );
% Full usage: [ params, exitflag ] = psfFit_Image( img, param_init, param_optimizeMask, useIntegratedGauss, hWinSize, sigma_init )
%  Fit multiple spot candidates in a given image with a gaussian point spread function model.
%
% Coordinate convention: Integer coordinates are centered in pixels. I.e.
% the position xpos=3, ypos=5 is exactly the center of pixel img(5,3). Thus
% pixel center coordinates range from 1 to size(img,2) for x and 1 to
% size(img,1) for y coordinates.
%
% Use empty matrix [] for parameters you don't want to specify.
%
% Input:
%   img        - Image to fit to. (internally converted to double)
%   param_init - Initial values PxN to fit N spot candidates in the given image with
%                initial conditions for P parameters. At least position [xpos; ypos]
%                must be specified (P>=2). You can specify up to [xpos;ypos;A;BG;sigma].
%                If negative values are given, the fitter estimates a value for that parameter.
%   param_optimizeMask - Must be true(1)/false(0) for every parameter [xpos,ypos,A,BG,sigma].
%                Parameters with value 'false' are not fitted. | default: ones(5,1) -> 'optimize all'
%   useIntegratedGauss - Wether to use pixel integrated gaussian or not | default: 'false'
%   useMLErefine - Use Poissonian noise based maximum likelihood estimation after
%                  least squares fit. Make sure to input image intensities in photons
%                  for this to make sense. | default: false
%   hWinSize   - window around each candidate is (2*hWinsize+1)x(2*hWinsize+1) | default: 5
%   sigma_init - For convenience sigma can also be given as an extra parameter.
%               This simply sets all candidates initial sigma to sigma_init.
%               This overwrites the value given in param_init.
%
% Output
%   params     -  Fitted parameters 8xN. Columns are in order
%                 [xpos; ypos; A; BG; q1; q2; q3; exitflag].
%
%           The last row 'exitflag' returns the state of optimizer.
%           Positive = 'good'. Negative = 'bad'.
%             1 - CONVERGENCE
%            -1 - NO_CONVERGENCE
%            -2 - FAILURE
%
% Author: Simon Christoph Stein
% Date:   June 2015
% E-Mail: scstein@phys.uni-goettingen.de

% Copyright (c) 2015, Simon Christoph Stein
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
% This software is provided by the copyright holders and contributors “AS IS” and any express or
% implied warranties, including, but not limited to, the implied warranties of merchantability and
% fitness for a particular purpose are disclaimed. In no event shall Google Inc. be liable for any
% direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited
% to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption)
% however caused and on any theory of liability, whether in contract, strict liability, or tort (including
% negligence or otherwise) arising in any way out of the use of this software, even if advised of the
% possibility of such damage.

% Make sure logicals are passed as correct datatype
if numel(varargin) >= 3;  varargin{3} = logical(varargin{3});  end
if numel(varargin) >= 4;  varargin{4} = logical(varargin{4});  end

% Convert img to double if neccessary
[ params ] = mx_psfFit_ImageNEW( double(img), varargin{:} );

end


function [fittingOptions] = consolidateOptions(fittingOptions)
%initializer function

fitBothSigma = false;

if ~isempty(fittingOptions.astigmaticCalibrationFile)
    cal_struct = load(fittingOptions.astigmaticCalibrationFile,'calibrationData');
    if isfield(cal_struct,'calibrationData')
        fittingOptions.fitPsfsigma = true;
        fitBothSigma = true;
        fittingOptions.calibrationData = cal_struct.calibrationData;
    else
        error('%s is not a valid calibration file. Aborting.',fittingOptions.astigmaticCalibrationFile);
    end
end

if fittingOptions.fitAngle
    fittingOptions.fitPSFsigma = true;
    fitBothSigma = true;
    if fittingOptions.usePixelIntegratedFit
        fittingOptions.usePixelIntegratedFit = false;
        warning(sprintf('Fitting a rotation angle is not possible while using pixel-integrated Gaussian model. \nSwitching to sampled Gaussian.'));
    end
end

%always fit x,y,A,BG, determine if one has to fit sigma_x,
%sigma_y,theta_rot
fittingOptions.halfw = round(4*fittingOptions.PSFsigma);
fittingOptions.varsToFit = [ones(4,1);fittingOptions.fitPSFsigma;fitBothSigma;fittingOptions.fitAngle];

end %consolidateOptions


function [fitData,fittingOptions] = calculateZRot(fitData,fittingOptions)
%post-processing function

calibrationFileExists = isfield(fittingOptions,'calibrationData');
emptyFrames = cellfun('isempty',fitData);

% without astigmatism or rotation, there is only one sigma value sigma =
% 1/sqrt(2*q1) which has to be corrected
if ~calibrationFileExists && ~fittingOptions.fitAngle
    % back calculate sigma, delete q2, q3
    fitData(~emptyFrames) = cellfun(@(var) [var(:,1:5),1./sqrt(2*var(:,6)),var(:,end)],fitData(~emptyFrames),'UniformOutput',false);
    return
end

% otherwise, the calculation is more complicated. Go through frame by
% frame...
for iFrame = 1:numel(fitData)
    if emptyFrames(iFrame)
        continue
    end
    
    fitData_frame = fitData{iFrame};
    if fittingOptions.fitAngle
        q1 = fitData_frame(:,6); q2 = fitData_frame(:,7); q3 = fitData_frame(:,8);
        
        angle = 0.5*atan(2*q3./(q2-q1)); %angle towards largest eigenvector axis in radian
        sigma_x = 1./sqrt(q1+q2-2*q3./sin(2*angle));
        sigma_y = 1./sqrt(q1+q2+2*q3./sin(2*angle));
        
        % if the angle is 0 or very close to 0, calculating sigma will fail
        idxFaultyValues = (angle == 0 | sum(isnan([sigma_x,sigma_y]),2)>0);
        sigma_x(idxFaultyValues) = 1./sqrt(2*q1(idxFaultyValues));
        sigma_y(idxFaultyValues) = 1./sqrt(2*q2(idxFaultyValues));
        angle(idxFaultyValues) = 0;
        
        fitData_frame(:,6:8) = [sigma_x,sigma_y,angle];
    else
        fitData_frame = [fitData_frame(:,1:5),1./sqrt(2*fitData_frame(:,6:7)),fitData_frame(:,end)]; %delete q3, calc. sigma_x, sigma_y
    end
    
    if calibrationFileExists
        fitData_frame(:,3) = interp1(fittingOptions.calibrationData.aspectRatioSigmaSmooth(:,2),...
            fittingOptions.calibrationData.aspectRatioSigmaSmooth(:,1),...
            fitData_frame(:,6)./fitData_frame(:,7),'linear','extrap')-repmat(fittingOptions.calibrationData.zMidpoint,size(fitData_frame,1),1);
        %we have sigma_x/sigma_y(z) curve -> interpolate from inverse function
        %alternative: evaluate polynom and find minimum of curve difference
    end
    
    fitData{iFrame} = fitData_frame;
end

end %calculateZRot
