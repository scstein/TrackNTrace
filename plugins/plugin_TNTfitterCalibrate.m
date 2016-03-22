function [plugin] = plugin_TNTfitterCalibrate()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'TNT z-Calibration';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 2;

% The functions this plugin implements
mainFunc =  @refinePositions_psfFitCeres_calibrate;

% Description of output parameters
outParamDescription = {'x';'y';'Amp (Peak)'; 'Background'; 'sigma_x'; 'sigma_y'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Additional functions to call before and after main function
plugin.initFunc = @refinePositions_psfFitCeres_consolidateOptions;
plugin.postFunc = @refinePositions_psfFitCeres_createCalibrationTable;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Calibration version of TNT fitter for obtaining astigmatic z-calibration curve.\n\n', ...
               'This function creates a calibration Table [z,sigma_x/sigma_y] from an astigmatic calibration movie where a number of diffraction-limited fluorescent emitters were recorded at different axial positions with a cylindrical lens in the imaging path and using a scanning stage. \n', ...
               'It is assumed that drastic outliers (particles in focus before stage moves to defocused position) were removed by editing the movie beforehand and that the z-stage moves smoothly at linear intervals from a position below (or above) the focus to another position above (or below) the focus. Parts of the curve where the aspect ratio gradually saturates (derivative close to 0) are removed from the smoothed interpolation curve. \n \n',...
               'The functions fits a polynom of 4th degree to the aspect ratio (sigma_x/sigma_y)(z) and also finds the midpoint z-pixel where the aspect ratio is equal to 1. A smooth interpolation curve with a resolution of 0.1 z-pixels is calculated and everything is saved to a calibration file in the same folder as the movie. \n \n',... 
               'The fitting code utilizes the ceres-solver library for optimization currently developed by Google (2015).'];
% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('PSFsigma',...
      'float',...
      {1.3, 0,inf},...
      'Standard deviation of the PSF in pixels. \nsigma = FWHM/(2*sqrt(2*log(2))) ~ 0.21*lambda/NA where lambda is the emission wavelength in pixels and NA is the numerical aperture of the objective.');
plugin.add_param('usePixelIntegratedFit',...
          'bool',...
          true,...
          'Use a pixel integrated Gaussian PSF for fitting which emulates the discrete pixel grid of a camera chip. \nThis is automatically disabled when fitting a rotated Gaussian..');
plugin.add_param('useMLE',...
          'bool',...
          false,...
          'Use Maximum Likelihood Estimation in addition to Least-squares optimization. \nMLE fitting absolute requires photon conversion!.');  
plugin.add_param('zInterval',...
          'float',...
          {15,0,inf},...
          'z-scan interval used for obtaining calibration stack in [nm].');  
end


%   -------------- User functions --------------

function [fittingData] = refinePositions_psfFitCeres_calibrate(img,candidatePos,options,currentFrame)
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
%     fittingData: 2D double array of fitted parameters
%     [x,y,A,B,sigma_x,sigma_y] used for creating calibration
%     file. Refer to locateParticles.m or to TrackNTrace manual for more
%     information.



[params] = psfFit_Image( img, candidatePos.',options.varsToFit,options.usePixelIntegratedFit,options.useMLE,options.halfw,options.PSFsigma);
fittingData = params(1:end-2,params(end,:)==1).'; %delete exitflag, don't save angle

end

function [fittingOptions] = refinePositions_psfFitCeres_consolidateOptions(fittingOptions)

fittingOptions.varsToFit = [ones(6,1);0]; %don't fit angle
fittingOptions.halfw = ceil(3*fittingOptions.PSFsigma);

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


function [fittingData,fittingOptions] = refinePositions_psfFitCeres_createCalibrationTable(fittingData,fittingOptions)
% This function creates a calibration Table [z,sigma_x/sigma_y] from an
% astigmatic calibration movie where a number of diffraction-limited
% fluorescent emitters was recorded at different axial positions with a
% scanning stage.
% 
% It's assumed that drastic outliers (particles in focus before stage moves
% to defocused position) were removed by editing the movie beforehand and
% that the z-stage moves smoothly at linear intervals from a position below
% (or above) the focus to another position above (or below) the focus.
% Parts of the curve where the aspect ratio gradually saturates (derivative
% close to 0) are removed from the smoothed interpolation curve.
% 
% The functions fits a polynom of 4th degree to the aspect ratio
% (sigma_x/sigma_y)(z) and also finds the midpoint z-pixel where the aspect
% ratio is equal to 1. A smooth interpolation curve with a resolution of
% 0.1 z-pixels is calculated and everything's saved to a calibration file
% in the same folder as the movie.


% post-processing function
global filename_movie

nrFrames = size(fittingData,1);
emptyFrames = cellfun('isempty',fittingData);

% calculate aspect ratio vs z pixel
sigma_xy = NaN(nrFrames,2);
sigma_values = vertcat(cellfun(@(var) [mean(var(:,5),1),mean(var(:,6),1)],fittingData(~emptyFrames),'UniformOutput',false));
sigma_xy(~emptyFrames,:) = vertcat(sigma_values{:});
aspectRatioSigma = sigma_xy(:,1)./sigma_xy(:,2); %simga_x/sigma_y, empty frames show as NaN

% find calibration curve and z-middle of curve
p = polyfit((1:nrFrames).',aspectRatioSigma,4);
p_root = p; p_root(end) = p_root(end)-1;
z_root = roots(p_root);
z_root = z_root(imag(z_root)==0);
if ~isempty(z_root)
    z_root = z_root(z_root>nrFrames*0.1 & z_root<nrFrames*0.9); %middle point must be in the inner 80% of the whole z-interval
    if numel(z_root)~=1
        z_root = 0; %if there's no unique solution, disable correction by middle point
    end
end

% cut off part where derivative of interpolation function gets close to
% zero or where derivative changes sign. Instead, only pick part of the
% calibration curve which is monotonous
smooth_interval = (1:0.1:nrFrames).';
numder = polyval(polyder(p),smooth_interval); 
[~,idxMax] = max(abs(numder));
numder = numder/numder(idxMax); %numder is now mostly positive definite and normalized to max = 1.0
idxGood = numder>0.15;
%now find the longest interval of appropriate values. This ensures monotonicity
idx_delta = diff([0;idxGood;0]);
idx_starts = find(idx_delta > 0);
idx_ends = find(idx_delta < 0) - 1;
idx_lengths = idx_ends-idx_starts+1;
startIdx = idx_starts(idx_lengths == max(idx_lengths));
endIdx = idx_ends(idx_lengths == max(idx_lengths));

% save calibration file next to movie and TNT file
calibrationData.aspectRatioSigma = aspectRatioSigma;
calibrationData.aspectRatioSigmaSmooth = [smooth_interval(startIdx:endIdx),polyval(p,smooth_interval(startIdx:endIdx))];
calibrationData.sigmaTable = sigma_xy;
calibrationData.zPixel = fittingOptions.zInterval;
calibrationData.polynomFunc = p;
calibrationData.zMidpoint = z_root;

[movie_path,movie_name,~] = fileparts(filename_movie);
save([movie_path,filesep,movie_name,'_calibrationData_TNT.mat'],'calibrationData');

end
