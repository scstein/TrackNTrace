function [plugin_name, plugin_type, plugin_info] = plugin_TNTfitter(h_panel, inputOptions)
%    -------------- TNT core code, not to change by user --------------    
    if nargin < 2
        inputOptions = [];
    end
    
    % This stores the setup of all parameters
    param_specification = cell(0,4);
    

%    -------------- User definition of plugin --------------

    % Name of the component these options are for
    plugin_name = 'TNT Fitter';

    % Type of plugin.
    % 1: Candidate detection
    % 2: Spot fitting
    % 3: Tracking
    plugin_type = 2;
    
    % Description of plugin, supports sprintf format specifier like '\n' for a newline
    plugin_info = ['Fast Gaussian PSF fitting implemented in C++.\n\n', ...
                   'The fitting code utilizes the ceres-solver library for optimization currently developed by Google (2015).'];
    
    % The function this plugin implements
    plugin_function =  @fitPositions_psfFitCeres;
    
    % Add parameters
    % read comments of function subfun/add_plugin_param for HOWTO
    add_param('PSFsigma',...
          'float',...
          {1.2, 0,inf},...
          'Standard deviation of the PSF in [pixels]. sigma = FWHM/(2*sqrt(2*log(2)))');
    add_param('fitPSFsigma',...
              'bool',...
              false,...
              'Controls if the PSF standard deviation is optimized by the fitting routine (true) or kept fixed (false).');
    add_param('usePixelIntegratedFit',...
              'bool',...
              true,...
              'Use a pixel integrated Gaussian PSF for fitting (true, recommended due to higher accuracy) or not.');
    add_param('useMLE',...
              'bool',...
              false,...
              'Use Maximum Likelihood Estimation in addition to Least-squares optimization (true) or not (false).');
    
          
%   -------------- TNT core code, not to change by user --------------
%               
    % Calling the plugin function without arguments just returns its name, type and info
    if (nargin == 0); return; end
    
    % Create the panel for this plugin
    createOptionsPanel(h_panel, plugin_name, param_specification, inputOptions);

    % Save handle of the plugins function
    options = getappdata(h_panel,'options');
    options.functionHandle = plugin_function;
    setappdata(h_panel,'options',options);
    
    function add_param(par_name, par_type, par_settings, par_tooltip)
        param_specification = add_plugin_param(param_specification, par_name, par_type, par_settings, par_tooltip);
    end
end


%   -------------- User functions --------------

function [fitData] = fitPositions_psfFitCeres(img,candidatePos,options)
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
%     [x,y,A,B,sigma,flag]. Refer to locateParticles.m or to TrackNTrace
%     manual for more information.  

varsToFit = [1,1,1,1,options.fitPSFsigma];
halfw = round(4*options.PSFsigma);

[params] = psfFit_Image( img, candidatePos.',varsToFit,options.usePixelIntegratedFit,options.useMLE,halfw,options.PSFsigma);
fitData = {params(:,params(end,:)==1).'};

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
%   params     -  Fitted parameters 6xN. Columns are in order
%                 [xpos; ypos; A; BG; sigma; exitflag].
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
[ params ] = mx_psfFit_Image( double(img), varargin{:} );

end

