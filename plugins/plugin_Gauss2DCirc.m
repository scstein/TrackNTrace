function [plugin] = plugin_Gauss2DCirc()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Gauss-2D-Circ';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 2;

% The functions this plugin implements
mainFunc =  @fitParticles_gauss2dcirc;

% Description of output parameters
outParamDescription = {'x';'y';'z';'Amp (Peak)'; 'Background'; 'width'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Add initial function
plugin.initFunc = @fitParticles_gauss2dcirc_init;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Refine candidate positions by assuming a Gaussian PSF and finding the center by matrix inversion/linear regression using gauss2dcirc. \n\nThe algorithm was published in ''Anthony, S.M. & Granick, S. Image analysis with rapid and accurate two-dimensional Gaussian fitting. Langmuir 25, 8152–8160 (2009)''.';

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('PSFSigma',...
    'float',...
    {1.3, 0, inf},...
    'Standard deviation of the PSF in pixels. \nsigma = FWHM/(2*sqrt(2*log(2))) ~ 0.21*lambda/NA where lambda is the emission wavelength in pixels and NA is the numerical aperture of the objective.');
plugin.add_param('backgroundNoiseLevel',...
    'float',...
    {10, 0, inf},...
    'Standard deviation of background noise in [ADU] (camera counts). \nThis determines how the pixel intensities are weighted and should match the experiment.');
end


function [fittingData] = fitParticles_gauss2dcirc(img,candidatePos,options,currentFrame)
% Wrapper function for gauss2dcirc function (see below). Refer to tooltips
% above and to gauss2dcirc help to obtain information on input and output
% variables. gauss2dcirc.m was released as part of the following
% publication under the GNU public licence: 
% Anthony, S.M. & Granick, S. Image analysis with rapid and accurate
% two-dimensional Gaussian fitting. Langmuir 25, 8152–8160 (2009).
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
%     [x,y,z,A,B,particle width]. 

fittingData = zeros(size(candidatePos,1),6);
[x_mesh,y_mesh] = meshgrid(1:size(img,1),1:size(img,2));


for iCand = 1:size(candidatePos,1)
    candPos = round(candidatePos(iCand,:));
    idx_x = max(1,candPos(1)-options.windowHalfSize):min(size(img,1),candPos(1)+options.windowHalfSize);
    idx_y = max(1,candPos(2)-options.windowHalfSize):min(size(img,1),candPos(2)+options.windowHalfSize);
    
    [xc,yc,Amp,width] = gauss2dcirc(img(idx_y,idx_x),x_mesh(idx_y,idx_x),y_mesh(idx_y,idx_x),options.backgroundNoiseLevel);
    fittingData(iCand,:) = [xc,yc,0,Amp,0,width];
end
fittingData = fittingData(fittingData(:,1)>0,:);

end


function [fittingOptions] = fitParticles_gauss2dcirc_init(fittingOptions)

fittingOptions.windowHalfSize = ceil(3*fittingOptions.PSFSigma);

end


function [xc,yc,Amp,width]=gauss2dcirc(z,x,y,noiselevel)
% function downloaded from: http://groups.mrl.uiuc.edu/granick/software.html
% 
%[xc,yc,Amp,width]=gauss2dcirc(arr,x,y,noiselevel)
%
%GAUSS2DCIRC.m attempts to find the best 2D circular Gaussian fit for the
%selected region of the image. 
%
%%Copyright 2008 Stephen M. Anthony, U. Illinois Urbana-Champaign
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program. If not, see <http://www.gnu.org/licenses/>
%
%INPUTS:    Z:          The subregion of the image to be fit, with each
%                       element giving the value of one pixel. 
%           X:          This matrix or vector should be the same size as
%                       arr, with each element giving the x-coordinate of
%                       the corresponding pixel. 
%           Y:          This matrix or vector should be the same size as
%                       arr, with each element giving the y-coordinate of
%                       the corresponding pixel. 
%           NOISELEVEL: The standard deviation of the background noise. 
%
%OUTPUTS:   XC,YC:      The center of the Gaussian in x and y
%           Width:      The width of the Gaussian. 
%           Amp:        The amplitude of the gaussian. 
%
%Written by Stephen Anthony 1/2007 U. Illinois Urbana-Champaign
%Last Modified by Stephen Anthony on 12/16/2008

%Convert to column form, in case it is not already
x=x(:); y=y(:);
Z=z(:)+1e-15;
%The miniscule offset on Z has no significant effect, but is a shortcut to
%ensure that no element of Z equals 0. The probability that an element of
%arr equals 0 is not insignificant (as we may be working with integers, the
%probability than an element equals exactly -1e-15 is negligible. 

%Compute the weighting. This is approximate, as we are working in the log
%scale, where the noise will be assymmetric. Bounds were placed on this to
%avoid logs of negative numbers, but for any case the bound triggers, the
%weight will be quite low anyway. 
noise=log(max(Z+noiselevel,1e-10))-log(max(Z-noiselevel,1e-20));
wght=(1./noise);
wght(Z<=0)=0;

n=[x y log(Z) ones(size(x))].*(wght*ones(1,4));
d=-(x.^2+y.^2).*wght;
a=n\d;
%In the least squares sense, a was selected such that 
%   n*a = d
%or the best possible match thereof. 

%Extract the desired values
xc = -.5*a(1);
yc = -.5*a(2);
width=sqrt(a(3)/2);
Amp=exp((a(4)-xc^2-yc^2)/(-2*width^2));
end