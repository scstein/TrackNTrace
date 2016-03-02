function [experimParam,trajParam,fitParam,printParam] = setDefaultOptions()
% [experimParam,trajParam,fitParam,printParam] = setDefaultOptions()
% Set all options for diffusion analysis.
% 
% OUTPUT:
%     experimParam: Struct of options for experimental parameters.
%     
%     trajParam: Struct of options for trajectory handling.
%     
%     fitParam: Struct of options for fitting parameters.
%     
%     printParam: Struct of options for printing/plotting parameters.


%Experimental parameters
experimParam.diffusionCoeffInitial = [12]; %double, [µm^2/s] expected diffusion coefficient as starting point of fit. Make array to fit additional D coeff.
experimParam.timestep = 1/934.58; %double, camera acquisition time in s
experimParam.pixel = 24; %double, width of actual camera pixel in µm
experimParam.magnification = 200/3*100/30; %double, total magnification of microscope
experimParam.temperature = 20; %double, temperature in °C
experimParam.viscosityTable = 'viscosity_water_dyn.mat'; %string, path to viscosity table [temperature,viscosity] used to renormalize D to 20°C.

%Trajectory parameters
trajParam.trajMaxLength = 60; %integer, maximum allowed length for trajectory
trajParam.frameMaxDistance = 10; %double, maximum distance a molecule travels in pixels per frame
trajParam.diffusionGuess = 0.05; %double, minimum sqrt(MSD) of trajectory
trajParam.inclusionZone = [38,32,21/2]; %[24,24,5]; %double, exclude all trajectories outside of circle with [middlePoint(1:2),radius]. Just leave as empty to not consider this

%Fitting parameters
fitParam.method = 'histogram'; %string, either 'histogram' (fits displacement histogram, then MSD curve), 'msd' (fit average msd curve), or 'msd-single' (fit msd curves for all trajectories, then create D histogram)
fitParam.maximumSkip = 10; %integer, maximum time step for trajectory, cannot be longer than the length of longest trajectory minus 1
fitParam.minimumSkip = 10; %integer, minimum time step later used to fit MSD curve. Cannot be smaller than 4
fitParam.correlated = false; %boolean, calculate (un)correlated displacements (false: more data points, true: independent displacements, see Trajectory.m)
fitParam.distanceMethod = 'disp'; %string, choose between fitting displacement x,y 'disp' or jump distance sqrt(x.^2+y.^2) 'jump'
fitParam.histogramBinWidth = 0; %double, manual override of histogram bin width. Choose 0 for automated mode (recommended)
fitParam.fitVelocity = false; %boolean, enables/disables velocity/flow fitting
fitParam.isotropic = true; %boolean, set to false if diffusion is non-isotropic. Automatically false if vswitch = true
fitParam.plotFit = false; %boolean, plot switch for fitted histograms
fitParam.breakCondition = [24,9]; %integer, stop fitting if [# of displacements, # of bins] < [m,n]
fitParam.fitAlpha = 0.0; %double, accept higher order if rejecting null hypothesis at confidence level 1-alpha
fitParam.weightedFit = true; %boolean, take weights into account when MSD curve is fitted
fitParam.fitTrueMSDValues = false; %boolean, when fitting an MSD curve e.g. from displacement histograms and multiple species are allowed, only consider MSD values which are the best outcome of a multivariate histogram fit.

%Printing in latex
printParam.latex = false; %boolean, print Latex figure afterwards, requires matlabfrag and ps-tools
printParam.size = [8,6,17,8]; %double, width and height of printed msd and histogram figure
printParam.filename = 'F:\filename'; %string, filename of printed picture
