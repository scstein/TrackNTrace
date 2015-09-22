function plot_drift_correction(trajectoryData, driftCorr_trajectoryData, drift, pixSize )
% Plots the original localizations in comparison to drift corrected data
% and (optionally) the drift.
%
%   Input:
%         trajectoryData -  2D double array, list of trajectories with
%                           columns [id,frame,xpos,ypos, ... ]. 
%                           (The first possible frame must be 1 not 0)
%         driftCorr_trajectoryData -  drift corrected trajectories 2D double
%                           array, list of trajectories with columns 
%                           [id,frame,xpos,ypos, ... ]. 
%                           (The first possible frame must be 1 not 0)
%         drift - The computed drift with columns [frame, x,y]. (can be left empty)
%         pixSize - pixel size in nanometers (used for labeling axis if given)
%
% Author: Simon Christoph Stein

if nargin<3 || isempty(drift)
    drift = [];
end

if nargin<4 || isempty(pixSize)
   pixSize = 1; 
   unit = ' [pix]';
else
   unit = ' [nm]';
end

% Plot the original and corrected datasets
h = figure;

hs1 = subplot(1,2,1);
plot(trajectoryData(:,3)*pixSize, trajectoryData(:,4)*pixSize, '.','markersize',1);
axis image;
axis ij; % flip y-axis upside down
title(hs1,'original data');
xlabel(['xPos',unit], 'FontSize',14);
ylabel(['yPos',unit], 'FontSize',14);
set(gca,'FontSize',14)

hs2 = subplot(1,2,2);
plot(driftCorr_trajectoryData(:,3)*pixSize, driftCorr_trajectoryData(:,4)*pixSize, '.','markersize',1);
axis image;
axis ij; % flip y-axis upside down
title(hs2,'corrected data');
xlabel(['xPos',unit], 'FontSize',14);
ylabel(['yPos',unit], 'FontSize',14);
set(gca,'FontSize',14)

linkaxes([hs1,hs2],'xy');

tightfig; % Remove space around plot


% Plot the drift
if ~isempty(drift)
    figure;
    plot(drift(:,2)*pixSize, drift(:,3)*pixSize);
    title('Drift with respect to first frame');
    xlabel(['xPos',unit], 'FontSize',14);
    ylabel(['yPos',unit], 'FontSize',14);
    set(gca,'FontSize',14)
    axis image;
    axis ij; % flip y-axis upside down
end

end

