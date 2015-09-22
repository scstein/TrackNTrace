function [msd_result,use_v] = diffusionAnalysisMSD(tracks_processed,fitParam,experimParam)
% [msd_result,use_v] = diffusionAnalysisMSD(tracks_processed,fitParam,experimParam)
% Calculate MSD values for trajectories according to setup in RunDiffusionAnalysis.m
%
% INPUT:
%     tracks_processed: Cell row array of trajectories with one trajectory
%     row array per cell. Trajectory segments have to be separated by a NaN
%     line.
%
%     fitParam: Struct of fitting options, see RunDiffusion Analysis.m for
%     details.
%
%     experimParam: Struct of experimental parameters, see RunDiffusion
%     Analysis.m for details.
%
%
% OUTPUT:
%     msd_result: Cell column array of msd curves. Every cell contains a
%     double array [t,sigma,error_sigma,weight,error_weight,outcome]. Here,
%     t is the time in [s], sigma are the MSD values of each diffusing
%     subspecies in [µm^2] and weight are the weights (1 and 0 respectively
%     as one cannot consider subspecies for
%
%     The column marks the dimension.
%
%     use_v: Switch for fitting velocity later. If input options don't
%     match, fitting the velocity will be disabled.


% Parse initial parameters
px = experimParam.pixel/experimParam.magnification;
dt = experimParam.timestep;
max_frame = fitParam.maximumSkip;
use_corr = fitParam.correlated;
fit_curve = fitParam.distanceMethod;
use_v = fitParam.fitVelocity;
use_iso = fitParam.isotropic;

if ~use_iso && strcmp(fit_curve,'jump');
    fprintf('Cannot consider non-isotropicity for jump distance MSD. Aborting. \n');
    return;
end

if use_v
    if strcmp(fit_curve,'jump')
        fprintf('Cannot fit velocity with jump distance histogram. Disabling velocity fit. \n');
        use_v = false;
    end
    
    if use_iso
        fprintf('Cannot fit velocity for one dimension. Disabling velocity fit. \n');
        use_v = false;
    end
end

dim = 1+~use_iso;

msd_result = cell(1,dim);
for iFrame=1:max_frame
    
    %format all trajectories with appropriate skip frame length
    [traj] = tracksSkip(tracks_processed,iFrame,use_corr);
    
    if strcmp(fit_curve,'disp')
        if use_iso
            traj = traj(:);
        end
        for jDim = 1:dim
            msd_result(1,jDim) = {[msd_result{1,jDim};[iFrame*dt, mean(traj(:,jDim).^2)*px^2,std(traj(:,jDim).^2)*px^2,[1,0],true]]}; %different species cannot be considered here but weights are added nonetheless to keep code compatible
        end
    else
        msd_result(1) = {[msd_result{1};[iFrame*dt, mean(traj(:,1).^2+traj(:,2).^2)*px^2,std(traj(:,1).^2+traj(:,2).^2)*px^2,[1,0],true]]};
    end
end
