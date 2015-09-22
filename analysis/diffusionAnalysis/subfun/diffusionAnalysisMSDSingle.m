function [D_final,velocity_final] = diffusionAnalysisMSDSingle(tracks_processed,fitParam,experimParam)
% [D_final,velocity_final] = diffusionAnalysisMSDSingle(tracks_processed,fitParam,experimParam)
% Obtain histogram of diffusion and velocity values from single-trajectory MSD fits.
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
%     D_final: Struct of histogram parameters for diffusion distribution.
%         .barBins: Cell array of bins (x-axis), one cell per dimension.
%         .barFreq: Cell array of frequencies (y-axis), one cell per
%         dimension.
%         
%     velocity_final: Struct of histogram paramters for velocity
%     distribution.


% Parse initial parameters
px = experimParam.pixel/experimParam.magnification;
dt = experimParam.timestep;
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


traj_start_idx = [1;find(isnan(tracks_processed(:,1)))+1]; %index of start of trajectories
nrTraj = size(traj_start_idx,1)-1;
data_D = zeros(nrTraj,2*dim);

if use_v
    data_v = zeros(nrTraj,2*dim);
end

use_disp = strcmp(fit_curve,'disp');


for iTraj=1:nrTraj
    traj_temp = tracks_processed(traj_start_idx(iTraj):traj_start_idx(iTraj+1)-2,1:2);
    
    if use_v
        [data_D(iTraj,:),data_v(iTraj,:)] = fitMSDSingle(traj_temp,use_v,use_iso,use_disp,fitParam);
    else
        [data_D(iTraj,:),~] = fitMSDSingle(traj_temp,use_v,use_iso,use_disp,fitParam);
    end
end

data_D = checkOutput(data_D)*px^2/dt;
bins_out_D = cell(1,dim);
hist_out_D = cell(1,dim);

if use_v
    data_v = checkOutput(data_v)*px/dt;
    bins_out_v = cell(1,dim);
    hist_out_v = cell(1,dim);
end


for jDim=1:dim
    [bins_out_D{jDim},hist_out_D{jDim},~] = Hist1D(data_D(:,1+(jDim-1)*2),0,false,'freedman',true,false,false);
    if use_v
        [bins_out_v{jDim},hist_out_v{jDim},~] = Hist1D(data_v(:,1+(jDim-1)*2),0,false,'freedman',true,false,false);
    end
end

if fitPara.plotFit
    for jDim=1:dim
        figure;
        bar(bins_out_D{jDim},hist_out_D{jDim});
        xlabel('D [µm^2/s]');
        ylabel('Frequency');
        if jDim==1
            if use_iso
                title_string = 'Diffusion histogram, one dimension.';
            else
                title_string = 'Diffusion histogram, x-dimension.';
            end
        else
            title_string = 'Diffusion histogram, y-dimension.';
        end
        title(title_string);
        
        if use_v
            figure;
            bar(bins_out_v{jDim},hist_out_v{jDim});
            xlabel('v [µm/s]');
            ylabel('Frequency');
            if jDim==1
                title_string = 'Diffusion histogram, x-dimension.';
            else
                title_string = 'Diffusion histogram, y-dimension.';
            end
            title(title_string);
        end
    end
    
    
end

D_final.barBins = bins_out_D;
D_final.barFreq = hist_out_D;
velocity_final.barBins = bins_out_v;
velocity_final.barFreq = hist_out_v;


end %END MAINFUN


function [data] = checkOutput(data)
data = data(any(data(:,1),2),:);
data = data(~isinf(sum(data,2)),:);

if size(data,2)/2==1
    data = data(data(:,1)>0);
else
    idx = data(:,1)>0 & data(:,3)>0;
    data = data(idx,:);
end
end

