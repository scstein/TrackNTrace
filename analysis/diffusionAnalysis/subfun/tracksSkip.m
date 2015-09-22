function [traj] = tracksSkip(tracks_processed,delta_f,use_corr)
% [traj] = tracksSkip(tracks_processed,delta_f,use_corr) Create array of
% displacements from trajectories for a timepoint of Delta_t = delta_f
% frames.
% 
% INPUT:
%     tracks_processed: [Nxd] double array of trajectories, separated by
%     one NaN row each. Currently, only d=2 is considered as we're handling
%     displacements here.
%     
%     delta_f: Double, calculate displacement for a time shift of delta_f
%     frames.
%     
%     use_corr: If particle position correlation is considered, only
%     independet displacements are calculated with obviously much fewer
%     data points.
%     
%     
% OUTPUT:
%     traj: [Nxd] double array of displacements, d=2.

traj_start_idx = [1;find(isnan(tracks_processed(:,1)))+1]; %index of start of trajectories
nrTraj = size(traj_start_idx,1)-1;
traj = cell(nrTraj,1);

if ~use_corr
    for iTraj=1:nrTraj
        traj_temp = tracks_processed(traj_start_idx(iTraj):traj_start_idx(iTraj+1)-2,:);
        traj_temp_length = size(traj_temp,1);
        
        if traj_temp_length>1
            traj(iTraj) = {traj_temp(1+delta_f:1:end,1:2)-traj_temp(1:1:end-delta_f,1:2)};
%             traj(iTraj) = {traj_temp(1+delta_f:1+delta_f:end,1:2)-traj_temp(1:1+delta_f:end-delta_f,1:2)};
        end
    end
    
else
    for iTraj=1:nrTraj
        traj_temp = tracks_processed(traj_start_idx(iTraj):traj_start_idx(iTraj+1)-2,:);
        traj_temp_length = size(traj_temp,1);
        
        if traj_temp_length>1
            traj(iTraj) = {traj_temp(1+delta_f:1+delta_f:end,1:2)-traj_temp(1:1+delta_f:end-delta_f,1:2)};
%             traj(iTraj) = {traj_temp(1+delta_f:end,1:2)-traj_temp(1:end-delta_f,1:2)};
        end
    end
end

traj = vertcat(traj{:});
