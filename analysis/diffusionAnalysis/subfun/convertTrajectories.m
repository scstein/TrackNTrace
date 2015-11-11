function [traj_result] = convertTrajectories(traj,trajParam,trmethod)
% [traj_traj_result,pnum] = convertTrajectories(traj,trajParam,trmethod)
% Convert trajectories from several tracking engines to format suitable for
% diffusionAnalysis
% 
% INPUT:
%     traj: Can be filename of file where trajectory information is stored
%     or trajectory variable itself.
%     
%     trajParam: Struct of options used to parse trajectories, see
%     RunDiffusionAnalysis.m for details.
%     
%     trmethod: String, name of tracking engine. Currently supported:
%     TrackNTrace 'tnt', MOSAIC (through Fiji) 'mosaic', Icy 'icytrack',
%     and u-Track 'utrack'.
%  
%     
% OUTPUT:
%     traj_result: [Nx1] cell array of trajectories where N is the number
%     of trajectories. Every trajectory is divided into segments, row
%     arrays of continously recorded particles, which are separated by an
%     NaN line.

switch lower(trmethod)
    
    case('tnt')
        if ischar(traj)
            traj = load(traj,'-mat','trajectoryData');
            traj = traj.trajectoryData;
        end
        
        nrTraj = max(traj(:,1)); %number of trajectories
        traj_result = cell(nrTraj,1);
        
        for iTraj=1:nrTraj
            traj_temp = traj(traj(:,1)==iTraj,2:4); %only positions -> this could be changed
            segment_end_frames = [0;find(traj_temp(2:end,1)-traj_temp(1:end-1,1)-1);size(traj_temp,1)]; %contains end of all trajectory segments
            
            for jSeg=1:numel(segment_end_frames)-1
                seg_temp = traj_temp(1+segment_end_frames(jSeg):segment_end_frames(jSeg+1),2:3);
                seg_size = size(seg_temp,1);
                
                sigma_guess = mean([std(seg_temp(:,1)),std(seg_temp(:,2))]); %rough guess for per-frame diffusion
                
                is_in_zone = true;
                if ~isempty(trajParam.inclusionZone) %check if particle are outside of considered circular region and discard it if applicable
                    is_in_zone = sum(sqrt(sum((seg_temp-repmat(trajParam.inclusionZone(1:2),seg_size,1)).^2,2))<=trajParam.inclusionZone(3))==seg_size;
                end
                
                %check if particle is moving, if the segment has the right length, if it's outside the inclusion zone
                if sigma_guess>=trajParam.minimalDiffusionGuess && seg_size>1 && seg_size<=trajParam.trajMaxLength && is_in_zone
                    traj_result(iTraj) = {[traj_result{iTraj};seg_temp;NaN(1,2)]};
                end
            end
        end
 
        
    case('mosaic') %Sbalzarini Fiji particle tracker, x,y,I (also pattern, but not loaded)
        if ischar(traj)
            try
                traj = dlmread(traj,'\t',1,0);
                traj = traj(:,2:5);
            catch
                traj = read_in_traj(traj,1); %MEX file
            end
        end
        traj = traj(:,1:4);
        traj(:,2) = traj(:,2)+1; %mosaic starts at t=0. Correct to t=1 frame
        
        [traj_result] = convertTrajectories(traj,trajParam,'TNT');
        
        
    case('icytrack') %ICY spot detector+probabalistic particle tracker
        traj = xlsread(traj);
        
        idx = find(~isnan(traj(:,1)))+2;
        idx(1:end-1,2) = idx(2:end,1)-4; idx(end,2) = size(traj,1);
        nrTraj = size(idx,1);
        
        traj_result = cell(nrTraj,1);
        for iTraj=1:nrTraj
            seg_temp = traj(idx(iTraj,1):idx(iTraj,2),3:4);
            
            D_guess = mean([std(seg_temp(:,1)),std(seg_temp(:,2))]);
            seg_size = size(seg_temp,1);
            
            is_in_zone = true;
            if ~isempty(trajParam.inclusionZone)
                is_in_zone = sum(sqrt(sum((seg_temp-repmat(trajParam.inclusionZone(1:2),seg_size,1)).^2,2))<=trajParam.inclusionZone(3))==seg_size;
            end
            
            if D_guess>=trajParam.minimalDiffusionGuess && seg_size>1 && seg_size<=trajParam.trajMaxLength && is_in_zone
                traj_result(iTraj) = {[traj_result{iTraj};seg_temp;NaN(1,2)]};
            end
            
        end
        
        
    case('utrack') %u-track position+tracking, x,y,I,sigma_x,sigma_y,radius of gyration, ...
        traj = load(traj,'-mat','tracksFinal'); traj = traj.tracksFinal;
        
        nrTraj = numel(traj);
        %seqofEvents: has a 0 where track has gap (NaN in
        %tracksCoordAmp) tracksCoordAmpCG:
        %[x,y,z,amp,xerr,yerr,zerr,amperr,x,y,...] tracksFeatIndxCG:
        %contains start and end frame
        
        traj_result = cell(nrTraj,1);
        for iTraj = 1:nrTraj
            nrFrames = size(tracksFinal(iTraj).tracksFeatIndxCG,2); %number of frames
            frameStart = tracksFinal(iTraj).seqOfEvents(1,1);
            xyamp = [tracksFinal(iTraj).tracksCoordAmpCG(1:8:end).',tracksFinal(iTraj).tracksCoordAmpCG(2:8:end).',tracksFinal(iTraj).tracksCoordAmpCG(4:8:end).']; %[x,y,amp]
            
            
            traj_temp = [(frameStart:frameStart+nrFrames-1).', xyamp]; traj_temp = traj_temp(~isnan(traj_temp(:,2)),:); %[frame,x,y,amp]
            segment_end_frames = [0;find(traj_temp(2:end,1)-traj_temp(1:end-1,1)-1);size(traj_temp,1)]; %contains end of all trajectory segments
            for jSeg=1:numel(segment_end_frames)-1
                seg_temp = traj_temp(1+segment_end_frames(jSeg):segment_end_frames(jSeg+1),2:3);
                D_guess = mean([std(seg_temp(:,1)),std(seg_temp(:,2))]);
                seg_size = size(seg_temp,1);
                
                is_in_zone = true;
                if ~isempty(trajParam.inclusionZone)
                    is_in_zone = sum(sqrt(sum((seg_temp-repmat(trajParam.inclusionZone(1:2),seg_size,1)).^2,2))<=trajParam.inclusionZone(3))==seg_size;
                end
                
                if D_guess>=trajParam.minimalDiffusionGuess && seg_size>1 && seg_size<=trajParam.trajMaxLength && is_in_zone
                    traj_result(iTraj) = {[traj_result{iTraj};seg_temp;NaN(1,2)]};
                end
            end
        end

        
    otherwise
        fprintf('No method selected or unknown tracker, aborting. \n');
        traj_result = [];
        return
end
