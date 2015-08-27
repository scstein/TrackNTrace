function [trajData] = trackParticles(pos_file,trackingOptions)
% trajData = trackParticles(posFilename,trackingOptions)
% This function loads particle positions within movies obtained by the
% TrackNTrace routine "locateParticles" and passes them to a particle
% tracking rountine to obtain trajectories. 
% Supported trackers are: simpletracker, uTrack, track.m by Crocker and
% Grier
% 
% The position file is loaded and an array suitable for the respective
% tracker is created. If needed, the file is split up into several parts
% and the parts are tracked separately to prevent memory and performance
% problems. Afterwards, the parts are spliced together at the seam by
% tracking the two frames adjacent to the split "border"
% 
% INPUT:
%     posFilename: string, Full path of .mat file where result of
%     locateParticles routine is stored trackingOptions: struct, options for
%     trackers
%     
% OUTPUT:
%     trajData: 2D double array, list of trajectories in a suitable format
%     [id,frame,xpos,ypos,amp]. Every trajectory is given an id, starting
%     at 1, after which the list is sorted. Frame number starts with 1,
%     positions are given in pixels and the amplitude is the gaussian
%     amplitude value A in A*exp(...) given by locate Particles

%Parse inputs
method = trackingOptions.method;
probDim = trackingOptions.probDim;
r_max = trackingOptions.maxRadius;
gap = trackingOptions.maxGap;
tlen_min = trackingOptions.minTrackLength;
n_split = trackingOptions.splitMovieParts;
verbose = trackingOptions.verbose;

%Read positions file and convert to array appropriate for respective
%tracker
[pos] = convertPositions(pos_file,method,probDim);

traj_id = 0; %global traj_id, keeps track of every slice

trajData = [];


if ~verbose
    fprintf('Tracking particles using %s ..\n',method);
end

switch method
    case 'simpletracker'
        n_frames = size(pos,1); %number of frames
        stack_slice = round(n_frames/n_split); %size of slice if tracking is divided

        %collect complete list of positions with amplitude and frame as
        %simpletracker does not return either
        global_pos = convertPositions(pos_file,'track_cg',5); %automatically gives 2D array of global positions [x y amp b sigma frame]
        
        
        for iDiv = 1:n_split %slice position array if memory not large enough
            slice = [1+(iDiv-1)*stack_slice,min(iDiv*stack_slice,n_frames)]; %start and end of this slice
            global_pos_slice = global_pos(global_pos(:,end)>=slice(1)&global_pos(:,end)<=slice(2),:);
            
            %call tracker
            [~, adjacency_tracks] = simpletracker(pos(slice(1):slice(2)), 'Method', trackingOptions.linkingMatrix, 'MaxLinkingDistance',r_max, 'MaxGapClosing',gap,'Debug',verbose);
            %remove trajectories deemed too short
            adjacency_tracks = adjacency_tracks(find(cellfun(@(var) numel(var)>(tlen_min-1),adjacency_tracks))); %#ok<FNDSB>
            n_tracks = numel(adjacency_tracks);
            
            cell_coord_tracks = cell(n_tracks, 1);
            for iTrack = 1:n_tracks
                track = adjacency_tracks{iTrack};
                track_points = global_pos_slice(track,:);
                track_len = size(track_points,1);
                
                %index with relative id but write global id
                cell_coord_tracks{iTrack} = [repmat(traj_id+iTrack,track_len,1),track_points(:,[end,1:5])]; %[id, [frame,x,y,amp,b,sigma]]
            end
            
            %finally save data
            trajData = [trajData; vertcat(cell_coord_tracks{:})]; %#ok<AGROW>
            
            %connect slices, it's impossible to consider gaps between
            %slices, so only connect adjacent frames
            if iDiv>1
                %get trajectories in frames adjacent to slice border
                cell_track_slice = [{trajData(trajData(:,2)==slice(1)-1,:)};{trajData(trajData(:,2)==slice(1),:)}];
                track_slice = vertcat(cell_track_slice{:});
                
                %track with same settings. gaps cannot happen
                if ~isempty(cell_track_slice{1}) && ~isempty(cell_track_slice{2})
                    [~, adjacency_tracks] = simpletracker([{cell_track_slice{1}(:,(3:3+probDim-1))},{cell_track_slice{2}(:,(3:3+probDim-1))}], 'Method', 'Hungarian', 'MaxLinkingDistance',r_max, 'MaxGapClosing',gap);
                    adjacency_tracks = adjacency_tracks(find(cellfun(@(var) numel(var)>1,adjacency_tracks))); %#ok<FNDSB>
                    n_tracks_sliced = numel(adjacency_tracks);
                else
                    n_tracks_sliced = 0;
                end
                
                if n_tracks_sliced>0
                    for iTrack = 1:n_tracks_sliced
                        track = adjacency_tracks{iTrack};
                        
                        id_pair = track_slice(track,1); %[old id;new id]
                        
                        connected_track_idx = (trajData(:,1)==id_pair(2)); %get all entries of new trajectory to be connected to older one
                        trajData(connected_track_idx,1) = repmat(id_pair(1),sum(connected_track_idx),1); %set to old id
                    end
                    
                    %some of the ids were reverted to the smaller values of
                    %trajectories in old slice. all unconnected or new
                    %tracks will have an id which is too high. we'll
                    %correct this now
                    traj_update = 0;
                    for jId = traj_id+1:traj_id+n_tracks %go through all tracks in newly added slice
                        idx_trajData_update = trajData(:,1)==jId;
                        n_to_update = sum(idx_trajData_update);
                        
                        if n_to_update>0 %is this an unconnected or new track?
                            trajData(idx_trajData_update,1) = repmat(traj_id+traj_update+1,n_to_update,1); %then set back the id...
                            traj_update = traj_update+1; %increment the updater id...
                        end
                    end
                    
                    traj_id = traj_id+traj_update; %...and finally update the global id
                else
                    traj_id = traj_id+n_tracks; %we are not in the first slice and there's nothing to reconnect. it's time to update the global id
                end %if n_tracks_sliced>0
            else 
                traj_id = traj_id+n_tracks; %if in very first slice, there's nothing to reconnect. then just update the global id
            end %if iDiv>1
        end %for iDiv=1:n_split
        
        % if the movie was split, result array has to be sorted by
        % trajectory id again
        if n_split>1
            [~,idx] = sort(trajData(:,1));
            trajData = trajData(idx,:);
        end
        %FINISHED simpletracker
        
        
    case 'utrack'
        
        n_frames = size(pos,1); %number of frames
        stack_slice = round(n_frames/n_split); %size of slice if tracking is divided
        
        %get option structs for utrack
        [gapCloseParam,costMatrices,kalmanFunctions] = parseUtrackOptions(trackingOptions);
        trackingOptions_slice = trackingOptions; trackingOptions_slice.minTrackLength = 2;
        [gapCloseParam_slice,costMatrices_slice,kalmanFunctions_slice] = parseUtrackOptions(trackingOptions_slice);
        
        for iDiv = 1:n_split %slice position array if memory not large enough
            slice = [1+(iDiv-1)*stack_slice,min(iDiv*stack_slice,n_frames)]; %start and end of this slice
            
            %call utrack with probDim=2 as it can only deal with spatial
            %probDims plus amplitude. we only track 2D here
            [tracksFinal,~,~] = trackCloseGapsKalmanSparse(pos(slice(1):slice(2)),costMatrices,gapCloseParam,kalmanFunctions,2,0,verbose);
            n_tracks = numel(tracksFinal);
            %seqofEvents: has a 0 where track has gap (NaN in
            %tracksCoordAmp) tracksCoordAmpCG:
            %[x,y,z,amp,xerr,yerr,zerr,amperr,x,y,...] tracksFeatIndxCG:
            %contains start and end frame
            
            cell_coord_tracks = cell(n_tracks,1);
            for iTrack = 1:n_tracks
                nFrames = size(tracksFinal(iTrack).tracksFeatIndxCG,2); %number of frames
                frameStart = tracksFinal(iTrack).seqOfEvents(1,1)+(iDiv-1)*stack_slice; %get relative frame and correct for global frame
                xyamp = [tracksFinal(iTrack).tracksCoordAmpCG(1:8:end).',tracksFinal(iTrack).tracksCoordAmpCG(2:8:end).',tracksFinal(iTrack).tracksCoordAmpCG(4:8:end).']; %[x,y,amp]
                
                cell_coord_tracks(iTrack) = {[repmat(traj_id+iTrack,nFrames,1), (frameStart:frameStart+nFrames-1).', xyamp]}; %[id,frame,x,y,amp]
            end
            
            %finally save data
            trajData = [trajData; vertcat(cell_coord_tracks{:})]; %#ok<AGROW>
            
            %connect slices, it's impossible to consider gaps between
            %slices, so only connect adjacent frames
            if iDiv>1
                %get trajectories in frames adjacent to slice border
                track_slice = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),2,1); %1D struct array with nrFrames lines, inner arrays have two columns [value,error]
                for iFrame=1:2
                    pos_frame_now = trajData(trajData(:,2)==(slice(1)-2+iFrame),:);
                    nCand = size(pos_frame_now,1);
                    track_slice(iFrame).xCoord = [pos_frame_now(:,3),zeros(nCand,1)]; %careful, check if error should be >0!
                    track_slice(iFrame).yCoord = [pos_frame_now(:,4),zeros(nCand,1)];
                    track_slice(iFrame).amp = [pos_frame_now(:,5),zeros(nCand,1)];
                    if iFrame==1
                        pos_frame_first = pos_frame_now; %need this for getting trajectory ids back
                    end
                end
                
                %track those two frames
                if ~isempty(pos_frame_now) && ~isempty(pos_frame_first)
                    try
                        [tracksFinal_sliced,~,~] = trackCloseGapsKalmanSparse(track_slice,costMatrices_slice,gapCloseParam_slice,kalmanFunctions_slice,2,0,verbose);
                        n_tracks_sliced = numel(tracksFinal_sliced);
                    catch
                        %very rarely, utrack tries to index a feature track
                        %matrix with a wrong relaive index starting at 0
                        n_tracks_sliced = 0;
                    end
                else
                    n_tracks_sliced = 0;
                end
                
                if n_tracks_sliced>0
                    for iTrack = 1:n_tracks_sliced
                        x_pos = tracksFinal_sliced(iTrack).tracksCoordAmpCG(1:8:end).'; %[x_old,x_new]
                        
                        %we cannot get back the respective ids, so we have
                        %to search for positions instead (risky - use fuzzy
                        %equal)
                        id_pair = [pos_frame_first(abs(pos_frame_first(:,3)-x_pos(1))<1e-6,1);pos_frame_now(abs(pos_frame_now(:,3)-x_pos(2))<1e-6,1)];
                        
                        %now correct ids
                        connected_track_idx = (trajData(:,1)==id_pair(2)); %get all entries of new trajectory to be connected to older one
                        trajData(connected_track_idx,1) = repmat(id_pair(1),sum(connected_track_idx),1);
                    end
                    
                    %some of the ids were reverted to the smaller values of
                    %trajectories in old slice. all unconnected or new
                    %tracks will have an id which is too high. we'll
                    %correct this now
                    traj_update = 0;
                    for jId = traj_id+1:traj_id+n_tracks %go through all tracks in newly added slice
                        idx_trajData_update = trajData(:,1)==jId;
                        n_to_update = sum(idx_trajData_update); %is this an unconnected or new track?
                        if n_to_update>0
                            trajData(idx_trajData_update,1) = repmat(traj_id+traj_update+1,n_to_update,1); %then set back the id...
                            traj_update = traj_update+1; %increment the updater id...
                        end
                    end
                    traj_id = traj_id+traj_update; %...and finally update the global id
                else
                    traj_id = traj_id+n_tracks; %we are not in the first slice and there's nothing to reconnect. it's time to update the global id
                end %if n_tracks_sliced>0
            else 
                traj_id = traj_id+n_tracks; %if in very first slice, there's nothing to reconnect. then just update the global id
            end %if iDiv>1
        end %for iDiv=1:n_split
        
        %remove NaNs resulting from closed gaps
        trajData = trajData(~isnan(trajData(:,end)),:);
        
        % if the movie was split, result array has to be sorted by
        % trajectory id again
        if n_split>1
            [~,idx] = sort(trajData(:,1));
            trajData = trajData(idx,:);
        end
        %FINISHED utrack
    
        
    case 'track_cg'
        
        n_frames = max(pos(:,end)); %number of frames
        stack_slice = round(n_frames/n_split); %size of slice if tracking is divided
        
        %Prepare option struct
        track_cgParam.mem = gap-1;
        track_cgParam.dim = probDim;
        track_cgParam.good = tlen_min;
        track_cgParam.quiet = ~verbose;
        
        pos = pos(:,[1:1+probDim-1,end]); %track_cg cannot handle an array larger than it's supposed to be juding from the probDim
        
        
        for iDiv = 1:n_split %slice position array if memory not large enough
            slice = [1+(iDiv-1)*stack_slice,min(iDiv*stack_slice,n_frames)]; %start and end of this slice

            %call tracker
            tracks_cg = track_cg(pos(pos(:,end)>=slice(1)&pos(:,end)<=slice(2),:),r_max,track_cgParam); %[x,y,...,t,id], user-friendly format!
            n_tracks = max(tracks_cg(:,end));
            
            %save data, updating relative id with global id
            trajData = [trajData; [tracks_cg(:,end)+traj_id,tracks_cg(:,[end-1,1:probDim])]]; %#ok<AGROW>
            
            %connect slices, it's impossible to consider gaps between
            %slices, so only connect adjacent frames
            if iDiv>1
                %get trajectories in frames adjacent to slice border
                pos_frame_first = trajData(trajData(:,2)==(slice(1)-1),:);
                pos_frame_now = trajData(trajData(:,2)==slice(1),:);
                
                %track those two frames
                try %track_cg is sometimes terrible at tracking so no frames remain. in this case, it throws an error message which is also terrible
                    tracks_cg_sliced = track_cg([pos_frame_first(:,[3:end,2]);pos_frame_now(:,[3:end,2])],r_max,track_cgParam);
                    n_tracks_sliced = max(tracks_cg_sliced(:,4));
                catch
                    n_tracks_sliced=0;
                end
                
                
                if n_tracks_sliced>0
                    for iTrack = 1:n_tracks_sliced
                        %we cannot get back the respective ids, so we have
                        %to search for positions instead (risky - use fuzzy
                        %equal)
                        x_pos = tracks_cg_sliced(find(tracks_cg_sliced(:,4)==iTrack,2),1);
                        id_pair = [pos_frame_first(abs(pos_frame_first(:,3)-x_pos(1))<1e-6,1);pos_frame_now(abs(pos_frame_now(:,3)-x_pos(2))<1e-6,1)];
                        
                        %now correct ids
                        connected_track_idx = (trajData(:,1)==id_pair(2)); %get all entries of new trajectory to be connected to older one
                        trajData(connected_track_idx,1) = repmat(id_pair(1),sum(connected_track_idx),1);
                    end
                    
                    %some of the ids were reverted to the smaller values of
                    %trajectories in old slice. all unconnected or new
                    %tracks will have an id which is too high. we'll
                    %correct this now
                    traj_update = 0;
                    for jId = traj_id+1:traj_id+n_tracks %go through all tracks in newly added slice
                        idx_trajData_update = trajData(:,1)==jId;
                        n_to_update = sum(idx_trajData_update); %is this an unconnected or new track?
                        if n_to_update>0
                            trajData(idx_trajData_update,1) = repmat(traj_id+traj_update+1,n_to_update,1); %then set back the id...
                            traj_update = traj_update+1; %increment the updater id...
                        end
                    end
                    traj_id = traj_id+traj_update; %...and finally update the global id
                else
                    traj_id = traj_id+n_tracks; %we are not in the first slice and there's nothing to reconnect. it's time to update the global id
                end %if n_tracks_sliced>0
            else 
                traj_id = traj_id+n_tracks; %if in very first slice, there's nothing to reconnect. then just update the global id
            end %if iDiv>1
        end %for iDiv=1:n_split
        
        % if the movie was split, result array has to be sorted by
        % trajectory id again
        if n_split>1
            [~,idx] = sort(trajData(:,1));
            trajData = trajData(idx,:);
        end
        %FINISHED track_cg
        
    case('nn_cpp')       
        trajData = nn_tracker_cpp(pos,tlen_min,r_max,r_max,gap,[],verbose).';
        
    otherwise
        error('Unknown tracker ''%s''',method);
end

if ~verbose
    fprintf('\b done\n');
end

end


% function [trajData,traj_update] = correctTrajIds(trajData,traj_id,n_tracks)
% traj_update = 0;
% for jId = traj_id+1:traj_id+n_tracks %go through all tracks in newly added slice
%     idx_trajData_update = trajData(:,1)==jId;
%     n_to_update = sum(idx_trajData_update); %is this an unconnected or new track?
%     if n_to_update>0
%         trajData(idx_trajData_update,1) = repmat(traj_id+traj_update+1,n_to_update,1); %then set back the id...
%         traj_update = traj_update+1; %increment the updater id...
%     end
% end
% end



