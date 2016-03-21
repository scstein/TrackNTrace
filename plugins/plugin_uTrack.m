function [plugin] = plugin_uTrack()

%    -------------- Definition of plugin --------------
% Name of the component these options are for
name = 'u-Track';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 3;

% The functions this plugin implements
mainFunc =  @trackParticles_uTrack;

% Description of output parameters
outParamDescription = {'Track-ID';'Frame';'x';'y';'z';'Amplitude'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Track particles using u-Track. \n\n',...
    'u-Track is a very complex, memory intensive particle tracking suite with a vast array of options of which only a small part can be set here in the GUI. Advanced users can edit the file ''plugin_uTrack'' to access all options. \nu-Track was published in ''Jaqaman et al, NatMet 5(8), 695-702 (2008), doi:10.1038/nmeth.1237''.'];

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
plugin.add_param('minTrajLength',...
    'int',...
    {2, 1, inf},...
    'Minimum length of trajectories in [frames].');
plugin.add_param('maxTrackRadius',...
    'float',...
    {6, 0, inf},...
    'Maximum allowed linking distance between two localizations in [pixels].');
plugin.add_param('maxFrameGap',...
    'int',...
    {0, 0, inf},...
    'Maximum allowed time gap between two segments used for gap closing in [frames]. \nWhen fluorophores blink or disappear, they obviuosly cannot be tracked anymore. If they reappear, however, tracked segments can be stitched back together. \n0 = no gap closing.');
plugin.add_param('splitMovieIntervals',...
    'int',...
    {5, 0, inf},...
    'Number of slices to divide movie into, useful when processing large number of particles and frames. \nTracks are stitched back together at the split frames.');
plugin.add_param('track3D',...
    'bool',...
    false,...
    'Enable 3D tracking.');
plugin.add_param('verbose',...
    'bool',...
    false,...
    'Switch on to see tracking progress in command window.');
end


%   -------------- User functions --------------

function [trackingData] = trackParticles_uTrack(fittingData,options)
% Wrapper function for u-Track (see below). Refer to tooltips above, to
% parseUtrackOptions function and to u-Track manual to obtain information
% on input and output variables.
% u-Track was programmed in the lab of Gaudenz Danuser, see:
% Jaqaman et al, Nature Methods - 5 (8), 695 - 702 (2008), doi:10.1038/nmeth.1237
% Refer to GPL-License.txt for licensing information. u-Track can be downloaded here:
% http://lccb.hms.harvard.edu/software.html
% To install, extract all files into 'external' folder
%
% INPUT:
%     fittingData: Cell array of localizations created by locateParticles.m
%     Refer to that function or to TrackNTrace manual for more information.
%
%     options: Struct of input parameters provided by GUI.
%
% OUTPUT:
%     trackingData: 2D double array of trajectories in format
%     [id,frame,xpos,ypos,zpos,amp]. Refer to trackParticles.m or to TrackNTrace
%     manual for more information.


% convert fittingData cell array to adjust to tracker function input

nrFrames = size(fittingData,1);
if nrFrames<200
    options.splitMovieIntervals=1; %no need to split small movies
end

pos = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'sigma',[]),nrFrames,1); %1D struct array with nrFrames lines, inner arrays have two columns [value,error]


for iFrame=1:nrFrames
    if(isempty(fittingData{iFrame})); continue; end; % Jump empty frames
    pos_frame_now = fittingData{iFrame};
    valid_pos = pos_frame_now(:,1)>=0; %error flag is 1?
    nCand = sum(valid_pos);
    pos_frame_now(pos_frame_now==0) = 1e-6; %this is a dirty hack for particles which run out of the frame
    
    pos(iFrame).xCoord = [pos_frame_now(valid_pos,1),zeros(nCand,1)]; %careful, check if error should be >0!
    pos(iFrame).yCoord = [pos_frame_now(valid_pos,2),zeros(nCand,1)];
    
    if options.track3D
        pos(iFrame).zCoord = [pos_frame_now(valid_pos,3),zeros(nCand,1)];
    else
        pos(iFrame).zCoord = [zeros(nCand,2)];
    end
    
    pos(iFrame).amp = [pos_frame_now(valid_pos,4),zeros(nCand,1)];
    pos(iFrame).sigma = [pos_frame_now(valid_pos,6),zeros(nCand,1)];
end


% call main function
[trackingData] = uTrackMain(pos,options);

end


function [trackingData] = uTrackMain(pos,trackingOptions)

trackingData = [];
traj_id = 0; %global trajectory idx

% probDim = 2;
probDim = 2+trackingOptions.track3D;
nrSplit = trackingOptions.splitMovieIntervals;
verbose = trackingOptions.verbose;
n_frames = size(pos,1); %number of frames
stack_slice = round(n_frames/nrSplit); %size of slice if tracking is divided

%get option structs for utrack
[gapCloseParam,costMatrices,kalmanFunctions] = parseUtrackOptions(trackingOptions);
trackingOptions_slice = trackingOptions; trackingOptions_slice.minTrackLength = 2;
[gapCloseParam_slice,costMatrices_slice,kalmanFunctions_slice] = parseUtrackOptions(trackingOptions_slice);

for iDiv = 1:nrSplit %slice position array if memory not large enough
    slice = [1+(iDiv-1)*stack_slice,min(iDiv*stack_slice,n_frames)]; %start and end of this slice
    
    %call utrack with probDim=2+track3D
    [tracksFinal,~,~] = trackCloseGapsKalmanSparse(pos(slice(1):slice(2)),costMatrices,gapCloseParam,kalmanFunctions,probDim,0,verbose);
    nrTracks = numel(tracksFinal);
    %seqofEvents: has a 0 where track has gap (NaN in
    %tracksCoordAmp) tracksCoordAmpCG:
    %[x,y,z,amp,xerr,yerr,zerr,amperr,x,y,...] tracksFeatIndxCG:
    %contains start and end frame
    
    cell_coord_tracks = cell(nrTracks,1);
    for iTrack = 1:nrTracks
        nrFrames = size(tracksFinal(iTrack).tracksFeatIndxCG,2); %number of frames
        frameStart = tracksFinal(iTrack).seqOfEvents(1,1)+(iDiv-1)*stack_slice; %get relative frame and correct for global frame
%         xyamp = [tracksFinal(iTrack).tracksCoordAmpCG(1:8:end).',tracksFinal(iTrack).tracksCoordAmpCG(2:8:end).',tracksFinal(iTrack).tracksCoordAmpCG(4:8:end).']; %[x,y,amp]
        xyzamp = [tracksFinal(iTrack).tracksCoordAmpCG(1:8:end).',tracksFinal(iTrack).tracksCoordAmpCG(2:8:end).',tracksFinal(iTrack).tracksCoordAmpCG(3:8:end).',tracksFinal(iTrack).tracksCoordAmpCG(4:8:end).']; %[x,y,amp]
        
%         cell_coord_tracks(iTrack) = {[repmat(traj_id+iTrack,nrFrames,1), (frameStart:frameStart+nrFrames-1).', xyamp]}; %[id,frame,x,y,amp]
        cell_coord_tracks(iTrack) = {[repmat(traj_id+iTrack,nrFrames,1), (frameStart:frameStart+nrFrames-1).', xyzamp]}; %[id,frame,x,y,z,amp]
    end
    
    %finally save data
    trackingData = [trackingData; vertcat(cell_coord_tracks{:})]; %#ok<AGROW>
    
    %connect slices; it's impossible to consider gaps between
    %slices, so only connect adjacent frames
    if iDiv>1
        %get trajectories in frames adjacent to slice border
%         track_slice = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),2,1); %1D struct array with nrFrames lines, inner arrays have two columns [value,error]
        track_slice = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),2,1); %1D struct array with nrFrames lines, inner arrays have two columns [value,error]
        for iFrame=1:2
            pos_frame_now = trackingData(trackingData(:,2)==(slice(1)-2+iFrame),:);
            nCand = size(pos_frame_now,1);
            track_slice(iFrame).xCoord = [pos_frame_now(:,3),zeros(nCand,1)]; %careful, check if error should be >0!
            track_slice(iFrame).yCoord = [pos_frame_now(:,4),zeros(nCand,1)];
%             track_slice(iFrame).amp = [pos_frame_now(:,5),zeros(nCand,1)];
            
            track_slice(iFrame).zCoord = [pos_frame_now(:,5),zeros(nCand,1)];
            track_slice(iFrame).amp = [pos_frame_now(:,6),zeros(nCand,1)];
            if iFrame==1
                pos_frame_first = pos_frame_now; %need this for getting trajectory ids back
            end
        end
        
        %track those two frames
        if ~isempty(pos_frame_now) && ~isempty(pos_frame_first)
            try
                [tracksFinal_sliced,~,~] = trackCloseGapsKalmanSparse(track_slice,costMatrices_slice,gapCloseParam_slice,kalmanFunctions_slice,probDim,0,false);
                nrTracks_sliced = numel(tracksFinal_sliced);
            catch
                %very rarely, utrack tries to index a feature track
                %matrix with a wrong relative index starting at 0
                nrTracks_sliced = 0;
            end
        else
            nrTracks_sliced = 0;
        end
        
        if nrTracks_sliced>0
            for iTrack = 1:nrTracks_sliced
                x_pos = tracksFinal_sliced(iTrack).tracksCoordAmpCG(1:8:end).'; %[x_old,x_new]
                
                %we cannot get back the respective ids, so we have
                %to search for positions instead (risky - use fuzzy
                %equal)
                id_pair = [pos_frame_first(abs(pos_frame_first(:,3)-x_pos(1))<1e-6,1);pos_frame_now(abs(pos_frame_now(:,3)-x_pos(2))<1e-6,1)];
                
                %now correct ids
                connected_track_idx = (trackingData(:,1)==id_pair(2)); %get all boolean entries of new trajectory to be connected to older one
                trackingData(connected_track_idx,1) = repmat(id_pair(1),sum(connected_track_idx),1);
            end
            
            %some of the ids were reverted to the smaller values of
            %trajectories in old slice. all unconnected or new
            %tracks will have an id which is too high. we'll
            %correct this now
            traj_update = 0;
            for jId = 1+traj_id:nrTracks+traj_id %go through all tracks in newly added slice
                idx_trackingData_update = trackingData(:,1)==jId;
                n_to_update = sum(idx_trackingData_update); %is this an unconnected or new track?
                if n_to_update>0
                    trackingData(idx_trackingData_update,1) = repmat(traj_id+traj_update+1,n_to_update,1); %then set back the id...
                    traj_update = traj_update+1; %increment the updater id...
                end
            end
            traj_id = traj_id+traj_update; %...and finally update the global id
        else
            traj_id = traj_id+nrTracks; %we are not in the first slice and there's nothing to reconnect. it's time to update the global id
        end %if nrTracks_sliced>0
    else
        traj_id = traj_id+nrTracks; %if in very first slice, there's nothing to reconnect. then just update the global id
    end %if iDiv>1
end %for iDiv=1:nrSplit

if ~isempty(trackingData)
    %remove NaNs resulting from closed gaps
    trackingData = trackingData(~isnan(trackingData(:,end)),:);
    
    % if the movie was split, result array has to be sorted by
    % trajectory id again
    if nrSplit>1
        [~,idx] = sort(trackingData(:,1));
        trackingData = trackingData(idx,:);
    end
end

end



function [gapCloseParam,costMatrices,kalmanFunctions] = parseUtrackOptions(trackingOptions)
%taken from utrack: scriptTrackGeneral.m

%%%%%DO NOT EVER CHANGE THIS VALUE
gapCloseParam.timeWindow = trackingOptions.maxFrameGap+1; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
%%%%%

gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.

%%%%DO NOT EVER CHANGE THIS VALUE
gapCloseParam.minTrackLen = trackingOptions.minTrajLength; %minimum length of track segments from linking to be used in gap closing.
%%%%%


%optional input:
gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

%parameters
parameters.linearMotion = 0; %use linear motion Kalman filter.
parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.

%%%%%DO NOT EVER CHANGE THIS VALUE
parameters.maxSearchRadius = trackingOptions.maxTrackRadius; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
%%%%%


parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
% parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.

%optional input
parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

%parameters

%needed all the time
parameters.linearMotion = costMatrices(1).parameters.linearMotion; %use linear motion Kalman filter.

parameters.minSearchRadius = costMatrices(1).parameters.minSearchRadius; %minimum allowed search radius.
parameters.maxSearchRadius = costMatrices(1).parameters.maxSearchRadius; %maximum allowed search radius.
parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

parameters.brownScaling = [0 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
% parameters.timeReachConfB = 3; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
parameters.timeReachConfB = gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

parameters.ampRatioLimit = [0.7 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.linStdMult = 1*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

parameters.linScaling = [0.25 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
% parameters.timeReachConfL = 4; %similar to timeReachConfB, but for the linear part of the motion.
parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).

%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';


end


