function [ drift ] = calc_drift_from_tracks( filepathOrData, minLength, maxLength, minParticles )
% USAGE: calc_drift_from_tracks( filepathOrData, minLength, maxLength, minParticles )
%This calculates the drift from the particle tracks. The
%drift is simply calculated as the average of particle displacements in
%each frame.
%
% Note: If a path to a TNT file is given as input, the result is
% automatically appended to that TNT file.
%   Input:
%         filepathOrData -  path to TNT file 
%                               OR
%                           trackingData: 2D double array, list of trajectories with
%                           columns [id,frame,xpos,ypos, ... ]. 
%                           (The first possible frame must be 1 not 0)
%         minLength - Minimum length of trajectories to use for computation | default = 1         
%         maxLength - Maximum length of trajectories to use for computation | default = inf
%         minParticles - Minimum number of particles (/displacements) which | default = 1
%                        must be present in a frame to be used 
%    Note: Assumes no drift when no fitting trajectories are present 
%          (or less then the user specified).
%
%   Output:
%      drift - The computed drift with columns [frame, x,y].
%
% Author: Simon Christoph Stein

% Parse inputs
if ischar(filepathOrData)
   load(filepathOrData,'trackingData','-mat');
   if ~exist('trackingData','var')
       error('No trackingData found in specified file');
   end
else
   trackingData = filepathOrData; 
end

if nargin<2 || isempty(minLength)
   minLength = 1; 
end
if nargin<3 || isempty(maxLength)
   maxLength = inf; 
end
if nargin<4 || isempty(minParticles)
   minParticles = 1; 
end

%% Get number of frames
min_frame = min(trackingData(:,2));
max_frame = max(trackingData(:,2));
n_frames =  max_frame-min_frame+1;


%% Convert data to cell, drop particle identifier
id_tracks = unique(trackingData(:,1));

tracks = {};
cnt = 1;
for iTrack = id_tracks.'
    trackData = trackingData( trackingData(:,1)== iTrack, 2:4);
    if (size(trackData,1)< minLength) || (size(trackData,1) > maxLength)
        continue
    end
    tracks{cnt} = trackData;
    cnt = cnt+1;
end

%%
% NOTE@Simon: What follows is almost 100% the calc_drift_from_tracks function 
% once used in the software I wrote for the cryo

% Variables
n_tracks = numel(tracks);

drift = zeros(n_frames-1,2);
nr_non_nans = zeros(n_frames-1,1);

% Here we sum up the frame-to-frame displacement of all particles present 
% in both frame i and i+1. If a particle is not present in one of these
% frames, its position should be [NaN NaN] and it is not counted.
for i_track = 1:n_tracks-1
    
    % To ease further computation, we fill every time intervall not present
    % for a particle with NaNs
    filled_track = NaN(n_frames,2);
    filled_track(tracks{i_track}(:,1)-min_frame+1,:) = tracks{i_track}(:,2:3);
    
    displacement = diff(filled_track);
    non_nans = ~isnan( displacement(:,1) );
    
    % Set NaNs to 0 in the displacement vector (no contribution from non-
    % present particles.
    displacement(~non_nans,:) = 0;
    
    drift = drift + displacement;
    nr_non_nans = nr_non_nans + non_nans;
end    
    % Average over the number of contributing particles (non-nans)
    drift = drift ./ nr_non_nans(:,ones(1,size(drift,2))); 
    missing_data_frames = isnan(drift(:,1)) | isnan(drift(:,2)); % frames where no data is available
    nr_bad_frames = sum(missing_data_frames);
    
    if( nr_bad_frames > 0)       
       fprintf('  Missing drift data for %i frames! Assuming zero drift during this time.\n', nr_bad_frames);
    end    
    drift(missing_data_frames, :) = zeros(nr_bad_frames,2);
        
    % Set drift to zero if not enough particles were present during the computation
    untrustworthy_frames = logical(nr_non_nans> 0) & logical(nr_non_nans <  minParticles);
    nr_untrustworthy_frames = sum(untrustworthy_frames);
    if( nr_untrustworthy_frames > 0)
       fprintf('  Not enough particles for %i frames! Assuming zero drift during this time.\n', nr_untrustworthy_frames);
    end    
    drift(untrustworthy_frames, :) = zeros(nr_untrustworthy_frames,2);
    
    fprintf('  Successful drift calculation for %i frames.\n', n_frames-nr_bad_frames-nr_untrustworthy_frames);
    
    % Output some helpful numbers
    used_nr_particles = nr_non_nans(find(nr_non_nans>=minParticles)); % Filter out the entries used for computation
    fprintf('   -> Average number of tracks per frame used for estimation: %.0f. (Min: %i, Max: %i)\n', mean(used_nr_particles), min(used_nr_particles), max(used_nr_particles));
    
    
    % The first frame is the reference, such there is no drift.
    drift = [0 0; 
            drift];
    
    % Up to now, drift contains the frame-to-frame displacement.
    % To get the position with respect to the first frame, we have 
    % to sum up the displacements of all preceding frames.
    drift = cumsum(drift,1);
    
    % Add the frame number to the calculated drift curve
    drift = horzcat( (min_frame:max_frame)', drift);
    
    
    % Save result to TNT file if given
    if ischar(filepathOrData)
        driftOptions.minLength = minLength;
        driftOptions.maxLength = maxLength;
        driftOptions.minParticles = minParticles; %#ok<STRNU>
        
        save(filepathOrData, 'driftOptions','drift','-append');
    end
end


% end

