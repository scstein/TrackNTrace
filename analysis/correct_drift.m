function [ driftCorr_trajectoryData ] = correct_drift( filepathOrData, drift)
% Subtracts the drift present in the given trajectories using the supplied
% drift curve and returns the corrected data
%
% Note: If a path to a TNT file is given as input, the result is
% automatically appended to that TNT file.
%   Input:
%         filepathOrData -  path to TNT file 
%                               OR
%                           trajectoryData: 2D double array, list of trajectories with
%                           columns [id,frame,xpos,ypos, ... ]. 
%                           (The first possible frame must be 1 not 0)
%   Output:
%         driftCorr_trajectoryData -  2D double array, corrected positions list of 
%                          trajectories with columns [id,frame,xpos,ypos, ... ]. 
%
% Author: Simon Christoph Stein

% Parse input
if ischar(filepathOrData)
    load(filepathOrData,'trajectoryData','-mat');
    if ~exist('trajectoryData','var')
        error('No trajectoryData found in specified file.');
    end
    
    if nargin<2 || isempty(drift)
        load(filepathOrData,'drift','-mat');
        if ~exist('drift','var')
            error('No drift data found in specified file. Compute drift first.');
        end
    end
else
    trajectoryData = filepathOrData;
end

%% Subtract the drift
 min_frame = min(trajectoryData(:,2));
 min_frame_drift = min(drift(:,1));
 
 % Sanity check
 if(min_frame ~= min_frame_drift)
    error('drift data and trajectory data don''t seem to match!') 
 end
 
  driftCorr_trajectoryData = trajectoryData; 
 % We use indexing to subtract the right drift component from each
 % localization. To do this we offset the frameID by subtracting
 % min_frame+1, such that data taken at min_frame is mapped to the first
 % row of the drift vector.
 driftCorr_trajectoryData(:,3:4) = driftCorr_trajectoryData(:,3:4)-drift(trajectoryData(:,2)-min_frame+1,2:3);
 
 
 % Save result to TNT file if given
if ischar(filepathOrData)
    save(filepathOrData, 'driftCorr_trajectoryData','-append');
end
    
end