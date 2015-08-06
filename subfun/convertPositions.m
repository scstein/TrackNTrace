function [pos] = convertPositions(pos_file,method,probDim)
% [posArray] = convertPositions(posFilename,method)
% Reads file created by locateParticles.m and converts the array to a data format suitable for the respective tracker
%
% INPUT:
%     posFilename: string, Full path to .mat file created with locateParticles.m
%
%     method: string, name of tracker
%
% OUTPUT:
%     posArray: N-D array suited to a specific tracker


%main input is fitData array containing 6 lines for x,y,amp,bckground,sigma,exitflag, N_max columns with (N_max maximum overall number of candidates in any frame) and nrFrames indices in 3D
% if pos_file is not a path but an array, it's a fitData array and we're in
% test mode, so we don't need to read anything
if ischar(pos_file)
    load(pos_file);
else
    fitData = pos_file;
end
nrFrames = size(fitData,3);


switch(method)
    case('simpletracker') %NNT, simpletracker by Jean-Yves Tinevez, released 2011 on Matlab file exchange
        %init cell array
        pos = cell(nrFrames,1);
        for i=1:nrFrames
            pos_frame_now = squeeze(fitData(1:5,:,i)).'; %get x,y,amp,backround,sigma
            pos_frame_now(pos_frame_now==0) = 1e-6; %this is a dirty hack for particles which run out of the frame
            valid_pos = any(pos_frame_now,2); %get valid positions, zero positions count!
            pos(i) = {pos_frame_now(valid_pos,1:probDim)}; %1D cell array with nrFrames lines, each cell containing 2D array [x,y] or bigger array [x,y,amp], [x,y,amp,background], ...
        end
        
    case('utrack') %complex tracker searching for global optimum, VERY memory intensive for large number of particles
        pos = repmat(struct('xCoord',[],'yCoord',[],'amp',[],'sigma',[]),nrFrames,1); %1D struct array with nrFrames lines, inner arrays have two columns [value,error]
        for i=1:nrFrames
            pos_frame_now = squeeze(fitData(1:5,:,i)).';
            valid_pos = any(pos_frame_now,2);
            nCand = sum(valid_pos);
            pos_frame_now(pos_frame_now==0) = 1e-6; %this is a dirty hack for particles which run out of the frame
            if nCand>0
                pos(i).xCoord = [pos_frame_now(valid_pos,1),zeros(nCand,1)]; %careful, check if error should be >0!
                pos(i).yCoord = [pos_frame_now(valid_pos,2),zeros(nCand,1)];
                pos(i).amp = [pos_frame_now(valid_pos,3),zeros(nCand,1)];
                pos(i).sigma = [pos_frame_now(valid_pos,5),zeros(nCand,1)];
            end
        end
        
    case('track_cg') %NNT, Crocker and Grier, 1994
        nrPos = zeros(nrFrames,1);
        for i=1:nrFrames
            nrPos(i) = sum(any(squeeze(fitData(3,:,i)),1));
        end
        pos = zeros(sum(nrPos),4);
        
        idx_start = 1;
        for i=1:nrFrames
            idx_end = sum(nrPos(1:i));
            
            %2D double array [x,y,amp,frame] with data from several frames stacked vertically. probDim is ignored here and dealt with later
            pos(idx_start:idx_end,:) = [fitData(1:3,1:nrPos(i),i).',repmat(i,nrPos(i),1)];
            idx_start = idx_end+1;
        end
        pos(pos==0) = 1e-6; %this is a dirty hack for particles which run out of the frame
end
