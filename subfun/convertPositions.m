function [pos] = convertPositions(pos_file,method)
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


%main input is fitData cell array containing 1 cell per frame with 6 lines for x,y,amp,bckground,sigma,exitflag, N_max columns with (N_max maximum overall number of candidates in any frame)
% if pos_file is not a path but an array, it's a fitData array
if ischar(pos_file)
    load(pos_file, 'fitData', '-mat'); % '-mat' forces loading as MAT-file
else
    fitData = pos_file;
end
nrFrames = size(fitData,1);


switch(method)        
    case('utrack') %complex tracker searching for global optimum, VERY memory intensive for large number of particles
        pos = repmat(struct('xCoord',[],'yCoord',[],'amp',[],'sigma',[]),nrFrames,1); %1D struct array with nrFrames lines, inner arrays have two columns [value,error]
        for iFrame=1:nrFrames
            pos_frame_now = fitData{iFrame}.';
            valid_pos = pos_frame_now(:,6)==1; %error flag is 1?
            nCand = sum(valid_pos);
            pos_frame_now(pos_frame_now==0) = 1e-6; %this is a dirty hack for particles which run out of the frame
            if nCand>0
                pos(iFrame).xCoord = [pos_frame_now(valid_pos,1),zeros(nCand,1)]; %careful, check if error should be >0!
                pos(iFrame).yCoord = [pos_frame_now(valid_pos,2),zeros(nCand,1)];
                pos(iFrame).amp = [pos_frame_now(valid_pos,3),zeros(nCand,1)];
                pos(iFrame).sigma = [pos_frame_now(valid_pos,5),zeros(nCand,1)];
            end
        end
        
    case('nn_cpp')
        nrPos = zeros(nrFrames,1);
        for i=1:nrFrames
            nrPos(i) = sum((fitData{i}(6,:)==1)); %error flag is 1?
        end
        pos = zeros(6,sum(nrPos));
        
        nrPos_cs = [0;cumsum(nrPos)];
        
        for jFrame = 1:nrFrames           
            pos_frame_now = fitData{jFrame};
            valid_pos = pos_frame_now(6,:)==1; %error flag is 1?
            pos_frame_now(pos_frame_now==0) = 1e-6; %this is a dirty hack for particles which run out of the frame
            pos(:,nrPos_cs(jFrame)+1:nrPos_cs(jFrame+1)) = [repmat(jFrame,1,nrPos(jFrame));pos_frame_now(1:5,valid_pos)];
        end
        
    otherwise 
        error('Unknown tracker ''%s''',method);
end
