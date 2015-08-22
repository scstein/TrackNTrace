function fitData = filterPositionsSNB(fitData,SNB_THRSH)
% [fitDataFiltered] = filterPositionsSNB(fitData,SNB) This function is used
% to filter fitted positions obtained from TrackNTrace by their respective
% signal-to-background ratio defined as A/B where A is the amplitude and B
% is the background.
% 
% INPUT:
%     fitData: 3D double array of fitted positions obtained in the
%     locateParticles step of TrackNTrace
%     
%     SNR: double, signal-to-background ratio
%     
% OUTPUT:
%     fitDataFiltered: filtered array, all positions A/B<=SNB have
%     their error flag set to -3 resulting in their removal during the
%     tracking step



for iFrame = 1:size(fitData,3)
    fitData_frame = fitData(:,:,iFrame);
    pos_not_kept = fitData_frame(3,:)./fitData_frame(4,:) <= SNB_THRSH;
    
    fitData(end,pos_not_kept,iFrame) = repmat(-3,1,sum(pos_not_kept));
end
    
    
    