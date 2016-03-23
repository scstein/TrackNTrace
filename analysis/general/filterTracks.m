function [ filterData ] = filterTracks( trackingDa, minLength, maxLength )
% Filter out tracks between minLength and maxLength

if nargin<2 || isempty(minLength)
    minLength = 0;
end
if nargin <3 || isempty(maxLength)
    maxLength = inf;
end

trackIDs = unique(trackingDa(:,1));

filterData = {};
cnt = 1;
for iTrack = trackIDs.'
    trackLength = numel( find(trackingDa(:,1)==iTrack) );
    if(trackLength >= minLength && trackLength<=maxLength)
        filterData{cnt} = trackingDa(trackingDa(:,1)==iTrack,:);
        cnt = cnt +1;
    end
end

filterData = vertcat(filterData{:});

end



