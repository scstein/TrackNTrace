function trackLengthHist(trackingData)    
% Compute a histogram of the length of all trajectory

id_tracks = unique(trackingData(:,1));    
n_tracks = numel(id_tracks);

traj_length = size(n_tracks,1);
cnt = 1;
for iTrack = 1:n_tracks
    traj_length(cnt) = size(trackingData( trackingData(:,1)== id_tracks(cnt) , 2:4),1);
    cnt = cnt+1;
end

xcenters = min(traj_length(:)):max(traj_length(:));

figure;
hist(traj_length,xcenters);
xlabel('Track length');
ylabel('frequency')

end