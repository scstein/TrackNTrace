%% exclude emitters that move to fast

fitData_cp = fitData;

%% mark invalid based on sigma
fitData(:,squeeze(fitData(5,:,:))<0.8) = NaN;
fitData(:,squeeze(fitData(5,:,:))>1.2) = NaN;

%% mark invalid based on amplitude
fitData(:,squeeze(fitData(3,:,:))<50) = NaN;

%% Exclude emitters with short traces
min_valid_fits = 1065;
NonNaN_check = ~isnan(squeeze(sum(fitData,1))); % was any parameter set to NaN?
exclude_mask = sum(NonNaN_check,2)<min_valid_fits;
fitData(:,exclude_mask,:) = NaN;

%% - Exclude emitters based on displacement - (works well but very restrictive)
% fitData = fitData_cp;

xpos = squeeze(fitData(1,:,:));
ypos = squeeze(fitData(2,:,:));

% Integrate displacement over 5 frames
displace = sqrt(diff(xpos,1,2).^2 + diff(ypos,1,2).^2);
%   displace = (displace(:,1:2:end) + displace(:,2:2:end))/2;
delete_mask = max(displace,[],2)>1;
fitData(:,delete_mask,:) = [];


%% Plot movie and fit results
figure;
for iF = 1:size(movie,3)
    if(iF ==1)
        xl = [0.5,size(movie,2)+0.5];
        yl = [0.5,size(movie,1)+0.5];
        frame = movie(:,:,1);
        zl = [min(frame(:)), max(frame(:))];
    else
       xl = xlim;   
       yl = ylim;
    end
   
   
   mim(movie(:,:,iF)); title(sprintf('Frame %i/%i',iF,size(movie,3)));
   hold on;
   plot( squeeze(fitData(1,:,iF)),squeeze(fitData(2,:,iF)),'go');
   hold off;
   
   if(iF == 1)
       imcontrast;
       pause;
       zl = caxis;
   end
   
   
   xlim(xl);
   ylim(yl);
   caxis(zl);
   
   pause(0.05)
end

%% Show particle traces
figure; plot(squeeze(fitData(1,:,:))', squeeze(fitData(2,:,:))','.-'); 

