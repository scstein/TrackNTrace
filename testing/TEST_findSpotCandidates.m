%%
movie = read_tiff('U:\Weixing\weixing\101 polar origami\bleaching-3.7exci-roi2.tif');
avg = mean(movie(:,:,end-200:end),3);

%% Candidate selection by normalized cross correlation
model = gaussian2d(10,1);
pattRows = size(model,1);
pattCols = size(model,2);

CORR_THRESH = 0.05;
MIN_DIST = 10;
[ match_data, match_img ] = matchSpot( avg, model, CORR_THRESH, MIN_DIST);
match_data(match_data(:,3)<0.8,:) = [];

[~,idx] = sort(match_data(:,3),'descend');
match_data = match_data(idx,:);

%% Post processing selection
amplitudes = avg(logical(match_img));
[amplitudes, idx] = sort(amplitudes);
fitPos = (round(size(amplitudes)/2)-round(size(amplitudes)*0.25) : round(size(amplitudes)/2)+round(size(amplitudes)*0.25))';
ampli_fit = amplitudes(fitPos);

p = polyfit(fitPos,ampli_fit,1);
yPoly = polyval(p,(1:numel(amplitudes))');
margin = 0.025;

residuals = abs(amplitudes-yPoly);
inliers = residuals < margin*yPoly;
xInl = find(inliers);
ampInl = amplitudes(inliers);

figure;
plot( 1:numel(amplitudes),amplitudes,'r.', 1:numel(amplitudes), yPoly,'b', 1:numel(amplitudes), yPoly*(1+margin),'b--',1:numel(amplitudes), yPoly*(1-margin),'b--', xInl,ampInl,'g');

match_data = match_data(inliers,:);

%% Plot result
figure;
mim(avg);
hold on
plot(match_data(:,2)+ceil(pattCols/2),match_data(:,1)+ceil(pattRows/2),'co');
hold off; 

%% Plot result (to test match img)
cand_img = avg;
cand_img(match_img>0) = 1.05*max(avg(:));
figure; mim(cand_img);

%% Plot candidates individually with local contrast
half_wsize = 3;
figure;
for i=1:size(match_data,1)
    
    row_pos = match_data(i,1)+ceil(pattRows/2);
    col_pos = match_data(i,2)+ceil(pattCols/2);
    
   imagesc(avg); axis image; colormap hot;
   hold on
    plot(match_data(:,2)+ceil(pattCols/2),match_data(:,1)+ceil(pattRows/2),'co');
    
   title(sprintf('%i/%i\n', i,size(match_data,1)));
   % scale local contrast
   cutout = avg(max(1,row_pos-4*half_wsize):min(size(avg,1),row_pos+4*half_wsize), max(1,col_pos-4*half_wsize):min(size(avg,2),col_pos+4*half_wsize));
   caxis([min(cutout(:)) max(cutout(:))])
   ylim([row_pos-4*half_wsize row_pos+4*half_wsize]);
   xlim([col_pos-4*half_wsize col_pos+4*half_wsize]);
   
   pause;
end



