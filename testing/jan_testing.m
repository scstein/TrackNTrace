match_data = cell(10,1);
match_img = match_data;
CORR_THRESH = 0.2;
img_bck = movie(64:end,64:end,113); img_bck = img_bck(:);
bck_mean = mean(img_bck);
bck_std = std(img_bck);
sigma_backup = 520/(2*1.49)/(2*sqrt(2*log(2)))/80;
pvalue = 0.05;
perc_cut = 10;


imax = 5;

% frame 113
for j=1:20
img = movie(:,:,1+10*j);

image_corr = zeros(128,128);
image_corr_bool = image_corr;

for i=1:imax
    sigma = sigma_backup*(1+(i-1)*0.1);
    model = gaussian2d(round(10),sigma,1);
    pattRows = size(model,1);
    pattCols = size(model,2);
    
    MIN_DIST = round(5*sigma);
    [ match_data_tmp, match_img_tmp ] = matchSpot( img, model, CORR_THRESH, MIN_DIST);
%     halfw = (round(3*sigma)-1)/2;
%     [img_bck,img_std] = calcBG(img,128,128,halfw);

    
    match_data_tmp(:,1) = match_data_tmp(:,1)+ceil(pattRows/2);
    match_data_tmp(:,2) = match_data_tmp(:,2)+ceil(pattCols/2);
    
    idx = sub2ind(size(img),match_data_tmp(:,1),match_data_tmp(:,2));
    match_data_tmp = [match_data_tmp,img(idx),idx];
    
    image_corr = image_corr + match_img_tmp; %correlation image 
    image_corr_bool(idx) = image_corr_bool(idx)+1;
    pValue = 1 - normcdf(match_data_tmp(:,4),bck_mean,bck_std);
    
    match_data(i) = {match_data_tmp};
    match_img(i) = {match_img_tmp};
%     figure;
%     mim(img);
%     hold on
%     plot(match_data_tmp(:,2)+0.5,match_data_tmp(:,1)+0.5,'go');
%     plot(match_data_tmp(pValue>0.05,2)+0.5,match_data_tmp(pValue>0.05,1)+0.5,'bo');
%     hold off
%     pause
end


%create images
image_pval = reshape(1 - normcdf(img(:),bck_mean,bck_std),128,128);
image_pval(image_pval>pvalue | ~image_corr_bool) = NaN;

img_sort = sort(img(:));
intens_cut = img_sort(round((1-perc_cut/100)*size(img_sort,1)));
image_rval = img;
image_rval(img<intens_cut | ~image_corr_bool) = NaN;

%normalize
image_corr = image_corr/imax;
image_pval = (0.05-image_pval)/0.05;
image_rval = (image_rval-min(image_rval(:)))/(max(image_rval(:)-min(image_rval(:))));

score = [0.01,0.01,0.95,0.15];
mov_avg = ones(3,3)/9;
image_corr = conv2(image_corr,mov_avg,'same');
image_corr_bool = conv2(image_corr_bool,mov_avg,'same');



idx_total = find((image_corr(:)>=score(1)) &  (image_pval(:)>=score(3)) & (image_rval(:)>=score(4)));

[pos_y,pos_x] = ind2sub(size(img),idx_total);

% figure; imagesc(image_corr);
% figure; imagesc(image_corr_bool/imax);
% figure; imagesc(image_pval);
% figure; imagesc(image_rval);


figure;
mim(img);
hold on
plot(pos_x+0.5,pos_y+0.5,'go');
hold off

end

% idea:
% - determine score: correlation score, pval score, rval score and normalize to [0,1] each
% - user gives score thresholds: individual score, total score has to exceed threshold
% - select high score candidates
% - check direct environment for other high score spots and fit if possible

