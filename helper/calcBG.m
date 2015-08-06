function [bg_mean,bg_std] = calcBG(im,ywidth,xwidth,halfw)
%adapted from spatialMovAveBG.m developed as part of the "utrack" particle tracking software suite

window = 2*halfw+1;
px_start = window+halfw;
px_end_y = max(ywidth - (px_start-1),px_start);
px_end_x = max(xwidth - (px_start-1),px_start);

bg_mean(ywidth,xwidth) = 0;
bg_std = bg_mean;

% go through every block
for i_x = px_start : window : px_end_x
    for i_y = px_start : window : px_end_y
        
        im_tmp = im(i_y-(px_start-1):min(i_y+px_start-1,ywidth),i_x-(px_start-1):min(i_x+px_start-1,xwidth)); %get tmp image
        
        [bg_mean_tmp,bg_std_tmp] = lms_fun(im_tmp(:)); %use LMS filter to extract outliers (maxima) from the image and then calculate mean and std of the remaining pixels which constitute the background
        bg_std_tmp = max(bg_std_tmp,eps);
        
        bg_mean(i_y-halfw:i_y+halfw,i_x-halfw:i_x+halfw) = bg_mean_tmp;
        bg_std(i_y-halfw:i_y+halfw,i_x-halfw:i_x+halfw) = bg_std_tmp;
    end
end

% Fill empty space with border values
px_first = px_start - halfw;
last_y = i_y + halfw;
last_x = i_x + halfw;

for i_x = px_first : last_x
    bg_mean(1:px_first-1,i_x) = bg_mean(px_first,i_x);
    bg_mean(last_y+1:end,i_x) = bg_mean(last_y,i_x);
    bg_std(1:px_first-1,i_x) = bg_std(px_first,i_x);
    bg_std(last_y+1:end,i_x) = bg_std(last_y,i_x);
end
for i_y = 1 : ywidth
    bg_mean(i_y,1:px_first-1) = bg_mean(i_y,px_first);
    bg_mean(i_y,last_x+1:end) = bg_mean(i_y,last_x);
    bg_std(i_y,1:px_first-1) = bg_std(i_y,px_first);
    bg_std(i_y,last_x+1:end) = bg_std(i_y,last_x);
end

% % upper/lower
% bg_mean([1:window-1;px_end_y+halfw+1:ywidth],window:px_end_x+halfw) = bg_mean([window:2*window-2;px_end_y-halfw+1:px_end_y+halfw],window:px_end_x+halfw);
% bg_std([1:window-1;px_end_y+halfw+1:ywidth],window:px_end_x+halfw) = bg_std([window:2*window-2;px_end_y-halfw+1:px_end_y+halfw],window:px_end_x+halfw);
% 
% % left/right
% bg_mean(window:px_end_y+halfw,[1:window-1,px_end_x+halfw+1:xwidth]) = bg_mean(window:px_end_y+halfw,[window:2*window-2,px_end_x-halfw+1:px_end_x+halfw]);
% bg_std(window:px_end_y+halfw,[1:window-1,px_end_x+halfw+1:xwidth]) = bg_std(window:px_end_y+halfw,[window:2*window-2,px_end_x-halfw+1:px_end_x+halfw]);
% 
% % edges
% bg_mean([1:window-1;px_end_y+halfw+1:ywidth],[1:window-1,px_end_x+halfw+1:xwidth]) = bg_mean([window:2*window-2;px_end_y-halfw+1:px_end_y+halfw],[window:2*window-2,px_end_x-halfw+1:px_end_x+halfw]);
% bg_std([1:window-1;px_end_y+halfw+1:ywidth],[1:window-1,px_end_x+halfw+1:xwidth]) = bg_std([window:2*window-2;px_end_y-halfw+1:px_end_y+halfw],[window:2*window-2,px_end_x-halfw+1:px_end_x+halfw]);

end %function


function [data_mean, data_std] = lms_fun(data)
% LMS calculates mean and standard deviation of a vector while discarding outliers.
% The algorithm is based on "Least Median Squares" (Danuser 1992, Rousseeuw & Leroy 1987) and this function is inspired by "robustMean.m" supplied by utrack.
% utrack dan be downloaded from <http://lccb.hms.harvard.edu/software.html>, see also "Jaqaman et al., Nature Methods 5, pp. 695-702 (2008)"

n_data = numel(data);

% Define magic numbers:
k=3; %cut-off is roughly at 3 sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987
magic_number2=1.4826^2; %see same publications

% Get data statistics
data_median = median(data);
res2 = (data-data_median(ones(n_data,1))).^2;
medRes2 = max(median(res2),eps);

% Get test statistics
test_div = 1/(medRes2*magic_number2);
test_vec = test_div*res2;

% Accept values within 3*sigma of test statistics as inliers
inlier_bool = test_vec<=k^2;
n_inlier = sum(inlier_bool);

% Calculate mean and standard deviation
data_mean = mean(data(inlier_bool));
% data_std=sqrt(sum(res2(inlier_bool))/(n_inlier-4));
data_std = std(data(inlier_bool));

end %function
