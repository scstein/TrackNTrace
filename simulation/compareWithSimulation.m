function [resultStatistics] = compareWithSimulation(fitData, fitData_truth, options)


if numel(fitData)~=numel(fitData_truth)
    error('Outcome and ground truth have to have identical number of frames. Aborting...\n');
end

nrFrames = numel(fitData);
minTrackLength = 2;
linkingDistance = options.linkingDistance;

nrTruePositive = 0;
nrFalsePositive = zeros(nrFrames,1);
nrFalseNegative = zeros(nrFrames,1);
RMSE = zeros(nrFrames,1);
emptyFrames = [];


for iFrame = 1:nrFrames
    if ~isempty(fitData{iFrame}) && ~isempty(fitData_truth{iFrame})
        pos = preparePosArray(fitData{iFrame},fitData_truth{iFrame});
        track = mx_nn_tracker(pos,minTrackLength,linkingDistance,[],[],minTrackLength,false).';
        if isempty(track)
            nrFalsePositive(iFrame) = size(fitData{iFrame},1);
            nrFalseNegative(iFrame) = size(fitData_truth{iFrame},1);
            emptyFrames = [emptyFrames;iFrame];
            continue
        end
        nrIds = max(track(:,1));
        
        nrTruePositive = nrTruePositive+nrIds;
        
        nrFalsePositive(iFrame) = size(fitData{iFrame},1)-nrIds;
        nrFalseNegative(iFrame) = size(fitData_truth{iFrame},1)-nrIds;

        RMSE(iFrame) = 1/max(track(:,1))*sum(sum((track(2:2:end,3:4)-track(1:2:end-1,3:4)).^2,2),1); %1/N*sum( (x-x_true)^2+(y-ytrue)^2)
    else
        nrFalsePositive(iFrame) = size(fitData{iFrame},1);
        nrFalseNegative(iFrame) = size(fitData_truth{iFrame},1);
        emptyFrames = [emptyFrames;iFrame];
    end
end
RMSE(emptyFrames) = [];


resultStatistics.accuracy = [sqrt(mean(RMSE)), 0.5*std(RMSE)*1/sqrt(mean(RMSE))]*options.pixelSize; %value, standard deviation
resultStatistics.truePositive = nrTruePositive;
resultStatistics.falsePositive = [sum(nrFalsePositive),std(nrFalsePositive)*nrFrames];
resultStatistics.falseNegative = [sum(nrFalseNegative),std(nrFalseNegative)*nrFrames];
JAC = nrTruePositive/(nrTruePositive+sum(nrFalsePositive)+sum(nrFalseNegative));
JAC = [JAC,JAC/(nrTruePositive+sum(nrFalsePositive)+sum(nrFalseNegative))*sqrt(resultStatistics.falsePositive(2)^2+resultStatistics.falseNegative(2)^2)];
resultStatistics.JAC = JAC;

resolutionFRC = [];
if options.calcFRC
    randomize = false;
    [q,frc] = calcFRCCurve(fitData, fitData_truth,options,randomize);
    significantP = find(frc<0.1,1,'first')/numel(frc); %smooth relevant part of curve...
    frcSmooth = smooth(frc,significantP/10,'rloess'); %in about 10 windows, robust lsq regression
    resolutionFRC = 1/q(find(frcSmooth<1/7,1,'first'));
    FRC_error = [1/q(find(frc<1/7,1,'first')), 1/q(find(frc<1/7,1,'last'))];
    FRC_error = abs(diff(FRC_error));
    
end
resultStatistics.resolutionFRC = [resolutionFRC,FRC_error/2];

end


function [pos,id_true,id_loc] = preparePosArray(frame1,frame2)

%give pos array for tracker, first frame is ground truth

pos_loc1 = [2*ones(1,size(frame1,1)); frame1(:,1:2).'; 1:size(frame1,1)];
pos_loc2 = [ones(1,size(frame2,1)); frame2(:,1:2).'; -(1:size(frame2,1))];
pos = [pos_loc2,pos_loc1]; %track from ground truth to localized data. in general, localization will be erroneous and contain false negatives
id_true = pos_loc2(end,:).';
id_loc = pos_loc1(end,:).';

pos = [pos(1:3,:);zeros(1,size(pos,2));pos(end,:)]; %insert z=0 as nntracker tracks in third dimension by default

end



function [q,frc_curve] = calcFRCCurve(fitData1, fitData2, options, randomize)
persistent idx_radial_search

px_size = options.pixelSize;
p = options.zoomFactor;
im_w = options.imageWidth; %we can only really use one width, so user can put in the maximum one or something


% get localization tables
if randomize || isempty(fitData2)
    % draw random samples from fitData
%     fitData = [vertcat(fitData1{:});vertcat(fitData2{:})];
    fitData = vertcat(fitData1{:});  
    fitData = fitData(:,1:2);
%     rand_bool = logical(randi([0,1],1,size(fitData,1))).';
    rand_bool = randperm(size(fitData,1),size(fitData,1));
    rand_bool = rem(rand_bool/2,1)==0;
    loc1 = fitData(rand_bool,:);
    loc2 = fitData(~rand_bool,:);
else
    loc1 = vertcat(fitData1{:}); loc1 = loc1(:,1:2);
    loc2 = vertcat(fitData2{:}); loc2 = loc2(:,1:2);
    
    loc1 = loc1(loc1(:,1)>=3*256/8&loc1(:,1)<5*256/8&loc1(:,2)>=3*256/8&loc1(:,2)<5*256/8,:)-3*256/8;
    loc2 = loc2(loc2(:,1)>=3*256/8&loc2(:,1)<5*256/8&loc2(:,2)>=3*256/8&loc2(:,2)<5*256/8,:)-3*256/8;
end


% bin localizations
ctrs = linspace(0+0.5,im_w-1+0.5,im_w*p).';
img1 = hist3(loc1,{ctrs,ctrs}).';
img2 = hist3(loc2,{ctrs,ctrs}).';


% apply 2D tukey window at 1/8, 7/8 edges and transform
mask = tukeywin(im_w*p,2*1/8);
mask = mask*mask.';
% DC component for even width is at [width/2+1,width/2+1]
img1 = fftshift(fftn(ifftshift(img1.*mask))); 
img2 = fftshift(fftn(ifftshift(img2.*mask)));


% perform radial projection
nrQ = im_w*p/2;
if mod(nrQ,1)>0
    q = -nrQ:1:nrQ;
else
    q = -nrQ:1:nrQ-1;
end
[x,y] = meshgrid(q,q);
q_radius = sqrt(x.^2+y.^2);
[q_radius,idx_img_sorted] = sort(q_radius(:)); %[q_radius,img(idx_img_sorted)] would consist of pairs [R,img(R)] which have to be summed 

if isempty(idx_radial_search)
    q_search = (1:1:sqrt(2)*nrQ).';
    if q_search(end)>=q_radius(end)
        %no values outside of image
        q_search = q_search(1:end-1);
    end

    idx_radial_search = zeros(numel(q_search),1);
    % we will now find the index of Q corresponding to the index of
    % cumsum(img(idx_img)) where the cumulative sum value corresponds to the
    % sum over all pixels inside the search circle of radius q_search
    for iQ=1:numel(q_search)
        idx_radial_search(iQ) = find(q_radius>=q_search(iQ),1,'first')-1; %where is the last index of q_radius<q_search(iQ)?
    end
    if idx_radial_search(end)<size(q_radius,1)
        % include very last value
        idx_radial_search = [idx_radial_search;size(q_radius,1)];
    end
end


%calculate frc
img_prod = real(img1.*conj(img2));
img_prod_cumsum = cumsum(img_prod(idx_img_sorted),1);
numerator = [img_prod_cumsum(1);diff(img_prod_cumsum(idx_radial_search))];

img1_cumsum = cumsum(abs(img1(idx_img_sorted)).^2);
img1_cumsum = [img1_cumsum(1);diff(img1_cumsum(idx_radial_search))];
img2_cumsum = cumsum(abs(img2(idx_img_sorted)).^2);
img2_cumsum = [img2_cumsum(1);diff(img2_cumsum(idx_radial_search))];

%frc = sum_circle(img*img_true)/sqrt(sum_circle(|img|^2)*sum_circle(|img_true|^2))
frc_curve = numerator./sqrt(img1_cumsum.*img2_cumsum);
q = (0:1:numel(frc_curve)-1).'/(im_w*px_size); %spatial frequency in units of 1/nm


if isfield(options,'estimateSpurious') && options.estimateSpurious
    L = im_w*px_size;
    v = 1./(2*pi*q*L).*numerator./(sinc(pi*q*L).^2);
    
    sigma_m = 0:50; %loc uncertainty mean in nm
    sigma_s = 0:1:10; %stdev in nm
    deriv_frc = zeros(numel(sigma_m),numel(sigma_s));
    %brute force search for widest plateau corresponding to smallest mean
    %derivative
    for iMean=sigma_m
        for jStdev=sigma_s
            sigma = 1+8*pi^2*jStdev^2*q.^2;
            h = 1./sqrt(sigma).*exp(-4*pi^2*iMean^2*q.^2./sigma);
            v_temp = log(abs(v./h));v_temp(1) = v_temp(2);
            v_temp = smooth(v_temp,0.05,'loess');
            
            deriv_frc(iMean-min(sigma_m)+1,jStdev-min(sigma_s)+1) = mean(abs(diff([v_temp(3:end),v_temp(1:end-2)],1,2))); %central difference
        end
    end
    [~,idxMin] = min(deriv_frc(:));
    [min_row,min_col] = ind2sub(size(deriv_frc),idxMin);
    sigma_m = sigma_m(min_row); sigma_s = sigma_s(min_col);
    
    sigma = 1+8*pi^2*sigma_s^2*q.^2;
    h = 1./sqrt(sigma).*exp(-4*pi^2*sigma_m.^2*q.^2./sigma);
    v_plot = smooth(log10(abs(v./h)),0.10,'rloess');
    
    %NQ/4 = minimum of curve, sigma_m and sigma_s known
    %now subtract from numerator, add to denumerator
    
end
    
    

end