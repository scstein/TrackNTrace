%images start in top left corner, fitting uses cartesian coordinates -> use
%im(end:-1:1,:) to get correct correlation with angles!

function [candidatePos,img,img_truth,ground_truth] = createPSFFitTest(optionsSim,optionsFit,plotFits)

amp_mean = optionsSim.ampMean;
amp_sigma = optionsSim.ampSigma;
PSFsigma = optionsSim.PSFsigma;
sigma_aspect = PSFsigma^2;
background = optionsSim.background;
sigma_range = optionsSim.sigmaRange;
halfw_mask = optionsSim.halfwMask;
theta_range = optionsSim.thetaRange;

usePixelIntegratedFit = optionsFit.usePixelIntegratedFit;
useMLE = optionsFit.useMLE;
halfw = optionsFit.halfw;
varsToFit = optionsFit.varsToFit;

img = zeros(64,64);
ground_truth = [];
load('testrun_positions.mat','candidatePos');

for iPos = 1:size(candidatePos,1); %#ok<NODEF>
    r = -0.5+rand(1,2);
    angle = theta_range*rand(1);
    sigma_x = sigma_range(1)+diff(sigma_range)*rand(1);
    sigma_y = sigma_aspect/sigma_x;
    amp = amp_mean + amp_sigma*randn(1);
    
    mask = gaussianMaskElliptical(r(1),r(2),amp,sigma_x,sigma_y,angle/360*2*pi,halfw_mask,usePixelIntegratedFit);
    img(candidatePos(iPos,2)-halfw_mask:candidatePos(iPos,2)+halfw_mask,candidatePos(iPos,1)-halfw_mask:candidatePos(iPos,1)+halfw_mask) = img(candidatePos(iPos,2)-halfw_mask:candidatePos(iPos,2)+halfw_mask,candidatePos(iPos,1)-halfw_mask:candidatePos(iPos,1)+halfw_mask)+mask;
    
    ground_truth = [ground_truth;[candidatePos(iPos,:)+r,amp,sigma_x,sigma_y,angle]]; %#ok<AGROW>
end
img_truth = img;
img = poissrnd(img+background);


[params] = psfFit_Image( img, candidatePos.',varsToFit,usePixelIntegratedFit,useMLE,halfw,PSFsigma);
fitData = {params(:,params(end,:)==1).'};

for iPos = 1:size(fitData{1},1)
    a = fitData{1}(iPos,5); b = fitData{1}(iPos,7); c = fitData{1}(iPos,6);
    angle_fit = 0.5*atan(2*b/(c-a));
    if abs(angle_fit)<=pi/4
        sigma_x_fit = 1/sqrt(a+c-2*b/sin(2*angle_fit));
        sigma_y_fit = 1/sqrt(a+c+2*b/sin(2*angle_fit));
    else
        sigma_y_fit = 1/sqrt(a+c-2*b/sin(2*angle_fit));
        sigma_x_fit = 1/sqrt(a+c+2*b/sin(2*angle_fit));
        disp('This case should not happen!')
    end
    
    if angle_fit==0 ||sum(isnan([sigma_x_fit,sigma_y_fit]))>0
        sigma_x_fit = 1/sqrt(2*a);
        sigma_y_fit = 1/sqrt(2*c);
    end
    
    fitData{1}(iPos,5:7) = [sigma_x_fit,sigma_y_fit,angle_fit];
end

moment_truth = determineMoments(candidatePos,halfw,img_truth);

ground_truth = [ground_truth, NaN(size(candidatePos,1),1), ...
    fitData{1}(:,[1,2,3,5,6]),fitData{1}(:,7)*180/pi, NaN(size(candidatePos,1),1),...
    moment_truth];

if plotFits
    figure; mim(img); hold on; plot(fitData{1}(:,1),fitData{1}(:,2),'bo'); hold off;
end

end


function moment_truth = determineMoments(candidatePos,halfw,img_truth)
% calculate second intensity moments of image to determine single value
% decomposition for parameter initial guess

moment_truth = [];
for iPos = 1:size(candidatePos,1) 
    bg = 0; %dealing with image without background here
    m_00 = 0;
    m_20 = 0;
    m_02 = 0;
    m_10 = 0;
    m_01 = 0;
    m_11 = 0;
    
    for iCol=-halfw:halfw
        for iRow=-halfw:halfw
            img_part = img_truth(candidatePos(iPos,2)+iRow,candidatePos(iPos,1)+iCol);
            m_00 = m_00+(img_part-bg);
            m_10 = m_10+(iCol+candidatePos(iPos,1))*(img_part-bg);
            m_01 = m_01+(iRow+candidatePos(iPos,2))*(img_part-bg);
            m_20 = m_20+(iCol+candidatePos(iPos,1))^2*(img_part-bg);
            m_02 = m_02+(iRow+candidatePos(iPos,2))^2*(img_part-bg);
            m_11 = m_11+(iRow+candidatePos(iPos,2))*(iCol+candidatePos(iPos,1))*(img_part-bg);
        end
    end
    
    x_bar = m_10/m_00;
    y_bar = m_01/m_00;
    cov_matrix = [[m_20/m_00-x_bar^2,m_11/m_00-x_bar*y_bar];[m_11/m_00-x_bar*y_bar,m_02/m_00-y_bar^2]];
    lambda_1 = trace(cov_matrix)/2+sqrt(trace(cov_matrix)^2/4-det(cov_matrix));
    lambda_2 = trace(cov_matrix)/2-sqrt(trace(cov_matrix)^2/4-det(cov_matrix));
    angle_cov = 0.5*atan(2*cov_matrix(1,2)/(cov_matrix(1,1)-cov_matrix(2,2)))*180/pi;
    
    if cov_matrix(1,1)>cov_matrix(2,2) %x-axis closer to eigenvector of largest eigenvalue?
        moment_truth = [moment_truth;[sqrt(lambda_1),sqrt(lambda_2),angle_cov]]; %#ok<AGROW>
        
    else % change sign of angle and axis
%         angle_cov = -angle_cov;
        moment_truth = [moment_truth;[sqrt(lambda_2),sqrt(lambda_1),angle_cov]]; %#ok<AGROW>
        
    end
    
end

end
