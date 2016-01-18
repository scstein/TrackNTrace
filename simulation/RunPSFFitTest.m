%RunPSFFitTest


%% Set parameters
optionsSim.ampMean = 300; optionsSim.ampSigma = 100; % peak amplitude above background in photons
optionsSim.background = 1e3; % background count in photons
optionsSim.PSFsigma = 1.25; % PSF stdev in pixels
optionsSim.halfwMask = 7; % half window size of PSF mask
optionsSim.sigmaRange = [optionsSim.PSFsigma-0.6,optionsSim.PSFsigma+0.6];
optionsSim.thetaRange = 0; % rotation of Gaussian function. Set to 0 to simulate perfect axis alignment

% fit parameters
optionsFit.usePixelIntegratedFit = false; %use pixel integrated fit, not possible when fitting angles
optionsFit.useMLE = false; % use MLE
optionsFit.halfw = 5; %half window size of fit mask
optionsFit.varsToFit = [ones(7,1)]; %fit parameters [x,y,A,BG,sigma_x,sigma_y,angle] true/false

plotFits = true;

%% Create simulated image with predetermined positions
[candidatePos,img,img_truth,ground_truth] = createPSFFitTest(optionsSim,optionsFit,plotFits);


%% Test performance
% [params] = psfFit_Image( img, candidatePos.',optionsFit.varsToFit,optionsFit.usePixelIntegratedFit,optionsFit.useMLE,optionsFit.halfw,optionsSim.PSFsigma);
% fitData = {params(:,params(end,:)==1).'};
% 
% for iPos = 1:size(fitData{1},1)
%     a = fitData{1}(iPos,5); b = fitData{1}(iPos,7); c = fitData{1}(iPos,6);
%     angle_fit = 0.5*atan(2*b/(c-a));
%     if abs(angle_fit)<=pi/4
%         sigma_x_fit = 1/sqrt(a+c-2*b/sin(2*angle_fit));
%         sigma_y_fit = 1/sqrt(a+c+2*b/sin(2*angle_fit));
%     else
%         sigma_y_fit = 1/sqrt(a+c-2*b/sin(2*angle_fit));
%         sigma_x_fit = 1/sqrt(a+c+2*b/sin(2*angle_fit));
%         disp('This case should not happen!')
%     end
%     
%     if angle_fit==0 ||sum(isnan([sigma_x_fit,sigma_y_fit]))>0
%         sigma_x_fit = 1/sqrt(2*a);
%         sigma_y_fit = 1/sqrt(2*c);
%     end
%     
%     fitData{1}(iPos,5:7) = [sigma_x_fit,sigma_y_fit,angle_fit];
% end
% figure; mim(img); hold on; plot(fitData{1}(:,1),fitData{1}(:,2),'bo'); hold off;

tic;
for i=1:1e4
[params] = psfFit_Image( img, candidatePos.',optionsFit.varsToFit,optionsFit.usePixelIntegratedFit,optionsFit.useMLE,optionsFit.halfw,optionsSim.PSFsigma);
end
toc
