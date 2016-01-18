
%% Set options
optionsSim.ampMean = 300; optionsSim.ampSigma = 40; % peak amplitude above background in photons
optionsSim.background = 1e3; % background count in photons
optionsSim.PSFsigma = 1.25; % PSF stdev in pixels
optionsSim.halfwMask = 7; % half window size of PSF mask
optionsFit.usePixelIntegratedFit = false; %use pixel integrated fit, not possible when fitting angles


%% Run script
amp_mean = optionsSim.ampMean;
amp_sigma = optionsSim.ampSigma;
PSFsigma = optionsSim.PSFsigma;
background = optionsSim.background;
halfw_mask = optionsSim.halfwMask;

usePixelIntegratedFit = optionsFit.usePixelIntegratedFit;


load('testrun_positions.mat','candidatePos');
zPos_range = -500:10:500;
movie = zeros(64,64,numel(zPos_range));
fitData_truth = cell(numel(zPos_range),1);

for zIdx = 1:numel(zPos_range)
    zPos = zPos_range(zIdx);
    sigma_x = 1.0+(zPos-250).^2*2.0/1e3^2;
    sigma_y = 1.0+(zPos+250).^2*2.0/1e3^2;
    
    img = zeros(64,64);
    ground_truth = zeros(size(candidatePos,1),7);
    for iPos = 1:size(candidatePos,1);
        r = 0.95+0.1*rand(1,2);
        rpos = -0.5+rand(1,2);
        angle = 0;
        sigma_x = sigma_x*r(1);
        sigma_y = sigma_y*r(2);
        amp = amp_mean + amp_sigma*randn(1);
        
        mask = gaussianMaskElliptical(rpos(1),rpos(2),amp,sigma_x,sigma_y,angle/360*2*pi,halfw_mask,usePixelIntegratedFit);
        img(candidatePos(iPos,2)-halfw_mask:candidatePos(iPos,2)+halfw_mask,candidatePos(iPos,1)-halfw_mask:candidatePos(iPos,1)+halfw_mask) = img(candidatePos(iPos,2)-halfw_mask:candidatePos(iPos,2)+halfw_mask,candidatePos(iPos,1)-halfw_mask:candidatePos(iPos,1)+halfw_mask)+mask;
        
        ground_truth(iPos,:) = [candidatePos(iPos,:)+rpos,amp,background,sigma_x,sigma_y,zPos];
    end
    movie(:,:,zIdx) = poissrnd(img+background);
    fitData_truth(zIdx) = {ground_truth};
    
end

save('calibrationFileAstigmatism.mat','movie','fitData_truth','optionsSim');