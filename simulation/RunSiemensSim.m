% RunSiemensSim: create a STORM simulation movie of a Siemens star structure

%% Set params
optionsSim.frames = 1000;
optionsSim.density = 750; %per micron
optionsSim.arms = 24; %siemens star arms
optionsSim.NA = 1.4;
optionsSim.fullActivationOnce = false; %if an emitter is active, it will have maximum signal (set below) for one single frame only. Useful for SNR/accuracy test

optionsCamera.pixelSize = 100; %in nm
optionsCamera.fov = [512;512];
optionsCamera.readNoise = 0; %ADU counts of CCD readout noise standard deviation
optionsCamera.gain = 1; %EM gain
optionsCamera.sensitivity = 1; %how many electrons per pixel
optionsCamera.bias = 0; %how many counts are added to image in the end to avoid negative numbers

optionsPhoton.lambda = 670; %[nm]
optionsPhoton.amp = 50; %average number of photons for a full frame
optionsPhoton.bg = 50;
optionsPhoton.kOff = 1; %per frame
optionsPhoton.kBleach = 0.1; %per frame
optionsPhoton.kOn = 0.25*optionsPhoton.kOff/optionsSim.density; %slow density situation


%% Run simulation
% for i=1:numel(amp)
%     optionsPhoton.amp = amp(i);
%     load('ground_truth_noactivationonce.mat');
%     fitData_truth = cellfun(@(var) [var(:,1:3),var(:,4)./(10+3*(20))*optionsPhoton.amp,var(:,5)/10*optionsPhoton.bg,var(:,6:end)],fitData_truth,'UniformOutput',false);
% %     [movie, fitData_truth] = SiemensSim(optionsSim,optionsCamera,optionsPhoton);
%     movie = createMovie(fitData_truth,optionsSim,optionsCamera,optionsPhoton);
%     
%     save(['simulation_star_snr_',num2str(optionsPhoton.amp/sqrt(optionsPhoton.bg),3),'.mat'],'fitData_truth','optionsSim','optionsCamera','optionsPhoton');
%     save_tiff(['simulation_star_snr_',num2str(optionsPhoton.amp/sqrt(optionsPhoton.bg),3),'.tif'],movie,'uint',16,true);
% end