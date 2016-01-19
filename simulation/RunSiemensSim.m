% RunSiemensSim: create a STORM simulation movie of a Siemens star structure

%% Set params
optionsSim.frames = 5000;
optionsSim.density = 750; %per micron
optionsSim.arms = 6; %siemens star arms
optionsSim.NA = 1.4;
optionsSim.fullActivationOnce = true; %if an emitter is active, it will have maximum signal (set below) for one single frame only. Useful for SNR/accuracy test

optionsCamera.pixelSize = 100; %in nm
optionsCamera.fov = [64;64];
optionsCamera.readNoise = 0; %ADU counts of CCD readout noise standard deviation
optionsCamera.gain = 1; %EM gain
optionsCamera.sensitivity = 1; %how many electrons per pixel
optionsCamera.bias = 0; %how many counts are added to image in the end to avoid negative numbers

optionsPhoton.lambda = 670; %[nm]
optionsPhoton.amp = 50; %average number of photons for a full frame
optionsPhoton.bg = 10;
optionsPhoton.kOff = 1; %per frame
optionsPhoton.kBleach = 0.1; %per frame
optionsPhoton.kOn = 0.25*optionsPhoton.kOff/optionsSim.density; %slow density situation


%% Run simulation
[movie, fitData_truth] = SiemensSim(optionsSim,optionsCamera,optionsPhoton);