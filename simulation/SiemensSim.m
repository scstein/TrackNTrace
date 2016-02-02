function [movie, fitData_truth] = SiemensSim(optionsSim,optionsCamera,optionsPhoton)
% [movie, fitData_truth] = SiemensSim(optionsSim,optionsCamera,optionsPhoton) 
% Create STORM localization movie simulation on a Siemens star.
% 
% INPUT:
%     optionsSim: Simulation options struct, see RunSiemensSim for details.
%     
%     optionsCamera: EMCCD options struct, see RunSiemensSim for details.
%     
%     optionsPhoton: Fluorophore photon yield options struct, see
%     RunSiemensSim for details.
%     
% 
% OUTPUT:
%     movie: [Y,X,F] stack of F images with width X and height Y. Format is
%     UINT16.

fov = optionsCamera.fov;
[x,y] = meshgrid(-floor(fov(1)/2):floor(fov(1)/2)-1,-floor(fov(2)/2):floor(fov(2)/2)-1);
imgMask = cos(atan2(y,x)*optionsSim.arms)>0;
psf_sigma = 0.3*optionsPhoton.lambda/optionsSim.NA/optionsCamera.pixelSize;
halfw = round(3*psf_sigma);
imgMaskBorder = zeros(size(imgMask));
imgMaskBorder(halfw+2:fov(2)-halfw-1,halfw+2:fov(1)-halfw-1) = 1;
imgMask = imgMask & imgMaskBorder; %cut off border of mask

[idxRow,idxCol] = ind2sub(size(imgMask),find(imgMask));
nrMaskPixels = size(idxRow,1);
nrPos = round(nrMaskPixels*(optionsCamera.pixelSize/1e3)^2*optionsSim.density); %calculate total number of emitters
posIdx = randi(nrMaskPixels,[nrPos,1]);
posRow = idxRow(posIdx); posCol = idxCol(posIdx);
posRow = posRow+rand(nrPos,1)-0.5; posCol = posCol+rand(nrPos,1)-0.5;

[fitData_truth] = createStates([posCol,posRow],optionsSim,optionsCamera,optionsPhoton);

movie = createMovie(fitData_truth,optionsSim,optionsCamera,optionsPhoton);

end


function [fitData_truth] = createStates(pos,optionsSim,optionsCamera,optionsPhoton)

nrParticles = size(pos,1);
nrFrames = optionsSim.frames;
psf_sigma = 0.3*optionsPhoton.lambda/optionsSim.NA/optionsCamera.pixelSize;

%Markov chain model according to "Probability-based particle detection that
%enables threshold-free and robust in vivo single molecule tracking", Smith
%et al, 2015

k_deactivate = optionsPhoton.kOff+optionsPhoton.kBleach; %rate of molecule deactivation, either by going into off state or bleaching
k_activate = 1/(1/k_deactivate+1/optionsPhoton.kOn); %rate of molecule transition to active cycle
p_bleach = optionsPhoton.kBleach/k_deactivate; %probability of bleaching event occuring

nr_events_activate = poissrnd(k_activate*nrFrames,[nrParticles,1]); %activation events drawn from poissonian distribution because time between events is exponentially distributed
nr_events_to_bleach = 1+geornd(p_bleach,[nrParticles,1]); %for normal organic fluorophores: bleaching distribution is geometric. Could also be exponential.
nr_events = min(nr_events_activate,nr_events_to_bleach);

particle_data = cell(nrParticles,1);
for iParticle = 1:nrParticles
    if nr_events(iParticle)<1
        continue
    end
    
    if optionsSim.fullActivationOnce
        % let particle "blink" once and die
        nr_events(iParticle) = 1;
        nrFrames_active = 1;
        active_time = 1;
    else
        % Get frames of activation events and the activation starting time plus
        % duration
        start_time = rand(nr_events(iParticle),1);
        active_time = exprnd(1/k_deactivate,[nr_events(iParticle),1]);
        
        % Determine [frame,time of activation in frame] vector
        nrFrames_active = ceil(start_time+active_time);
    end
    frames_activation = sort(randperm(nrFrames,nr_events(iParticle)));
    all_frames_active = []; %[x y frame amp]
    
    for jEvent = 1:nr_events(iParticle)
        if nrFrames_active(jEvent)>1
            time_vec = [(frames_activation(jEvent):frames_activation(jEvent)+nrFrames_active(jEvent)-1).', ...
            [1-start_time(jEvent);ones(nrFrames_active(jEvent)-2,1);rem(start_time(jEvent)+active_time(jEvent),1)]];
        else
            time_vec = [frames_activation(jEvent),active_time(jEvent)];
        end
        all_frames_active = [all_frames_active;[repmat(pos(iParticle,:),nrFrames_active(jEvent),1), time_vec]];
    end
    
    % sort out non-unique solutions (unlikely but possible) and save
    [~,idx,~] = unique(all_frames_active(:,3),'stable');
    particle_data(iParticle) = {all_frames_active(idx,:)};     
end

% finally save fitData array
particle_data = vertcat(particle_data{:});
fitData_truth = cell(nrFrames,1);

for iFrame = 1:nrFrames
    idx_frame = particle_data(:,3)==iFrame;
    nrParticle = sum(idx_frame);
    particle_data_frame = [particle_data(idx_frame,1:2), zeros(nrParticle,1), particle_data(idx_frame,4)*optionsPhoton.amp, repmat(optionsPhoton.bg,nrParticle,1), repmat(psf_sigma,nrParticle,1), ones(nrParticle,1)];
    [~,idx_sort] = sort(particle_data_frame(:,1),'ascend');
    fitData_truth(iFrame) = {particle_data_frame(idx_sort,:)};
    %[x y z amp bg sigma exitflag] where z = 0 and exitflag = 1
end

end %end of createStates
    

function movie = createMovie(fitData_truth,optionsSim,optionsCamera,optionsPhoton)

psf_sigma = 0.3*optionsPhoton.lambda/optionsSim.NA/optionsCamera.pixelSize;
halfw = round(3*psf_sigma);

fov = optionsCamera.fov;
movie = zeros(fov(2),fov(1),optionsSim.frames);
for iFrame=1:numel(fitData_truth)
    if isempty(fitData_truth{iFrame})
        continue
    end
    
    data_frame = fitData_truth{iFrame}(:,[1 2 4 6]);
    img_truth = zeros(fov(2),fov(1));
    for jPos=1:size(data_frame,1)
        pos_idx = floor(data_frame(jPos,1:2));
        [img] = gaussianMaskElliptical(rem(data_frame(jPos,1),1),rem(data_frame(jPos,2),1),data_frame(jPos,3),data_frame(jPos,4),data_frame(jPos,4),0,halfw,true);
        img_truth(pos_idx(2)-halfw:pos_idx(2)+halfw,pos_idx(1)-halfw:pos_idx(1)+halfw) = img_truth(pos_idx(2)-halfw:pos_idx(2)+halfw,pos_idx(1)-halfw:pos_idx(1)+halfw)...
            +img;
    end
    
%     img_truth = poissrnd(img_truth+optionsPhoton.bg) *optionsCamera.gain/optionsCamera.sensitivity + optionsCamera.bias + randn(size(img_truth))*optionsCamera.readNoise;
    movie(:,:,iFrame) = img_truth;
end

movie = poissrnd(movie+optionsPhoton.bg) *optionsCamera.gain/optionsCamera.sensitivity + optionsCamera.bias + randn(size(movie))*optionsCamera.readNoise;
movie(movie<0) = 0;
movie = uint16(movie);

end
        
        