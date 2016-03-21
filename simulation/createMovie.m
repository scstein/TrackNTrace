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

for iChunk=1:ceil(size(movie,3)/2e2)
    movie(:,:,1+(iChunk-1)*2e2:min(iChunk*2e2,size(movie,3))) = poissrnd(movie(:,:,1+(iChunk-1)*2e2:min(iChunk*2e2,size(movie,3)))+optionsPhoton.bg) *optionsCamera.gain/optionsCamera.sensitivity + optionsCamera.bias;% + randn(size(movie))*optionsCamera.readNoise;
end
movie(movie<0) = 0;
movie = uint16(movie);

end