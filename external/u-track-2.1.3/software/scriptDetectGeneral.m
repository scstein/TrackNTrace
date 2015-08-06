movieParam.directory = 'U:\Jan Thiart\dicty\2015-05-19\'; %string, directory where movie is situated
% movieParam.filename = 'bilayer1_fastest_1.spe'; %string, movie filename
% movieParam.dfilename = 'dark1.spe'; %string, movie filename
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 5000; %number of last image in movie
movieParam.format = 'tif';
files_list = dir(movieParam.directory);
files_list = files_list(find(cellfun(@(var) ~isempty(strfind(var,'tif')) && isempty(strfind(var,'txt')) && isempty(strfind(var,'xml')),{files_list(:).name})));

matlabpool('open','Skynet',80,'AttachedFiles',{'xrepmat.m','repmat.mexw64'});
%% detection parameters
detectionParam.border = 0; %integer, crop images horizontally by <border> pixels left and right
detectionParam.psfSigma = 1; %double, 0.21*(emission wavelength)/(numerical aperture) in [px], for maxima detection
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 14; %Camera bit depth
detectionParam.alphaLocMax = 0.05; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.calcMethod = 'g';
detectionParam.background = [];

saveResults.dir = 'U:\Jan Thiart\dicty\2015-05-19\'; %directory where to save input and output
% saveResults.filename = 'testing.mat'; %name of file where input and output are saved

%verbose state
verbose = 0;


%% Prepare movie
for i=1:size(files_list,1)
    
    movieParam.filename = files_list(i).name;
    saveResults.filename = [movieParam.filename(1:end-4),'_res'];
    if strcmp(movieParam.format,'tif')
        rawim = double(readTiff([movieParam.directory,movieParam.filename]));
        rawim = rawim(:,detectionParam.border+1:end-detectionParam.border,:);
        %     dark = readTiff([movieParam.directory,movieParam.dfilename]);
        %     dark = dark(:,detectionParam.border+1:end-detectionParam.border,:);
    elseif strcmp(movieParam.format,'h5')
        rawim = hdf5read([movieParam.directory,movieParam.filename],'data');
        rawim = permute(rawim,[2 1 3]);
        rawim = rawim(:,detectionParam.border+1:end-detectionParam.border,:);
        rawim = double(rawim);
        %     dark = hdf5read(dfilename,'data');
        %     dark = permute(dark,[2 1 3]);
        %     dark = dark(:,border+1:end-border,:);
    elseif strcmp(movieParam.format,'spe')
        rawim = double(readSPE([movieParam.directory,movieParam.filename]));
        dark = double(readSPE([movieParam.directory,movieParam.dfilename]));
        rawim = rawim(:,detectionParam.border+1:end-detectionParam.border,:);
        dark = dark(:,detectionParam.border+1:end-detectionParam.border,:);
    end

%calculate correction image totm from dark images and correct rawim and bgim
% totm = CalculateDark(dark);
% 
% for j=1:size(rawim,3)
%     rawim(:,:,j) = rawim(:,:,j)+totm;
% end
% clear dark

% run the detection function

% matlabpool local
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(rawim(:,:,movieParam.firstImageNum:movieParam.lastImageNum),movieParam,detectionParam,verbose);
% matlabpool close
%save results
if isstruct(saveResults)
    save([saveResults.dir saveResults.filename],'movieParam','detectionParam',...
        'movieInfo','exceptions','localMaxima','background','psfSigma');
end

end
