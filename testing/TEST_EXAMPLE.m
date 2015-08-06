movie = read_tiff('U:\Weixing\weixing\101 polar origami\bleaching-3.7exci-roi2.tif', false);

%% Configure options
candidateOptions.sigma = 1.;
candidateOptions.corrThresh = 0.8;
candidateOptions.windowSize = 10;
candidateOptions.fitForward = true;

fittingOptions.fitSigma = true;
fittingOptions.usePixelIntegratedFit = true;
fittingOptions.useMLErefine = true;

% Compute the drift
fitData = fitMovie(movie, candidateOptions, fittingOptions);
