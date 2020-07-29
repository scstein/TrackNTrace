function [plugin] = plugin_defocusedRefinement()

% Add neccessary path for plugin to run
warning off
addpath('DefocusedPatterns_Code/');
addpath('DefocusedRefinement_Code/');
warning on
%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Defocused refinement';

% Type of plugin.
% 1: Candidate detection
% 2: Refinement
% 3: Tracking
type = 2;

% The functions this plugin implements
mainFunc =  @refine_candidates;

% Description of output parameters
outParamDescription = {'x';'y';'Integrated intensity';'Background';'in-plane angle [deg]';'out-of-plane angle [deg]'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Additional function to call before and after main function
plugin.initFunc = @copy_options;
plugin.postFunc = @assemble_pattern_movie;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Refinement for defocused emitters based on a conjugate gradient algorithm minimizing a least squares cost function of the difference of the actual image to theoretical template patterns.\n \nTemplate patterns are calculated using wave optical calculations as described in M. B?hmer and J. Enderlein, "Orientation imaging of single molecules by wide-field epifluorescence microscopy" J. Opt. Soc. Am. B 20 , 554 (2003).';

% Deactivate TNT's parallel processing, the findCandidates function is parallelized itself
plugin.useParallelProcessing = false;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
plugin.add_text('This plugin requires the ''Defocused patterns'' plugin to be used for candidate search. All relevant setup parameters (like magnification, NA, etc.) are copied over from there.','center','bold');
plugin.add_param('iterations',...
    'int',...
    {20,1,inf},...
    'Number of iterations the optimizer uns.');
plugin.add_param('gradientStep',...
    'float',...
    {1e-12,0,inf},...
    'The gradient step used when searching minima along the gradient line. \nFor too low values the optimization will not work well, for too high values the parameters explode. \nUse this in conjunction with the ''show_optimization'' flag to check what is a good value.');
plugin.newRow();
plugin.add_param('show_optimization',...
    'bool',...
    false,...
    'Visualize the optimization procedure for every candidate.');
plugin.add_param('generate_pattern_movie',...
    'bool',...
    false,...
    'Compute a movie of the detected patterns. The result is displayed and saved to the ''options.pattern_movie'' variable.');
end


%   -------------- User functions --------------

function [options] = copy_options(options)
% Copy parameters from the "Defocused patterns" candidate search plugin,
% which is required for this plugin to work.
global candidateOptions;

% Setup parameters
mag = candidateOptions.mag; % magnification
NA = candidateOptions.NA; % objective NA
lamem = candidateOptions.lambdaEmission/1000;   % emission wavelength in ?m (options.lambda_emission in nm)
pixelsize = candidateOptions.pixelSize;  % camera pixelsize in ?m
nn = candidateOptions.patternRadius;         % pixel radius of disk pattern is calculated for
defocus = candidateOptions.defocus; % focus position
    
% -- Compute electromagnetic field used later to calculate the defocused patterns --
[options.SEPdata, options.SetupParams] = SEPDipoleSimple(nn, pixelsize, NA, lamem, mag, defocus);

% Size of generated pattern image is 2*field+1
% -> This should be the size of the area in the image we compare the
% pattern to. Its not the size of the pattern! (which is defined by nn)
options.field = [size(candidateOptions.patterns,1),size(candidateOptions.patterns,2)];
% ---
       
% Create a temporary folder to save images in
if options.generate_pattern_movie
    time = clock;
    timestamp = sprintf('%i-m%02i-d%02i-%02ih%02i-%02is',time (1),time(2),time(3),time(4),time(5), time(6));
    options.temp_dir = sprintf('temp_patterns_refine_%s',timestamp);
    mkdir(options.temp_dir);
end

end

function [refinementData] = refine_candidates(img,candidateData,options,currentFrame)
% Wrapper function for psfFit_Image function (see below). Refer to
% tooltips above and to psfFit_Image help to obtain information on input
% and output variables.
%
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidateData: output of the "Defocused patterns" candidate search
%     plugin
%
%     options: Struct of input parameters provided by GUI.
%
% OUTPUT:
%     refinementData: 1x1 cell of 2D double array of fitted parameters
%     [x,y,z,A,B,[other parameters]]. Other parameters can be q1, q2, q3
%     (refer to locateParticles.m or to TrackNTrace manual for more
%     information). q_i will be calculated back to sigma_x,sigma_y,
%     rotation angle and possibly z in post-processing function (see
%     below).

nrCandidates = size(candidateData,1);

% reserve memory for output
refinementData = zeros(size(candidateData));

% Edge length of rectangular window used for optimization is 2*hwSize+1
hwSize = floor(options.field(1)/2);

im_detect = zeros(size(img)); % image to save detected theoretical patterns to
for iCand = 1:nrCandidates    
    
    
    % Initial values
    x0 = zeros(6,1);
    x0(1) = 0;  % xcent (coordinate system 0 at the center pixel and we make a cutout around the candidate position)
    x0(2) = 0;  % ycent (coordinate system 0 at the center pixel and we make a cutout around the candidate position)
    x0(3) = candidateData(iCand,5)*pi/180; % be/phi (in plane angle). Why do we need a minus here ???
    x0(4) = candidateData(iCand,6)*pi/180; % al/theta (out of plane angle)
    x0(5) = candidateData(iCand,3);        % Amp
    x0(6) = candidateData(iCand,4);        % BG
    
    %% Pack Setup parameters into a struct for better handling   
    pars.img = img(candidateData(iCand,2)-hwSize:candidateData(iCand,2)+hwSize, candidateData(iCand,1)-hwSize:candidateData(iCand,1)+hwSize);
    pars.field = [hwSize,hwSize];
    pars.SetupParams = options.SetupParams;
    pars.SEPdata = options.SEPdata;
    
    %% Conjugate Gradient search    
    
    % Conjugate Gradient Options
%     fname = 'MLE_PATT_FUNC_SEP';
    fname = 'LSQ_PATT_FUNC_SEP';    
    gname = []; % Use numerical gradient
    
    FTOL = -1;
    ITMAX = options.iterations;
    
%     ITFUNC = [];
% ---- Uncomment to visualize the optimization process (only useful if plugin is not paralellelized)----
if options.show_optimization
    ITFUNC = 'ITFUNC_PATT_VIS_SEP'; % Visualization after every iteration    
    h = figure;
    mim(pars.img);
    set(h,'Position',[300,300,500,500]);    
    h2 = figure;
    set(h2,'Position',[800,300,500,500]);
end
% ----
    verbose = false;
    [ x_fin, ~ ] = conjugateGradient( x0, fname, gname, FTOL, options.gradientStep, ITMAX, ITFUNC, verbose, pars);
% ----
if options.show_optimization
    close(h);
    close(h2);
end
% ----

%     lsqOptions = optimoptions('lsqnonlin');
%     lsqOptions.MaxIter = 100;
%     lsqOptions.MaxFunEvals = 100;
%     x_fin = lsqnonlin(fname, x0,[],[],lsqOptions,pars);

    refinementData(iCand,1) =  x_fin(1) + candidateData(iCand,1);
    refinementData(iCand,2) =  x_fin(2) + candidateData(iCand,2);
    refinementData(iCand,3) =  x_fin(5); % Amp
    refinementData(iCand,4) =  x_fin(6); % BG
    refinementData(iCand,5) =  x_fin(3)*180/pi; %'in-plane angle [deg]'
    refinementData(iCand,6) =  x_fin(4)*180/pi; %'out of plane angle [deg]'


    im_detect(candidateData(iCand,2)-hwSize:candidateData(iCand,2)+hwSize, candidateData(iCand,1)-hwSize:candidateData(iCand,1)+hwSize) = ...
        im_detect(candidateData(iCand,2)-hwSize:candidateData(iCand,2)+hwSize, candidateData(iCand,1)-hwSize:candidateData(iCand,1)+hwSize) + ...
        x_fin(5)*PatternGenerateSimple_SEP( x_fin(1), x_fin(2), x_fin(3), x_fin(4), pars.field, options.SetupParams, options.SEPdata);
end
    
if options.generate_pattern_movie
    filename = [options.temp_dir, filesep, sprintf('%05d.tif',currentFrame)]; %#ok<NASGU>
    evalc('save_tiff(filename, im_detect)');
end

end

function [candidateData, options] = assemble_pattern_movie(candidateData, options)

if options.generate_pattern_movie
    fprintf(' Defocused Refinement: Generating movie of detected patterns.\n')
    
    img_files = dir([options.temp_dir, filesep, '*.tif']);
    
    global movie;
    options.pattern_movie = single(zeros(size(movie{1})));
    
    for iImg = 1:size(img_files,1)
        %Using evalc to suppress output to terminal
        [~,img] = evalc(['read_tiff(''',options.temp_dir, filesep, img_files(iImg).name,''',false)'] );
        options.pattern_movie(:,:,iImg) = single(img);
    end
    
    global filename_movie;
    [~,movie_name,ext] = fileparts(filename_movie);
    movie_name = [movie_name,ext];
    h = TNTvisualizer(options.pattern_movie);
    set(h, 'name', sprintf('''%s'' - Refined patterns',movie_name));
    
    rmdir(options.temp_dir,'s');
end
end

