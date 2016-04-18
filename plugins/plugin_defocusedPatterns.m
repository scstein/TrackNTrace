function [plugin] = plugin_defocusedPatterns()

% Add neccessary path for plugin to run
warning off
addpath('DefocusedPattern_Code/');
warning on
%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Defocused patterns';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
type = 1;

% The functions this plugin implements
mainFunc =  @findCandidates_defocusedPatterns;

% Description of output parameters
outParamDescription = {'x';'y';'Integrated intensity';'Background';'in-plane angle [deg]';'out-of-plane angle [deg]'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Additional function to call before and after main function
plugin.initFunc = @calculate_patterns;
plugin.postFunc = @assemble_pattern_movie;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Candidate detection for defocused emitters based on normalized cross correlation of the input images with theoretical patterns.\n \nTemplate patterns are calculated using wave optical calculations as described in M. Böhmer and J. Enderlein, "Orientation imaging of single molecules by wide-field epifluorescence microscopy" J. Opt. Soc. Am. B 20 , 554 (2003).';

% Deactivate TNT's parallel processing, the findCandidates function is parallelized itself
plugin.useParallelProcessing = false;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
plugin.add_param('mag',...
    'float',...
    {400, 0, inf},...
    'Objective magnification');
plugin.add_param('NA',...
    'float',...
    {1.4, 0, inf},...
    'Objective numerical aperture.');
plugin.add_param('lambdaEmission',...
    'float',...
    {680, 0, inf},...
    'Emission wavelength of imaged fluorophore in unit [nanometer].');
plugin.add_param('pixelSize',...
    'float',...
    {24, 0, inf},...
    'Camera sensor pixel size in unit [µm].');
plugin.add_param('patternRadius',...
    'int',...
    {30, 0, inf},...
    'Pixel radius of disk pattern is calculated for.');
plugin.add_param('defocus',...
    'float',...
    {1, 0, inf},...
    'Focus distance to the molecule plane patterns are calculated for in [µm]');
plugin.newRow();
plugin.add_param('inplaneAngleStep',...
    'float',...
    {5, 0, 360},...
    'Step size for in-plane angle patterns are calculated for in [degree]');
plugin.add_param('outplaneAngleStep',...
    'float',...
    {10, 0, 90},...
    'Step size for out-of-plane angle patterns are calculated for in [degree]');
plugin.add_text('--------------------------------------------------------------','center');
plugin.add_param('CorrThreshold',...
    'float',...
    {0.3, 0, 1},...
    'Correlation threshold for pixels to count as emitter candidates. \nMust be between 0 and 1.');
plugin.newRow();
plugin.add_param('show_patterns',...
    'bool',...
    false,...
    'Display the computed patterns.');
plugin.add_param('generate_pattern_movie',...
    'bool',...
    false,...
    'Compute a movie of the detected patterns. The result is displayed and saved to the ''options.pattern_movie'' variable.');
plugin.add_param('optimize_for',...
    'list',...
    {'Memory','Speed'},...
    ['Algorithm to use for pattern matching.\n',...
    '\t \t \t Memory (recommended): Low memory consumption, but each fft for each pattern is computed for each frame. Optimal for single/few frames. \n \n', ...
    'Speed: Pattern ffts are precomputed. Choose this only if you have longer movies with small image size or very few patterns. \n',...
    'Memory consumption is roughly 2*(img_size+2*patternRadius)*nrPatterns*64/1e9 Gigabytes. \n',...
    'This is 100 GB for a 1024x1024 movie with 730 patterns (5 degree inplane step, 10 degree outplane step).']);
end


%   -------------- User functions --------------

function [options] = calculate_patterns(options)

global movie;

% Input variables to check if they changed
% Otherwise we do not compute the patterns but take the old ones
persistent old_mag old_NA old_lamem old_pixelSize old_nn old_defocus old_inplaneAngleStep old_outplaneAngleStep %#ok<USENS>
% Output that we store as backup
persistent old_patterns old_patt_angle_inplane old_patt_angle_outplane old_patt_fft_rot90;

% Setup parameters
mag = options.mag; % magnification
NA = options.NA; % objective NA
lamem = options.lambdaEmission/1000;   % emission wavelength in µm (options.lambdaEmission in nm)
pixelSize = options.pixelSize;  % camera pixelSize in µm
nn = options.patternRadius;         % pixel radius of disk pattern is calculated for
defocus = options.defocus; % focus position
inplaneAngleStep = options.inplaneAngleStep;
outplaneAngleStep = options.outplaneAngleStep;

% Angle values patterns are generated for
angle_inplane = 0:inplaneAngleStep:360;% in plane angle (degree)
angle_outplane = 0:outplaneAngleStep:90; % out-of-plane angles (degree)
% Size of generated pattern image is 2*field+1
field = [nn,nn];
% ---

% Check if we can reuse the last computed patterns
use_old_patterns = false;
if ~isempty(old_mag)
    if (old_mag == mag) && (old_NA == NA) && (old_lamem == lamem) && (old_pixelSize == pixelSize) && (old_nn == nn) ...
            && (old_defocus == defocus) && (old_inplaneAngleStep == inplaneAngleStep) && (old_outplaneAngleStep == outplaneAngleStep)
        
        use_old_patterns = true;
    end
end

% Copy parameters for pattern generation
old_mag = mag;          old_NA = NA;
old_lamem = lamem;      old_pixelSize = pixelSize;
old_nn = nn;            old_defocus = defocus;
old_inplaneAngleStep = inplaneAngleStep;
old_outplaneAngleStep = outplaneAngleStep;

if use_old_patterns
    options.patterns = old_patterns;
    options.patt_angle_inplane = old_patt_angle_inplane;
    options.patt_angle_outplane = old_patt_angle_outplane;
    % Delete old ffts if still in memory
    if( strcmp(options.optimize_for,'Memory'))
        options.patt_fft_rot90 = [];
    else
        options.patt_fft_rot90 = old_patt_fft_rot90;
    end
    
else
    %
    nPattRows = numel(angle_inplane);
    nPattCols = numel(angle_outplane);
    nPatterns = nPattRows*nPattCols;
    pattRowSize = 2*nn+1;
    pattColSize = 2*nn+1;
    
    % Output
    patterns = zeros( pattRowSize, pattColSize, nPatterns);
    patt_angle_inplane = zeros(nPatterns, 1);
    patt_angle_outplane = zeros(nPatterns, 1);
    
    % Compute the angles for each pattern
    for iTheta = 1:numel(angle_outplane);
        for iBeta = 1:numel(angle_inplane)
            iPattern = iBeta + (iTheta-1)*numel(angle_inplane);
            patt_angle_inplane(iPattern) = angle_inplane(iBeta);
            patt_angle_outplane(iPattern) = angle_outplane(iTheta);
        end
    end
    
    % -- Compute the defocused patterns --
    [SEPdata, SetupParams] = SEPDipoleSimple(nn, pixelSize, NA, lamem, mag, defocus);
    
    
    global parallelProcessingAvailable
    if parallelProcessingAvailable
        % Check MATLAB version
        MATLAB_2013b_or_newer = false;
        MATLABversion = strsplit(version,'.');
        if(str2double(MATLABversion(1))>=8 && str2double(MATLABversion(2))>=2) % matlabpool -> parpool in MATLAB 2013b (8.2.x) and later
            MATLAB_2013b_or_newer = true;
        end
        
        % Get number of available workers
        if MATLAB_2013b_or_newer
            p = gcp('nocreate');
            nrRunningWorkers = p.NumWorkers;
        else
            nrRunningWorkers = matlabpool('size');
        end
        
        % Estimate time needed for computing the patterns
        tic
        PatternGenerateSimple_SEP( 0, 0, 0, 0, field, SetupParams, SEPdata);
        patternTime = toc;
        fprintf(' Defocused Patterns: Generating %i patterns (using parallel processing). Estimated time needed: %3.f seconds.\n', nPatterns, nPatterns*patternTime / (nrRunningWorkers-0.25));
        
        parfor iPatt = 1:numel(patt_angle_inplane)
            patterns(:,:,iPatt) = PatternGenerateSimple_SEP( 0, 0, patt_angle_inplane(iPatt)*pi/180, patt_angle_outplane(iPatt)*pi/180, field, SetupParams, SEPdata);
        end
    else
        % Estimate time needed for computing the patterns
        tic
        PatternGenerateSimple_SEP( 0, 0, 0, 0, field, SetupParams, SEPdata);
        patternTime = toc;
        fprintf(' Defocused Patterns: Generating %i patterns. Estimated time needed: %3.f seconds.\n', nPatterns, nPatterns*patternTime);
        for iPatt = 1:numel(patt_angle_inplane)
            patterns(:,:,iPatt) = PatternGenerateSimple_SEP( 0, 0, patt_angle_inplane(iPatt)*pi/180, patt_angle_outplane(iPatt)*pi/180, field, SetupParams, SEPdata);
        end
    end
    
    
    
    % Precompute the pattern fourier transforms needed for the pattern matching
    if( strcmp(options.optimize_for,'Speed'))
        Template_size = size(patterns);
        Template_size = Template_size(1:2);
        Img_size = size(movie);
        Img_size = Img_size(1:2);
        
        outsize = Img_size + Template_size - 1;
        patt_fft_rot90 = zeros(outsize(1),outsize(2), nPatterns);
        
        tic;
        fft2(rot90(patterns(:,:,1),2),outsize(1),outsize(2));
        fftTime = toc;
        
        if parallelProcessingAvailable
            fprintf(' Defocused Patterns: Computing pattern ffts (using parallel processing). Estimated time needed: %3.f seconds.\n', nPatterns*fftTime/(nrRunningWorkers-0.25));
            parfor iPatt = 1:nPatterns
                patt_fft_rot90(:,:,iPatt) = fft2(rot90(patterns(:,:,iPatt),2),outsize(1),outsize(2));
            end
        else
            fprintf(' Defocused Patterns: Computing pattern ffts. Estimated time needed: %3.f seconds.\n', nPatterns*fftTime);
            for iPatt = 1:nPatterns
                patt_fft_rot90(:,:,iPatt) = fft2(rot90(patterns(:,:,iPatt),2),outsize(1),outsize(2));
            end
        end
    else
        patt_fft_rot90 = [];
    end
    
    % Save output
    options.patterns = patterns;
    options.patt_angle_inplane = patt_angle_inplane;
    options.patt_angle_outplane = patt_angle_outplane;
    options.patt_fft_rot90 = patt_fft_rot90;
    
    old_patterns = options.patterns;
    old_patt_angle_inplane = options.patt_angle_inplane;
    old_patt_angle_outplane = options.patt_angle_outplane;
    old_patt_fft_rot90 = options.patt_fft_rot90;
end

% Visualize patterns
if(options.show_patterns)
    global filename_movie;
    [~,movie_name,ext] = fileparts(filename_movie);
    movie_name = [movie_name,ext];
    
    h = TNTvisualizer(options.patterns);
    set(h, 'name', sprintf('''%s'' - Computed patterns', movie_name));
end


% Create a temporary folder to save images in
if options.generate_pattern_movie
    time = clock;
    timestamp = sprintf('%i-m%02i-d%02i-%02ih%02i-%02is',time (1),time(2),time(3),time(4),time(5), time(6));
    options.temp_dir = sprintf('temp_patterns_%s',timestamp);
    mkdir(options.temp_dir);
end

end

function candidatePos = findCandidates_defocusedPatterns(img, options, currentFrame)
%Wrapper function for candidate finding.
%
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidateOptions: Struct of input parameters provided by GUI.
%
%     currentFrame: Integer, current movie frame in main loop.
%
% OUTPUT:
%     candidatePos - 2D array of Nx2 matrix of particle candidate positions
%     [column pixel, row pixel] without subpixel position. Middle of upper
%     left pixel would be [1,1].

nPatterns = size(options.patterns,3);
template_size = size(options.patterns);
template_size = template_size(1:2);

CORR_THRESH = options.CorrThreshold;
MIN_DIST = round(max(template_size)/2);


img_size = size(img);
outsize = img_size + template_size - 1;
img_FFT = fft2(img,outsize(1),outsize(2));

% Pattern matching
% If requested, we compute an image of the detected patterns and save it to
% a tif file. The tif files are loaded in the post processing function and
% assembled into a movie later on. If we would try to directly save to a
% MATLAB variable this is problematic for parallel processing.
compute_im_detect = false;
if options.generate_pattern_movie
    compute_im_detect = true;
end

global parallelProcessingAvailable
[match_data, im_detect] = pattern_matching(img, options.patterns, CORR_THRESH, MIN_DIST, img_FFT, options.patt_fft_rot90, compute_im_detect, parallelProcessingAvailable); %#ok<NASGU>

if options.generate_pattern_movie
    filename = [options.temp_dir, filesep, sprintf('%05d.tif',currentFrame)]; %#ok<NASGU>
    evalc('save_tiff(filename, im_detect)');
end

% Find first nonempty entry
first_found_pattern = -1;
candidatePos = [];
for iPatt = 1:nPatterns
    if isempty(match_data{iPatt})
        continue
    end
    m_data = match_data{iPatt};
    candidateData = zeros(size(m_data,1),6);
    % Convert to pixel coordinates and swap the coordinates (rows,cols to x,y)
    candidateData(:,1:2) = [m_data(:,2)+ceil(template_size(2)/2),m_data(:,1)+ceil(template_size(1)/2)];
    candidateData(:,3) = m_data(:,3); % Integrated intensity
    candidateData(:,4) = m_data(:,4); % Per pixel background
    candidateData(:,5) = options.patt_angle_inplane(iPatt);
    candidateData(:,6) = options.patt_angle_outplane(iPatt);
    
    candidatePos = candidateData;
    
    first_found_pattern = iPatt;
    break;
end

if (first_found_pattern >=1 && first_found_pattern<nPatterns)
    for iPatt = (first_found_pattern+1):nPatterns
        if isempty(match_data{iPatt})
            continue
        end
        m_data = match_data{iPatt};
        candidateData = zeros(size(m_data,1),4);
        % Convert to pixel coordinates and switch them
        candidateData(:,1:2) = [m_data(:,2)+ceil(template_size(2)/2),m_data(:,1)+ceil(template_size(1)/2)];
        candidateData(:,3) = m_data(:,3); % Integrated intensity
        candidateData(:,4) = m_data(:,4); % Per pixel background
        candidateData(:,5) = options.patt_angle_inplane(iPatt);
        candidateData(:,6) = options.patt_angle_outplane(iPatt);
        
        candidatePos = [candidatePos; candidateData];
    end
end


end

function [candidateData, options] = assemble_pattern_movie(candidateData, options)

% Clears the persistent variables used in the init function
clear plugin_defocusedPatterns

if options.generate_pattern_movie
    fprintf(' Defocused Refinement: Generating movie of detected patterns.\n')
    
    img_files = dir([options.temp_dir, filesep, '*.tif']);
    
    global movie;
    options.pattern_movie = single(zeros(size(movie)));
    
    for iImg = 1:size(img_files,1)
        %Using evalc to suppress output to terminal
        [~,img] = evalc(['read_tiff(''',options.temp_dir, filesep, img_files(iImg).name,''',false)'] );
        options.pattern_movie(:,:,iImg) = single(img);
    end
    
    global filename_movie;
    [~,movie_name,ext] = fileparts(filename_movie);
    movie_name = [movie_name,ext];
    h = TNTvisualizer(options.pattern_movie);
    set(h, 'name', sprintf('''%s'' - Detected patterns',movie_name));
    
    rmdir(options.temp_dir,'s');
end

end

