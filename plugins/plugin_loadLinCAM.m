function [plugin] = plugin_loadLinCAM()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Import PHOTONS file';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
% 4: Postprocessing
% 5: Import
% 6: Export
type = 5;

% The functions this plugin implements
mainFunc =  @read_Photons;

% Description of output parameters
outParamDescription = {'FLIM movie','head'}; % set in init function

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Function to execute before frame-by-frame processing starts
% plugin.initFunc = @updateOutParamDescription;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = struct(...
    'Description','Load a PHOTONS file as produced by the LinCAM.',...
    'supportedFormats',{{'*.photons','PHOTONS'}},...
    'hasFLIM',true,...
    'hasTCSPC',true,...
    'getTCSPC',@getTCSPC,...
    'setFile',@setFile);

% If hasTCPSC is true the plugin needs an accumulate function to retrieve
% the TCSPC of the localisations. This is returned by plugin.postFunc.
% If plugin.initFunc is set the plugin supports caching via the output 
% 'index', which is passed instead of the inputfile to subsequent calls of 
% plugin.postFunc.
% The function has to called like fhdl = plugin.postFunc(pluginOptions).
% The syntax of fhdl has to be like:
% PTU_accumulate(inputfile, outputs, framebinning, pixelbinning, timegate)
plugin.initFunc = @getAccumFunc;
plugin.postFunc = @getAccumFunc;

plugin.add_param('window_x',...
    'int',...
    {1,1,inf},...
    'Boundries of an rectangular ROI.');
plugin.add_param('window_dx',...
    'int',...
    {inf,1,inf},...
    'Boundries of an rectangular ROI.');
plugin.add_param('autoWindow',...
    'bool',...
    false,...
    'Searches the position of the window automatically in each file. Usefull for batch processing');
plugin.newRow();
plugin.add_param('window_y',...
    'int',...
    {1,1,inf},...
    'Boundries of an rectangular ROI.');
plugin.add_param('window_dy',...
    'int',...
    {inf,1,inf},...
    'Boundries of an rectangular ROI.');
plugin.add_param('findWindow',...
    'button',...
    @findWindow,...
    'Tries to find the position of the window automatically.');
plugin.newRow();
plugin.add_param('pixelbinning',...
    'int',...
    {10,1,inf},...
    'The pixel binning is done after the window is applied.');
plugin.add_param('pixelsize',...
    'float',...
    {0,0,inf},...
    ['The pixel size before pixel binning in ' char(181) 'm. If unknown, set to zero.']);
plugin.add_param('cacheMovie',...
    'int',...
    {2,0,10},...
    'Number of generated movie versions to save temporary. 0 disables cache. Oldest will be replaced if necessary.');


    function setFile(file)
        if ischar(file)
            % We cannot access the private properties directly
            opts = plugin.getOptions();
            opts.importFile = file;
            if opts.autoWindow
                opts = findWindow(opts);
            end
            plugin.setOptions(opts);
        end
    end
end


%   -------------- User functions --------------


function [movie,metadata] = read_Photons(pluginOptions,filename_movie, frame_range, frame_binning,fastLT,timegate)
    if endsWith(filename_movie,'photons','IgnoreCase',true)
        if nargin<4 || isempty(frame_binning)
            frame_binning = 1;
        end
        if nargin<5 || isempty(fastLT) 
            fastLT = false;
        end
        if nargin<6
            timegate = [];
        end
        window = [pluginOptions.window_x, pluginOptions.window_x+pluginOptions.window_dx, pluginOptions.window_y, pluginOptions.window_y+pluginOptions.window_dy];
        
        % Check cache
        movieArgs = {filename_movie, frame_range, frame_binning, fastLT, timegate, window, pluginOptions.pixelbinning};
        outArgs = getCache(movieArgs{:});
        if ~isempty(outArgs)
            [movie,metadata] = outArgs{1:2};
            if isfield(pluginOptions,'pixelsize') && pluginOptions.pixelsize>0
                metadata.pixelsize = pluginOptions.pixelsize*metadata.head.ImgAcc_Pixelbinning;
                metadata.pixelsize_unit = [char(181) 'm'];% um
            end
            return;
        end
        
        if nargin<5 || isempty(fastLT) || ~logical(fastLT)
            [head,movie] = LinCAM_accumulate(filename_movie,{'head','tag'},[frame_binning, frame_range],pluginOptions.pixelbinning,timegate,window);
            movie = {permute(movie,[1 2 4 3])}; % shift empty channel dimension to the back
        else
            [head,movie,tau] = LinCAM_accumulate(filename_movie,{'head','tag','tau'},[frame_binning, frame_range],pluginOptions.pixelbinning,timegate,window);
            movie = {permute(movie,[1 2 4 3]),permute(tau,[1 2 4 3])};
        end
        
        if nargout>1 || (isfield(pluginOptions,'cacheMovie') && pluginOptions.cacheMovie>0)
            metadata = struct(...
                'filename',filename_movie,...
                'tau_unit','ns',...
                'head',head...
                );
            if isfield(pluginOptions,'pixelsize') && pluginOptions.pixelsize>0
                metadata.pixelsize = pluginOptions.pixelsize*metadata.head.ImgAcc_Pixelbinning;
                metadata.pixelsize_unit = [char(181) 'm'];% um
            end
        end
        % Save cache
        if isfield(pluginOptions,'cacheMovie')
            try
                setCache({movie,metadata},pluginOptions.cacheMovie, movieArgs{:});
            catch err
                warning('Could not cache the movie: Path might be read-only or contain special characters.');
                disp( getReport( err, 'extended', 'hyperlinks', 'on' ) );
            end
        end
    else
        error('No photons file.');
    end
end

function fhdl = getAccumFunc(pluginOptions)
    window = [pluginOptions.window_x, pluginOptions.window_x+pluginOptions.window_dx, pluginOptions.window_y, pluginOptions.window_y+pluginOptions.window_dy];
    fhdl = @(varargin)accumFunc(window,pluginOptions.pixelbinning,varargin{:});
end

function varargout = accumFunc(window,pixelbinning,varargin)
    varargin{6} = window;
    if isempty(varargin{4})% the fourth parameter could also be a mask. In that case the pixelbinning is derived from the mask size.
        varargin{4} = pixelbinning;
    end
    [varargout{1:nargout}] = LinCAM_accumulate(varargin{:});
end

function [tcspcdata,resolution] = getTCSPC(inputfile,maxPhotons)

    ms = photonscore.file_read(inputfile, '/photons/ms');
    
    count = min(ms(end),maxPhotons);
    offset = 0;    
    tcspcdata = photonscore.file_read(inputfile, '/photons/dt', offset, count);
    head = photonscore.file_info(inputfile);
    resolution = head.dt_channel*1e-12;  % in s
end

function options = findWindow(options)
    maxPhotons = 1e7;
    if isfield(options,'importFile') && exist(options.importFile,'file') 
        inputfile = options.importFile;
        getThres = @(vec,thres)quantile(vec,[0.1 0.99])*[1 -1;0 1]*[1-1/thres;1/thres];
        findFirstLast = @(vec)[find(vec,1) find(vec,1,'last')];
        
        ms = photonscore.file_read(inputfile, '/photons/ms');
        count = min(ms(end),maxPhotons);
        offset = 0;

        lineint = photonscore.file_read(inputfile, '/photons/x', offset, count);

        img_filt_sum = accumarray(lineint,1);
        if quantile(img_filt_sum,0.9)>getThres(img_filt_sum,100) % We do have dark padding -> zoom in
            xwind = (findFirstLast(movmean(img_filt_sum,5)>getThres(img_filt_sum,50)));
            
            options.window_x = xwind(1);
            options.window_dx = diff(xwind);
        end
        
        lineint = photonscore.file_read(inputfile, '/photons/y', offset, count);
        img_filt_sum = accumarray(lineint,1);
        if quantile(img_filt_sum,0.9)>getThres(img_filt_sum,100) % We do have dark padding -> zoom in
            ywind = (findFirstLast(movmean(img_filt_sum,5)>getThres(img_filt_sum,50)));
            
            options.window_y = ywind(1);
            options.window_dy = diff(ywind);
        end
        
    end
end

%% Cacheing

function argout = getCache(filename_movie, varargin)
    % varargin = {frame_range, frame_binning,fastLT,timegate,window}
    argout = {};
    try
        mf = openCache(filename_movie,false);
        if isempty(mf)
            return;
        end
        if any(strcmp(who(mf),'argout'))
            ncache = size(mf,'argout',1);
        else
            ncache = 0;
        end
        for icache = 1:ncache
            argin =  mf.argin(icache,1);
            argin = argin{1};
            if ~varargin{3}
                % Extra FLIM data does not interfere
                varargin{5} = argin{6};
                varargin{3} = argin{4};
            end
            if isequaln(argin,{filename_movie,varargin{:}}) %#ok<CCAT> This is not identical with cat, since cat removes empty cells.
                % All arguments are equal.
                argout = mf.argout(icache,1);
                argout = argout{1};
                mf.time_read(icache,1) = now; % Update last read timestamp
                break;
            end
        end
    catch err %#ok<NASGU>
        warning('Invalid cache file. Clearing cache.');
        openCache(filename_movie,false,true);  
    end
end
function setCache(argout,maxcache, filename_movie, varargin)
% varargin = {frame_range, frame_binning,fastLT,timegate,fastLTmethod};
    if maxcache == 0
        openCache(filename_movie,false,true);       
        return
    end
    mf = openCache(filename_movie,true);
    varlist = who(mf);
    if any(strcmp(varlist,'argout'))
        ncache = size(mf,'argout',1);
    else
        ncache = 0;
    end
    if ncache < maxcache 
        icache = ncache+1;
    else
        % find longest unread entry
        tread = mf.time_read;
        [~,icache] = min(tread);
    end
    mf.argin(icache,1) = {{filename_movie,varargin{:}}}; %#ok<CCAT> This is not identical with cat, since cat removes empty cells.
    mf.argout(icache,1) = {argout};
    mf.time_read(icache,1) = now;
    mf.time_created(icache,1) = now;
end
function mf = openCache(filename_movie,createFlag,clearFlag)
    cacheSuffix = '_FLIM.cache';
    [fpath,fname] = fileparts(filename_movie);
    cachename = fullfile(fpath,makeValidMatfileName([fname cacheSuffix]));% Matlab cannot write to matfiles with special characters in the name
    if exist(cachename,'file') && nargin > 2 && clearFlag
        delete(cachename);
    end
    if exist(cachename,'file') 
        % creates a matfile if it does not exist.
        mf = matfile(cachename,'Writable',true);        
    elseif (nargin>1 && createFlag)
        mf = matfile(cachename);
        mf.filename_movie = filename_movie;
    else
        mf = [];
    end
end