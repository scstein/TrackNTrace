function [varargout] = PTU_accumulate(inputfile, outputs, framebinning, pixelbinning, timegate, BidirectShift, channelmap)
% function [varargout] = PTU_accumulate(inputfile, outputs, framebinning, pixelbinning, timegate,BidirectShift)
% PTU_accumulate Generates the requsted outputs by accumulating the photons from
% the inputfile using an intermediate index (see PTU_index).
%   
% inputfile    - PTU file to process 
%                OR matfile containing head, im_col, im_line and im_tcspc
% outputs      - cellarray with names of requested outputs. Options are:
%                 tag                   [x y chan frame]    Pixelwise intensity per frame
%                 tags                  [x y chan]          Pixelwise intensity over all frames
%                 tau                   [x y chan frame]    Fast lifetime (STD) per frame
%                 taus                  [x y chan]          Fast lifetime (STD) over all frames
%                 tcspc_pix             [x y bin chan]      Pixelwise TCSPC over all frames
%                 tcspc_pix_frame       [x y bin chan frame]Pixelwise TCSPC per frame (Warning: Matrix can become very large)
%                 head                  STRUCT              Head of the PTU file
%                 index                 STRUCT/matfile      Returns depending of the size the content of matfile or a handle of it.
%                 {'tau',afun,pfun}     [x y chan frame]    Like tau, but allows to specify a custom accumulation function afun and a postprocessing fun pfun which is applied to the matrix.
%                 {'arrival',afun,pfun} [x y chan frame]    Accumulation of the total arrival time (in units of the micro time) for dead time corretion or pixelwise timetraces.
%                                                           If afun accepts two inputs it is called like afun(head,mf) and should return a function handle for accumulating.
% framebinning - Gives the number of frames in the PTU which are summed to
%                one frame in the output. Default: 1
%                Can have the start and stop frame (unbinned, inclusive) as
%                second and third argument: [binning start stop] 
% pixelbinning - Gives the number of pixels in x and y which are summed to
%                one pixel in the output. Default: 1
%              - As alternative a mask with indices of molecules can be
%                given. In that case the binning is determined by the ratio
%                between the scan and mask. In the output x=1 and y=mol.
% timegate     - [min max] photons ariving not inbeween min and max are
%                rejected. Can be in ns or bin numbers (as int).
% BidirectShift- Correct the shift between foward and reverse lines for
%                bidirectional scans. Can be false/0 for no shift, a float
%                with the shift in pixels or true (logical) to determine
%                the shift automatically.
% channelmap   - A vector of logicals with true for the channels to include. 
%                Disabled channels are removed from the output. By default all 
%                non-empty channels are included.
%              - Alternatively, channels can be merged by providing a vector of
%                indeces. E.g. [1 0 2 2] excludes the second channel and merges
%                the third and fourth channel.
%
% (C) Christoph Thiele, 2020.

%%
% inputfile = 'W:\Christoph\190426_Silicarhodamine_beads\data.sptw\Vero-cells_Alexa647_11.ptu';
narginchk(2,7);
if nargin < 3 || isempty(framebinning)
    framebinning = 1;
end
if nargin < 4 || isempty(pixelbinning)
    pixelbinning = 1;
end
if nargin<5 || isempty(timegate) || (islogical(timegate)&&~timegate)
    timegate = [];
end
if nargin<6 || isempty(BidirectShift)
    BidirectShift = true;
end

% Determin the output indices and set the flags what to accumulate
varargout       = {};
if ~iscell(outputs)
    outputs         = cellstr(outputs);
end
outi_head       = find(strcmpi(outputs,'head'));
outi_tag        = find(strcmpi(outputs,'tag'));
outi_tags       = find(strcmpi(outputs,'tags'));
outi_tau        = find(strcmpi(outputs,'tau'));
outi_taus       = find(strcmpi(outputs,'taus'));
outi_tcspc_pix  = find(strcmpi(outputs,'tcspc_pix'));
outi_tcspc_pix_frame = find(strcmpi(outputs,'tcspc_pix_frame'));
outi_mf         = find(strcmpi(outputs,'index'));

outi_cell       = find(cellfun(@iscell,outputs));
outn_cell       = cellfun(@(c)c{1},outputs(cellfun(@iscell,outputs)),'UniformOutput',false);

accum_taus      = ~isempty(outi_taus);
accum_tcspc_f   = ~isempty(outi_tcspc_pix_frame);
accum_tcspc     = accum_taus || ~isempty(outi_tcspc_pix);
accum_tags      = accum_taus || ~isempty(outi_tags);
accum_tag       = accum_tags || ~isempty(outi_tag);
accum_tau       = ~isempty(outi_tau);
accum_tau_cust  = ~isempty(find(strcmpi(outn_cell,'tau'), 1));
accum_arrival   = ~isempty(find(strcmpi(outn_cell,'arrival'), 1));

if isstruct(inputfile)||strcmpi(class(inputfile),'matlab.io.MatFile')
    % When processing many diffrent options (e.g. frames or masks on the
    % same file, passing the content of the matfile reduces disk reads.
    % This bypasses the check whether the BidirectShift matches the
    % option and disables the waitbar.
    mf = inputfile;
    head = mf.head;
    batchFlag = true;
    
    % Add parameters to the head
    head.ImgAcc_Pixelbinning = pixelbinning;
    head.ImgAcc_Framebinning = framebinning;
    head.ImgAcc_Timegate     = timegate;
    
    if ~isempty(outi_head)
        varargout{outi_head} = head;
    end
else
    batchFlag = false;
    if endsWith(inputfile,'.mat')
        mfile = inputfile;
    else
        [inputpath,inputname] = fileparts(inputfile);
        mfile = fullfile(inputpath, makeValidMatfileName([inputname '_index.mat']));
    end
    
    if ~exist(mfile,'file')
        PTU_index(inputfile,BidirectShift);
        if ~exist(mfile,'file') % incase PTU_index fails
            error('Could not find PTU index file.');
        end
    end
    
    mf = matfile(mfile);
    head = mf.head;
    if head.ImgHdr_Dimensions == 1 % Error for point measurements.
        warning('Measurement is not a scan.');
        [varargout{1:nargout}] = deal([]);
        if ~isempty(outi_head)
            varargout{outi_head} = head;
        end
        return
    end
    % Recalculate index if the BidirectShift is not matching
    if endsWith(inputfile,'.ptu')&&head.ImgHdr_BiDirect&&((( isempty(BidirectShift) || (islogical(BidirectShift) && BidirectShift==false)) && isfield(head,'ImgHdr_BiDirectFittedOffset') && head.ImgHdr_BiDirectFittedOffset ~= 0)...
                                                          ||( ~isempty(BidirectShift) && ~(islogical(BidirectShift) && BidirectShift==false) && (~isfield(head,'ImgHdr_BiDirectFittedOffset') || isfloat(BidirectShift) && head.ImgHdr_BiDirectFittedOffset ~= BidirectShift)))
        PTU_index(inputfile,BidirectShift);
        mf = matfile(mfile);
        head = mf.head;
    end
    
    % Add parameters to the head
    head.ImgAcc_Pixelbinning = pixelbinning;
    head.ImgAcc_Framebinning = framebinning;
    head.ImgAcc_Timegate     = timegate;
    
    if ~isempty(outi_head)
        varargout{outi_head} = head;
    end
    if ~accum_taus && ~accum_tau && ~accum_tau_cust && ~accum_tags && ~accum_tag && ~accum_tcspc && ~accum_tcspc_f && ~accum_arrival
        % Return if only the head is requsted.
        return;
    end
    
    % If the variables take less than a quater of the available memory they
    % are loaded completly to reduce slow file reads.
    loadvars = {'im_frame','im_frame_index','im_chan','im_col','im_line','im_tcspc'}; % Do not load im_sync by default
    if accum_arrival
        loadvars{end+1} = 'im_sync';
    end
    if sum(feval(@(vars)[vars(cellfun(@(name)any(strcmp(name,loadvars)),{vars.name})).bytes],whos(mf)))<getFreeMem()/2
        mf = load(mfile,loadvars{:});
    end
end

lastframe = double(mf.im_frame(mysize(mf,'im_frame',1),mysize(mf,'im_frame',2)));
if numel(framebinning)>1
    firstframe = max(framebinning(2),1);
    if numel(framebinning)>2
        lastframe = min(lastframe,framebinning(3));
    end
else
    firstframe = 1;    
end
framebinning = max(min(framebinning(1),lastframe),1);

binfun = @(frame)ceil((frame-firstframe+1)/framebinning);
lastframe_binned = binfun(lastframe);

if ~batchFlag, h1 = waitbar(0,sprintf('Accumulating photons: Frame %d out of %d',0,lastframe_binned));end

if numel(pixelbinning)>1
    maskFlag = true;
    mask = pixelbinning;
    pixelbinning = [head.ImgHdr_PixY head.ImgHdr_PixX]./size(mask);
    if prod(pixelbinning)~=1
        warning('Masks do not work yet with pixelbinning. Check LinCAM_accumulate.m for a working implemenation.');
    end
    if diff(pixelbinning)~=0
        error('Mask not matching ptu header.');
    else
        pixelbinning = pixelbinning(1);
    end
    nx = max(mask(:));
    ny = 1;    
else
    maskFlag = false;
    if isinf(pixelbinning)
        pixelbinning = 1/eps; % use a very large number instead of inf to avoid zeros in the look up.
    end
    nx = ceil(head.ImgHdr_PixX/pixelbinning);
    ny = ceil(head.ImgHdr_PixY/pixelbinning);
end

Resolution = max(1e9*head.MeasDesc_Resolution);
% chDiv      = 1e-9*Resolution/head.MeasDesc_Resolution;
% SyncRate   = 1./head.MeasDesc_GlobalResolution;
Ngate   = ceil(1e9*head.MeasDesc_GlobalResolution./Resolution);

%
dind = double(unique(mf.im_chan(1:min(1e6,mysize(mf,'im_chan',1)),1)));
if nargin < 7 || isempty(channelmap)
    if isfield(head,'MeasDesc_Nchan')
        maxch_n = head.MeasDesc_Nchan;
        if maxch_n ~= numel(dind)
            warning('Channel number in header does not match with channel indecs.');
        end
    else
        maxch_n = numel(dind);
        head.MeasDesc_Nchan = maxch_n;
    end
    channelmap = zeros(max(dind),1);
    channelmap(dind) = 1:numel(dind);
else % determine number of output channels from channelmap argument
    if islogical(channelmap)
        % one channel for each true channelmap 
        maxch_n = sum(channelmap);
        channelmap = find(channelmap);
    else
        % number of channels is maxium index in channelmap
        maxch_n = max(channelmap);
    end
    % extend to number of channels with zeros
    channelmap(end+1:max(dind)) = 0;
end

if nargin>4 && ~isempty(timegate)
    % Interprete unit of timegate as ns if it is below the length of the
    % sync periode or a non integer.
    % Interpretation as bin indices can be enforced by the use of a
    % interger data type.
    if isfloat(timegate) && max(timegate)<Resolution*Ngate || max(mod(timegate,1))> eps
        timegate = timegate./Resolution;
    end
    
    timegate(1) = max(timegate(1),1);                   % The timegate is inclusive. Set it to 1 to make calculation of Ngate easier.
    timegate(2) = min([timegate(2:end), double(Ngate)]);% Ensure that second entry exists
    Ngate = min(Ngate,diff(timegate)+1);                % Is automatically casted to class of Ngate
else
    timegate = [];
end

% Set up custom accumumaltion functions
if ~isempty(outi_cell)
    % For now only accept known types
    cell_types = {'tau','arrival'};
    outi_cell_typ = cellfun(@(c)find(strcmpi(c,cell_types),1),outn_cell,'UniformOutput',false);
    ind = ~cellfun(@isempty,outi_cell_typ);
    if ~all(ind)
        warning('Some requested outputs are not recognised and ignored: %s\n',outn_cell{find(~ind,1)});
    end
    outi_cell_typ = cell2mat(outi_cell_typ(ind));
    outi_cell = outi_cell(ind);
    outn_cell = cell_types(outi_cell_typ);
    
    outv_cell = cell(size(outi_cell(:)));
    afun_cell = cell(size(outi_cell)); % empty defaults to sum
    pfun_cell = cell(size(outi_cell));
    for cidx = 1:numel(outi_cell)
        if numel(outputs{outi_cell(cidx)})>1 && ~isempty(outputs{outi_cell(cidx)}{2})
            afun_cell{cidx} = outputs{outi_cell(cidx)}{2};
            if nargin(afun_cell{cidx})==2
               afun_cell{cidx} = afun_cell{cidx}(head,mf);
            end
        end
        if numel(outputs{outi_cell(cidx)})>2 && ~isempty(outputs{outi_cell(cidx)}{3})
            pfun_cell{cidx} = outputs{outi_cell(cidx)}{3};
        elseif iscell(afun_cell{cidx}(1)) % Do not postprocess cells by default
            pfun_cell{cidx} = @(x)x;
        else % Convert bin to time
            pfun_cell{cidx} = @(x)Resolution.*x;
        end
    end
end

if accum_tcspc
    tcspc_pix = zeros(ny,nx,Ngate,maxch_n,intminclass(head.TTResult_NumberOfRecords/max(1,nx*ny))); %This class should be quite safe. The estimation is that all photons are in one TCSPC channel, but evenly distributed over the pixels
end
if accum_tcspc_f
    tcspc_pix_frame = zeros(ny,nx,Ngate,maxch_n,lastframe_binned,intminclass(head.TTResult_NumberOfRecords/max(1,nx*ny)));
end
if accum_tag
    tag = zeros(ny,nx,maxch_n,lastframe_binned);
end
if accum_tag
    tau = nan(ny,nx,maxch_n,lastframe_binned);
end

if accum_tcspc || accum_tcspc_f
    sub_class = intminclass([ny nx maxch_n lastframe_binned Ngate]);  % Usually 'uint16'
else
    sub_class = intminclass([ny nx maxch_n lastframe_binned]);  % Usually 'uint16'
end

frame_lookup = cast(binfun(1:lastframe),sub_class)';
frame_fun = @(frames)frame_lookup(frames)-frame_lookup(frames(1,:))+1; % Relative indices of the frames. The first frame in frames allways gets the index 1.
px_lookup = cast(ceil((1:head.ImgHdr_PixX)/pixelbinning),sub_class)';
py_lookup = cast(ceil((1:head.ImgHdr_PixY)/pixelbinning),sub_class)';

if maskFlag
    pyx_lookup = @(y,x)[ones(size(y)) mask(sub2ind(size(mask),y,x))];
else
    pyx_lookup = @(y,x)[py_lookup(y) px_lookup(x)];
end

c_lookup = cast(channelmap(:),sub_class);

% Add the end of the last frame, select the binned frame_indeces and apply 
% the selection of first and last frame.
im_frame_index = [mf.im_frame_index;(mysize(mf,'im_frame',1)+1)];
im_frame_index = im_frame_index([firstframe:framebinning:lastframe, lastframe+1]);

% Set the number to process each iteration. If the total number of photons
% is below 1e8 (equal max_photons for 4 GB) all photons are processed at
% once. This avoids the call to getFreeMem() which takes approx 0.1 s.
if im_frame_index(end)-im_frame_index(1)>1e8 
    max_photons = getFreeMem()/(4*intbyte(sub_class)+32);       % Size of subs + tcspcdata + 3 extra double for intermediate vars.
    % max_photons = 1e6;
else
    max_photons = inf;
end

cframe = 1;
tcspcdata = [];
n_it = 0; % Iteration couter.
while cframe<=lastframe_binned
    n_it = n_it +1;
    clastframe = max(cframe+1,sum(im_frame_index < max_photons+im_frame_index(cframe))); % Current last frame: take as many frames as possible while staying below maxphotons. Take at least one full frame.
                                                                                         % clastframe is the first frame of the NEXT iteration.
%     clastframe = cframe+1;
    ind = im_frame_index(cframe,1):(im_frame_index(clastframe,1)-1);
    if isempty(ind)
        % Can occure when frames without photons are present
        cframe = clastframe;
        continue;
    end
    
    subs = [
            pyx_lookup(mf.im_line(ind,1),...
                       mf.im_col(ind,1))...
            c_lookup(mf.im_chan(ind,1)),...
            frame_fun(mf.im_frame(ind,1))];
        
    if accum_arrival
        syncs = mf.im_sync(ind,1);
    end
    if ~isempty(timegate)
        tcspcdata = mf.im_tcspc(ind,1);
        ind = tcspcdata>=timegate(1) & tcspcdata<=timegate(2);
        subs = subs(ind,:);
        tcspcdata = double(tcspcdata(ind,:)-timegate(1)+1);
    elseif accum_tau || accum_tau_cust || accum_tcspc || accum_tcspc_f || accum_arrival
        tcspcdata = double(mf.im_tcspc(ind,1));
    end
    ind = all(subs,2) & subs(:,3) <= maxch_n;
    subs = subs(ind,:); % Remove zero indices (eg. out of frame, out of pixel mask)
    if ~isempty(tcspcdata)
        tcspcdata = tcspcdata(ind);
    end        
    if accum_arrival
        syncs = syncs(ind) - syncs(find(ind,1)); % This avoids that the precision of double eps(syncs(end)) decreases to more than 1 for long measurements. Might interfer if absolute timeing is critical.
    end
    clear ind;
    
    if accum_tag
        tag(:,:,:,cframe:clastframe-1) = accumarray(subs,1,[ny nx maxch_n clastframe-cframe]);
    end
    if accum_tau
        tau(:,:,:,cframe:clastframe-1) = Resolution*sqrt(accumarray(subs,tcspcdata,[ny nx maxch_n clastframe-cframe],@myvar));
    end
    if  accum_tcspc
        tcspc_pix = tcspc_pix + cast(accumarray([subs(:,[1 2]) tcspcdata subs(:,3)],1,[ny,nx,Ngate,maxch_n]),class(tcspc_pix));
    end
    if  accum_tcspc_f
        tcspc_pix_frame(:,:,:,:,cframe:clastframe-1) = cast(accumarray([subs(:,[1 2]) tcspcdata subs(:,[3 4])],1,[ny,nx,Ngate,maxch_n,lastframe_binned]),class(tcspc_pix_frame));
    end
    if ~isempty(outi_cell)
        ind = strcmpi(outn_cell,'tau');
        outv_cell(ind,n_it) = cellfun(@(a,p)p(accumarray(subs,tcspcdata,[ny nx maxch_n clastframe-cframe],a)),afun_cell(ind),pfun_cell(ind),'UniformOutput',false); %Number of iterations is unknown at begining. Storing in a growing cell array and concating at the end causes the least memomory load, since a cell array is a vector of pointers.
        
        if accum_arrival
            ind = strcmpi(outn_cell,'arrival');
            syncs = tcspcdata + syncs .* ceil(1e9*head.MeasDesc_GlobalResolution./Resolution); % Not Ngate, because Ngate gets modified by timegateing
            outv_cell(ind,n_it) = cellfun(@(a,p)p(accumarray(subs,syncs,[ny nx maxch_n clastframe-cframe],a)),afun_cell(ind),pfun_cell(ind),'UniformOutput',false); %Number of iterations is unknown at begining. Storing in a growing cell array and concating at the end causes the least memomory load, since a cell array is a vector of pointers.
        end
    end
    
    if ~batchFlag
        if isgraphics(h1)
            waitbar(cframe/lastframe_binned,h1,sprintf('Calculating fast lifetime: Frame %d out of %d',cframe,lastframe_binned));
        else % Waitbar was closed, abort
            return;            
        end
    end
    cframe = clastframe;
end

if ~batchFlag, waitbar(1,h1,'Calculating merged variables...');end

if ~isempty(outi_cell)
    outv_cell = arrayfun(@(c)cat(4,outv_cell{c,:}),1:size(outv_cell),'UniformOutput',false);
end

if accum_tags
    tags = sum(tag,4);
end
if accum_taus
    bin = shiftdim((1:Ngate)',-2)*Resolution; % 1D time axis, but in the 3rd dimension. Much smaller.
    % taus = real(sqrt((squeeze(sum(bin.^2.*double(tcspc_pix),3))./tags)-(squeeze(sum(bin.*double(tcspc_pix),3))./tags).^2));
    taus = real(sqrt(((squeeze(sum(bin.^2.*double(tcspc_pix),3)))-((squeeze(sum(bin.*double(tcspc_pix),3))).^2)./tags)./(tags-1))); % Normalise with (N-1) to make it consitent with (my)var
    taus(tags==1) = 0;
    taus(tags==0) = NaN;
end

if ~batchFlag, close(h1);end
% 
% to enable support for multiple times the same output use deal:
% [varargout{outi_tag}] = deal(tag);
if ~isempty(outi_tag)
    varargout{outi_tag} = tag;
end
if ~isempty(outi_tags)
    varargout{outi_tags} = tags;
end
if ~isempty(outi_tau)
    varargout{outi_tau} = tau;
end
if ~isempty(outi_taus)
    varargout{outi_taus} = taus;
end
if ~isempty(outi_tcspc_pix)
    varargout{outi_tcspc_pix} = tcspc_pix;
end
if ~isempty(outi_tcspc_pix_frame)
    varargout{outi_tcspc_pix_frame} = tcspc_pix_frame;
end
if ~isempty(outi_mf)
    if isstruct(mf) % Do not modify the head in case mf is a matfile.
        mf.head = head;
    end
    varargout{outi_mf} = mf;
end
if ~isempty(outi_cell)
    varargout(outi_cell) = outv_cell;
end
end


function freeMem = getFreeMem()
try % Works only on Windows
    [~,freeMem] =  memory();
    freeMem = freeMem.PhysicalMemory.Available;
catch
    freeMem =  4e9;% Use 4 GB as default.
end
end

function b = intbyte(x,arraySizeFlag)
% function b = intbyte(x) 
% Returns the size in byte of the class x. If x is an integer the class is
% determined by the class of x.

narginchk(1,2);

% Multiplier if arraySizeFlag is set.
if nargin>1 && arraySizeFlag
    m = numel(x);
else
    m = 1;
end

if isnumeric(x)
    x = class(x);
end

switch x
    case {'int8','uint8'}
        b=1;
    case {'int16','uint16'}
        b=2;
    case {'int32','uint32'}
        b=4;
    case {'int64','uint64'}
        b=8;
    case 'double'
        b=8;
    case 'single'
        b=4;
    otherwise
        error('Unknown class %s',class);
end

b=b*m;
end

function y = myvar(x)
% Simplified version of var for the use with accumarry. Only works on the
% first dimension. Does behave diffrent in some edge cases like inputs
% containing inf or NaN which should not occure in this use case.
    dim = 1;
    n = size(x,dim);
    if n > 1
        y = sum((x - sum(x,dim)./n).^2, dim) ./ (n-1);
    elseif n==1
        y = 0;
    else
        y = NaN;
    end
end

function y = mysize(x,name,dim)
% variant of size which allows equal inputs for matfiles and structs
    if nargin<2
        name = '';
    end
    if strcmpi(class(x),'matlab.io.MatFile')
        if nargin<3 || isempty(dim)
            y = size(x,name);
        else
            y = size(x,name,dim); 
        end
    else
        if ~isempty(name)
            x = x.(name);
        end
        if nargin<3 || isempty(dim)
            y = size(x);
        else
            y = size(x,dim);
        end
    end    
end
