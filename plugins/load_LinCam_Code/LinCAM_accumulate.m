function [varargout] = LinCAM_accumulate(inputfile, outputs, framebinning, pixelbinning, timegate, window)
%PTU_exportTIFF Converts a multiframe PTU file in a multiframe tiff
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
%                                                           If afun accepts two inputs it is called like afun(head,Nagte) and should return a function handle for accumulating.
%                Even though the LinCam has no channels the dimesion chan
%                is kept for compatibility.
% framebinning - Gives the time in ms which are summed to
%                one frame in the output. Default: 1000
%                Can have the start and stop frame (unbinned, inclusive) as
%                second and third argument: [binning start stop] 
% pixelbinning - Gives the number of pixels in x and y which are summed to
%                one pixel in the output. Default: 1
%              - As alternative a mask with indices of molecules can be
%                given. In that case the binning is determined by the ratio
%                between the scan and mask. In the output x=1 and y=mol.
% timegate     - [min max] photons ariving not inbeween min and max are
%                rejected. Can be in ns or bin numbers (as int).
%
% Example: Correlate the whole file:
%  LinCAM_accumulate(filename,{{'arrival',@(head,Ngate)feval(@(Nsub,Timeunit)@(t){{Timeunit*cumsum(reshape(repmat(2.^(0:ceil(log2(maxtime(1)/Timeunit/Nsub))-1),Nsub,1),Ncasc*Nsub,1)),tttr2xfcs(t,true(size(t)),ceil(log2(maxtime(1)/Timeunit/Nsub)),Nsub)}},10,head.MeasDesc_GlobalResolution)}},inf,inf,false)
%
%%
% inputfile = 'W:\Christoph\190426_Silicarhodamine_beads\data.sptw\Vero-cells_Alexa647_11.ptu';
narginchk(2,6);
if nargin < 3 || isempty(framebinning)
    framebinning = 1000;
end
if nargin < 4 || isempty(pixelbinning)
    pixelbinning = 1;
end
if nargin<5 || isempty(timegate)
    timegate = true;
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
    
    readOffset = 0;
    
    % Add parameters to the head
    head.ImgAcc_Pixelbinning = pixelbinning;
    head.ImgAcc_Framebinning = framebinning;
    head.ImgAcc_Timegate     = timegate;
else
    batchFlag = false;
    if endsWith(inputfile,'.mat')
        mf = load(inputfile);
        head = mf.head;
    else
        
        head = photonscore.file_info(inputfile);
        
        % Add parameters to the head
        head.ImgAcc_Pixelbinning = pixelbinning;
        head.ImgAcc_Framebinning = framebinning;
        head.ImgAcc_Timegate     = timegate;
    end
    
    if ~isempty(outi_head)
        varargout{outi_head} = head;
    end
    if ~accum_taus && ~accum_tau && ~accum_tau_cust && ~accum_tags && ~accum_tag && ~accum_tcspc && ~accum_tcspc_f && ~accum_arrival
        % Return if only the head is requsted. The head misses some fields
        % which are added based on the read photons.
        return;
    end
    if isempty(outi_mf) && numel(framebinning)>1
        % Only read the relvant part if min/max frame is given and no
        % output as struct is requested
        range_seconds = [max(framebinning(2)-1,0), min([inf, framebinning(3)])]*1e-3; % range in s
        mf = photonscore.read_photons(inputfile, range_seconds);
        readOffset = mf.ms(round(range_seconds(1)*1e3)+1); % round is necessary due to numerical precission
    else    
        mf = photonscore.read_photons(inputfile);
        readOffset = 0;
    end
    head.ImgHdr_PixX = 4096;%double(max(mf.x));
    head.ImgHdr_PixY = 4096;%double(max(mf.y));
    head.MeasDesc_Resolution = head.dt_channel*1e-12;
    head.MeasDesc_Ngate = max(mf.dt);
    head.MeasDesc_Nchan = 1;
    head.MeasDesc_GlobalResolution = 1/head.aat_frequencty;
    mf.head = head;
end

lastframe = numel(mf.ms)-1;
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

if ~batchFlag, h1 = waitbar(0,sprintf('Calculating fast lifetime: Frame %d out of %d',0,lastframe_binned));end

nx = head.ImgHdr_PixX;
ny = head.ImgHdr_PixY;

xrange = [1; inf];
yrange = [1; inf];
if nargin>5 && ~isempty(window)
    if islogical(window) && window
        window = LinCAM_findWindow(inputfile);
    end
    if numel(window)==1 % Also gets called for window==false
        xrange = [window; inf];
    elseif numel(window)>1
        xrange = window(1:2);
        if numel(window)==3
            yrange = [window(3); inf];
        else
            yrange = window(3:4);
        end
    end
end
xrange = min(max(xrange,1),nx);
yrange = min(max(yrange,1),ny);

if numel(pixelbinning)>1
    maskFlag = true;
    mask = pixelbinning;
    pixelbinning = [diff(yrange) diff(xrange)]./size(mask);
    if abs(diff(pixelbinning))>0.5
        error('Binning apears to be non isotropic.');
    else
        pixelbinning = round(pixelbinning(1));
    end
else
    maskFlag = false;
    if isinf(pixelbinning)
        pixelbinning = 1/eps; % use a very large number instead of inf to avoid zeros in the look up.
    end
end

px_lookup = zeros(nx,1);
px_lookup(xrange(1):xrange(2))=ceil((1:diff(xrange)+1)/pixelbinning)';
py_lookup = zeros(ny,1);
py_lookup(yrange(1):yrange(2))=ceil((1:diff(yrange)+1)/pixelbinning)';

if maskFlag
    nx = max(mask(:));
    ny = 1;    
else
    nx = px_lookup(xrange(2));
    ny = py_lookup(yrange(2));
end

if maskFlag
    % pad mask with zeros and modify the lookups to point out of ROI
    % photons to zero.
    mask = padarray(mask,[1 1],0);
    px_lookup = px_lookup + 1;
    py_lookup = py_lookup + 1;
    pyx_lookup = @(y,x)[ones(size(y)) mask(sub2ind(size(mask),py_lookup(y),px_lookup(x)))];
else
    pyx_lookup = @(y,x)[py_lookup(y) px_lookup(x)];
end



Resolution = 1e9*head.MeasDesc_Resolution;
Ngate   = head.MeasDesc_Ngate;
maxch_n = head.MeasDesc_Nchan; % Only one channel

if (islogical(timegate) && ~timegate) || isempty(mf.dt)
    % timegate === false OR no photons
    timegate = [];
elseif ~isempty(timegate) && ~(islogical(timegate) && ~timegate)
    if islogical(timegate) && timegate
        timegate = [min(mf.dt) max(mf.dt)];
    else
        % Interprete unit of timegate as ns if it is below the length of the
        % sync periode or a non integer.
        % Interpretation as bin indices can be enforced by the use of a
        % interger data type.
        if isfloat(timegate) && max(timegate)<Resolution*Ngate || max(mod(timegate,1))> eps
            timegate = timegate./Resolution;
        end
        timegate(end+1:2) = inf; % Ensure that second entry exists
    end
    timegate(2) = min(timegate(2),double(Ngate));
    Ngate = min(Ngate,diff(timegate)+1); %Is automatically casted to class of Ngate
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
    tcspc_pix_class = intminclass(numel(mf.x)/nx/ny);
    tcspc_pix = zeros(ny,nx,Ngate,tcspc_pix_class); %This class should be quite safe. The estimation is that all photons are in one TCSPC channel, but evenly distributed over the pixels
end
if accum_tcspc_f
    tcspc_pix_frame = zeros(ny,nx,Ngate,1,lastframe_binned,intminclass(numel(mf.x)/nx/ny));
end
if accum_tag
    tag = zeros(ny,nx,1,lastframe_binned);
end
if accum_tag
    tau = nan(ny,nx,1,lastframe_binned);
end

if ~accum_arrival && isempty(outi_mf)
    mf = rmfield(mf,'t');
elseif accum_arrival && ~isfield(mf,'t')
    error('Macrotime not available.');
end

% if accum_tcspc || accum_tcspc_f
%     sub_class = intminclass([ny nx lastframe_binned Ngate]);  % Usually 'uint16'
% else
%     sub_class = intminclass([ny nx lastframe_binned]);  % Usually 'uint16'
% end
sub_class = 'double';
if lastframe_binned>1
%     max_photons = getFreeMem()/(4*intbyte(sub_class)+32);       % Size of subs + tcspcdata + 3 extra double for intermediate vars.
    max_photons = getFreeMem()/(4*intbyte(sub_class)+40);       % Size of subs + tcspcdata + 4 extra double for intermediate vars.
else
    % For now we cannot split frames anyway, so this saves some time
    max_photons = inf;
end
% max_photons = 1e6;
% frame_lookup = cast(binfun(1:lastframe),sub_class)';
% frame_fun = @(frames)frame_lookup(frames)-frame_lookup(frames(1,:))+1; % Relative indices of the frames. The first frame in frames allways gets the index 1.


% Add the end of the last frame, select the binned frame_indeces and apply 
% the selection of first and last frame.
im_frame_index = mf.ms([firstframe:framebinning:lastframe, lastframe+1])+1-readOffset;

% Generate im_frame to allow multiframe processing
% im_frame = zeros(size(mf.t),intminclass(numel(im_frame_index)));
% im_frame(im_frame_index(1:end-1))=1;
% im_frame = cumsum(im_frame);

%cframe = 1;     % current start frame (inclusive)
clastframe = 1; % current end frame (exclusive) -> next start frame
tcspcdata = [];
n_it = 0; % Iteration couter.
while clastframe<=lastframe_binned
    n_it = n_it +1;
    cframe = clastframe;
    clastframe = max(cframe+1,sum(im_frame_index < max_photons+im_frame_index(cframe))); % Current last frame: take as many frames as possible while staying below maxphotons. Take at least one full frame.
    ind = im_frame_index(cframe,1):(im_frame_index(clastframe,1)-1);
    if isempty(ind)
        continue;
    end
    im_frame = zeros(numel(ind),1);
    im_frame(im_frame_index(cframe:(clastframe-1),1)-im_frame_index(cframe,1)+1)=1;
    im_frame = cumsum(im_frame);
        
    subs = [
            pyx_lookup(mf.y(ind,1),...
                       mf.x(ind,1)),...
            im_frame,...
           ];
        
    if accum_arrival
        macrotime = mf.t(ind,1);
    end        
    if ~isempty(timegate)
        tcspcdata = mf.dt(ind,1);
        ind = tcspcdata>=timegate(1) & tcspcdata<=timegate(2);
        subs = subs(ind,:);
        tcspcdata = double(tcspcdata(ind,:)-timegate(1)+1);
    elseif accum_tau || accum_tau_cust || accum_tcspc || accum_tcspc_f || accum_arrival
        tcspcdata = double(mf.dt(ind,1));
    end
    ind = all(subs,2);
    subs = subs(ind,:); % Remove zero indices (eg. out of frame, out of pixel mask)
    
    if isempty(subs)
        continue;
    end    
    if ~isempty(tcspcdata)
        tcspcdata = tcspcdata(ind);
    end        
    if accum_arrival
        macrotime = macrotime(ind) - macrotime(find(ind,1)); % This avoids that the precision of double eps(syncs(end)) decreases to more than 1 for long measurements. Might interfer if absolute timeing is critical.
    end
    clear ind;
    
    if accum_tag
        tag(:,:,1,cframe:clastframe-1) = accumarray(subs,1,[ny nx (clastframe-cframe)]);
    end
    if accum_tau
        tau(:,:,1,cframe:clastframe-1) = sqrt(accumarray(subs,tcspcdata*Resolution,[ny nx (clastframe-cframe)],@myvar));
    end
    if  accum_tcspc
        tcspc_pix = tcspc_pix + cast(accumarray([subs(:,[1 2]) tcspcdata],1,[ny,nx,Ngate]),tcspc_pix_class);
    end
    if  accum_tcspc_f
        tcspc_pix_frame(:,:,:,1,cframe:clastframe-1) = cast(accumarray([subs(:,[1 2]) tcspcdata subs(:,3)],1,[ny,nx,Ngate,(clastframe-cframe)]),class(tcspc_pix_frame));
    end
    if ~isempty(outi_cell)
        ind = strcmpi(outn_cell,'tau');
        outv_cell(ind,n_it) = cellfun(@(a,p)p(permute(accumarray(subs,tcspcdata,[ny nx clastframe-cframe],a),[1 2 4 3])),afun_cell(ind),pfun_cell(ind),'UniformOutput',false); %Number of iterations is unknown at begining. Storing in a growing cell array and concating at the end causes the least memomory load, since a cell array is a vector of pointers.
                                                                                                                                                                               %Permuting to introduce singleton dimension for channels
        
        if accum_arrival
            ind = strcmpi(outn_cell,'arrival');
%             macrotime = tcspcdata + macrotime .* ceil(1e9*head.MeasDesc_GlobalResolution./Resolution); % Not Ngate, because Ngate gets modified by timegateing
            macrotime = double(macrotime); % For now just the macro time. THIS IS NOT CONSTENT WITH PTU_accumulate! The macrotime has the head.MeasDesc_GlobalResolution (usually 10 ns).
            outv_cell(ind,n_it) = cellfun(@(a,p)p(permute(accumarray(subs,macrotime,[ny nx clastframe-cframe],a),[1 2 4 3])),afun_cell(ind),pfun_cell(ind),'UniformOutput',false); %Number of iterations is unknown at begining. Storing in a growing cell array and concating at the end causes the least memomory load, since a cell array is a vector of pointers.
        end
    end
    
    if ~batchFlag
        if isgraphics(h1)
            waitbar((clastframe-1)/lastframe_binned,h1,sprintf('Calculating fast lifetime: Frame %d out of %d',(clastframe-1),lastframe_binned));
        else % Waitbar was closed, abort
            return;            
        end
    end
end

if ~batchFlag, waitbar(1,h1,'Calculating overall fast lifetime...');end

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

% Write outputs at the requsted positions
if ~isempty(outi_head) % Update head
    varargout{outi_head} = head;
end
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
    mf.head = head;
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

function out = paren(x, varargin)
    out = x(varargin{:});    
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
            y = size(x,name,dim); %#ok<GTARG>
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
