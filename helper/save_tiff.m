function save_tiff(filename, imgdata, exist_flag, imgresolution, tiff_tags)
% Usage: save_Tiff(filename, imgdata, datatype, bitdepth, force_overwrite)
%
% Saves an image or image stack to a Tiff file with specified datatype.
%
% imgdata: should be a matrix in the form [x y frame] or [x y frame RGB]. 
%
% exist_flag: false           - does nothing and generates a warning
%             1 ('overwrite') - overwrites the existing file
%             2 ('append')    - appends the frames
% imgresolution: sets the X and Y Pixelresolution in dots/cm. Can be scalar
%                or vector of [xRes yRes]
% tiff_tags: tiff_tags can be a struct or cell array {Name, Value} of tags
%            which will be embedded into the tiff. They must be part of the
%            TIFF specifications.
% 
% Note: Parameters can be left unspecified by giving an emtpy matrix -> []

% By
% Simon Christoph Stein - August 2014
% E-Mail: scstein@phys.uni-goettingen.de
%
% Modified and extend by CT, 2019
% - detects varable types (double/single, uint8/16)
% - supports appending or overwriting existing files
% - can add the img resolution
% - supports RGB tiffs
%

varname = inputname(2); % Get name of the variable to save

% -- Input parsing --
if(nargin<3) || isempty(exist_flag)
    exist_flag = false;
elseif ischar(exist_flag)
    exist_flag = find(strcmpi(exist_flag,{'overwrite','append'}));
end
if(nargin<4) || isempty(imgresolution)
    imgresolution = [];
end

% Check filename, append '.tif' if neccessary
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext, '.tif') == 0
    if isempty(pathstr)
        filename = [name '.tif'];
    else
        filename = [pathstr filesep name '.tif'];
    end
end

% Check if output file exists
if( exist(filename,'file') && ~exist_flag)
    warning(' File ''%s'' already exists. Set force_overwrite=true if you want to replace it. ', filename);
    return
end


% -- Setup Tiff metadata --
tagstruct.Compression = Tiff.Compression.None;
tagstruct.ImageLength = size(imgdata,1);
tagstruct.ImageWidth = size(imgdata,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
if ~isempty(imgresolution)
    tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Centimeter;
    tagstruct.XResolution = imgresolution(1);
    tagstruct.YResolution = imgresolution(end);
end

% Convert data to requested datatype
datatype = class(imgdata);
switch class(imgdata)
   case {'single','double'}
        imgdata = single(imgdata);
        tagstruct.BitsPerSample = 32;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
   case {'uint8','uint16','uint32','uint64'}
        tagstruct.BitsPerSample = str2double(datatype(5:end));
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt; 
   case {'int8','int16','int32','int64'}
        tagstruct.BitsPerSample = str2double(datatype(4:end));
        tagstruct.SampleFormat = Tiff.SampleFormat.Int; 
    otherwise
      error('Unsupported datatype ''%s\''. Supported: ''float'', ''uint'', ''int'' ', datatype);
end
% handle colour. Permute to [x y RGB frame]/[x y 1 frame]
imgdata = permute(imgdata,[1 2 4 3]);

tagstruct.SamplesPerPixel = size(imgdata,3);
% tagstruct.RowsPerStrip = 8; % Has to do with the memory layout of the data, better not touch..
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % Seperate saves R, G, B components seperatly instead of contigous (RGB RGB RGB)
tagstruct.Software = 'MATLAB';
% tagstruct.SubIFD = 1  % required to create subdirectories inside the TIFF file

%Add extra Tags
if nargin>4 && ~isempty(tiff_tags)
    if iscell(tiff_tags)
        tiff_tags = struct(tiff_tags{:});
    end
    for tagname = fieldnames(tiff_tags)'
        tagstruct.(tagname{:}) = tiff_tags.(tagname{:});
    end
end
% -- Write the movie --
fprintf('Settings:\n');
fprintf('\t\t\t   Datatype: %s\n', datatype);
disp(tagstruct) % Show the options
fprintf('Saving variable ''%s'' to file ''%s'' ..\n', varname, filename);


% Write frames
if exist_flag==2
    t = Tiff(filename,'a'); %Append frames
else
    t = Tiff(filename,'w');
end
cleanupTrigger = onCleanup(@() cleanupFunc(t));

msgAccumulator = ''; % needed for one-line printing
startTime = tic;
lastElapsedTime = 0;

for iF = 1:size(imgdata,4)
    elapsedTime = toc(startTime);    
    if( (elapsedTime-lastElapsedTime)>0.25) % Output every 0.25 seconds
        rewindMessages()
        rewPrintf('Saving Frame %i/%i ..',iF, size(imgdata,4));
        lastElapsedTime = elapsedTime;
    end
    
    t.setTag(tagstruct)
    t.write(imgdata(:,:,:,iF));        
    t.writeDirectory();
end
rewindMessages()
rewPrintf('Saving Frame %i/%i .. done\n',size(imgdata,4), size(imgdata,4));

t.close;


    function rewPrintf(msg, varargin)
        % Rewindable message printing: Print msg and cache it.
        % Usage is analogous to sprintf.
        msg = sprintf(msg, varargin{:});
        msgAccumulator = [msgAccumulator, msg];
        fprintf(msg);
    end

    function rewindMessages()
        % Remove cached messages from command line, reset cache
        reverseStr = repmat(sprintf('\b'), 1, length(msgAccumulator));
        fprintf(reverseStr);
        
        msgAccumulator = '';
    end

end



function cleanupFunc(t)
%     fprintf('\n')
    t.close;
    %     fprintf('\nUSER EXIT');
end
