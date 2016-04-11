function save_tiff(filename, imgdata, datatype, bitdepth, force_overwrite)
% Usage: save_Tiff(filename, imgdata, datatype, bitdepth, force_overwrite)
%
% Saves an image or image stack to a Tiff file with specified datatype.
%
% datatype: 'float', 'uint', 'int' | default: 'float'
% bitdepth: 8/16/32/64 | default: 16 for int/uint (float is always 32 bit)
% forceoverwrite: overwrite output file if it already exists | default: false
% 
% Note: Parameters can be left unspecified by giving an emtpy matrix -> []

% By
% Simon Christoph Stein - August 2014
% E-Mail: scstein@phys.uni-goettingen.de
%

varname = inputname(2); % Get name of the variable to save

% -- Input parsing --
if(nargin < 3) || isempty(datatype)
    datatype = 'float';
end

if(nargin<4) || isempty(bitdepth)
    bitdepth = 16;
end

if(nargin<5) || isempty(force_overwrite)
    force_overwrite = false;
end

% Check filename, append '.tif' if neccessary
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext, '.tif') == 0
    filename = [name '.tif'];
end

% Check if output file exists
if( exist(filename,'file') && ~force_overwrite)
    error(' File ''%s'' already exists. Set force_overwrite=true if you want to replace it. ', filename);
end


% -- Setup Tiff metadata --
tagstruct.Compression = Tiff.Compression.None;
tagstruct.ImageLength = size(imgdata,1);
tagstruct.ImageWidth = size(imgdata,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;

% Convert data to requested datatype
switch datatype
   case 'float'
        imgdata = single(imgdata);
        tagstruct.BitsPerSample = 32;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
   case 'uint'
        switch bitdepth
            case 8
                 imgdata = uint8(imgdata);
            case 16            
                 imgdata = uint16(imgdata);
            case 32            
                 imgdata = uint32(imgdata);
            case 64
                 imgdata = uint64(imgdata);
            otherwise
                error('Unsupported bitdepth. Supported are 8/16/32/64')
        end
        tagstruct.BitsPerSample = bitdepth;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt; 
   case 'int'
        switch bitdepth
            case 8
                 imgdata = int8(imgdata);
            case 16            
                 imgdata = int16(imgdata);
            case 32            
                 imgdata = int32(imgdata);
            case 64
                 imgdata = int64(imgdata);
            otherwise
                error('Unsupported bitdepth. Supported are 8/16/32/64')
        end
        tagstruct.BitsPerSample = bitdepth;
        tagstruct.SampleFormat = Tiff.SampleFormat.Int; 
    otherwise
      error('Unsupported datatype ''%s\''. Supported: ''float'', ''uint'', ''int'' ', datatype);
end

tagstruct.SamplesPerPixel = 1;
% tagstruct.RowsPerStrip = 8; % Has to do with the memory layout of the data, better not touch..
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % Seperate saves R, G, B components seperatly instead of contigous (RGB RGB RGB)
tagstruct.Software = 'MATLAB';
% tagstruct.SubIFD = 1  % required to create subdirectories inside the TIFF file



% -- Write the movie --
fprintf('Settings:\n');
fprintf('\t\t\t Datatype: %s\n', datatype);
disp(tagstruct) % Show the options
fprintf('Saving variable ''%s'' to file ''%s'' ..\n', varname, filename);


% Write frames
t = Tiff(filename,'w');
cleanupTrigger = onCleanup(@() cleanupFunc(t));

msgAccumulator = ''; % needed for one-line printing
startTime = tic;
lastElapsedTime = 0;

for iF = 1:size(imgdata,3)
    elapsedTime = toc(startTime);    
    if( (elapsedTime-lastElapsedTime)>0.25) % Output every 0.25 seconds
        rewindMessages()
        rewPrintf('Saving Frame %i/%i ..',iF, size(imgdata,3));
        lastElapsedTime = elapsedTime;
    end
    
    t.setTag(tagstruct)
    t.write(imgdata(:,:,iF));        
    t.writeDirectory();
end
rewindMessages()
rewPrintf('Saving Frame %i/%i .. done\n',size(imgdata,3), size(imgdata,3));

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
        fprintf([reverseStr]);
        
        msgAccumulator = '';
    end

end



function cleanupFunc(t)
%     fprintf('\n')
    t.close;
    %     fprintf('\nUSER EXIT');
end
