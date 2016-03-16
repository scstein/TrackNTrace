function imgdata = read_BigTiff(filename, convert2double, minmax_frame)
% Usage: imgdata = read_tiff(filename, convert2double, minmax_frame)
%
% Loads an image or image stack from a Tiff file. It possible to read only
% a specified part of the file. File can be read in their original format
% (instead of converting to double).
%
% This is used by the read_tiff.m function for reading files larger then 4GB. 
% Implementation is somewhat shaky, reads ONLY 16-bit uint format at this point
% and will probabably not work for all Tiff files, but many imageJ files should work.
%
% Input:
%   filename: Tiff file to read
%   convert2double: Datatype conversion to double after reading | default: true
%   minmax_frame: [minframe,maxframe] to read. If only one number is specified,
%                 it is assumed to be maxframe and [1, maxframe] is read.
%                 If maxframe==inf the file is read from minframe to the end.| default: read whole movie
%
% NOTE: Deactivate double conversion if not neccessary, reading the same datatype used in the Tiff file.
%       This potentially saves a lot of memory (factor 4 for 16-bit Tiff movies) and thus increased performance.
%       You can still perform a cast to double (or other datatypes) later if desired.
%         Example: imgdata = double( read_tiff(filename) );

% By
% Simon Christoph Stein - August 2015
% E-Mail: scstein@phys.uni-goettingen.de
%

fprintf('Reading file ''%s''.. \n', filename)

if nargin < 2 || isempty(convert2double)
    convert2double = true;
end

if nargin < 3
    minmax_frame = [];
end

% -- Create Tiff object --
% Check filename, append '.tif' if neccessary
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext, '.tif') == 0
    filename = [filename '.tif'];
end

iminfo = imfinfo(filename);
if(numel(iminfo)>1)
   error('Expected only one tiff directory description for files > 4GB.')  
end
if( ~strcmp(iminfo.Format,'tif'))
   error('File is not a tiff file.') 
end
nrFramesOfFile = floor( (iminfo.FileSize-iminfo.StripOffsets)/iminfo.StripByteCounts);


% -- Extract header information --
rows = iminfo.Height; % 1
cols = iminfo.Width; % 2

bitdepth = iminfo.BitsPerSample;
if(bitdepth ~= 16)
   error('For movies>4GB only uint16 format is supported at the moment.');
end

sampleFormat = Tiff.SampleFormat.UInt; % Assume Format is 16bit uing ...  
samplesPerPixel = iminfo.SamplesPerPixel; % TODO: For multichannel Tiff. Reading not supported here atm.
if samplesPerPixel > 1
    warning(' Multichannel Tiff not supported! Reading first channel..');
end


% If no interval was given,read whole file
if isempty(minmax_frame)
        minmax_frame = [1, nrFramesOfFile];
% If one number was given, treat as maxframe
elseif numel(minmax_frame)==1 
        minmax_frame = [1,minmax_frame];
% If maxframe is inf, start at minframe and count to end of file
elseif (numel(minmax_frame)==2) 
    if (minmax_frame(1) > nrFramesOfFile) 
       error('Requesting to read from frame %i in a file with %i frames total.',minmax_frame(1), nrFramesOfFile) 
    end
    
    if (minmax_frame(2) == inf)
        minmax_frame(2) = nrFramesOfFile;
    end
    
    if (minmax_frame(2) > nrFramesOfFile) 
       warning('Requesting to read until frame %i in a file with %i frames total. Reading to end of file',minmax_frame(2), nrFramesOfFile) 
    end
end
% nr_frames
nr_frames = minmax_frame(2)-minmax_frame(1)+1;


% -- pre allocate memory according to datatype --
switch sampleFormat
    case Tiff.SampleFormat.IEEEFP
        imgdata = single(zeros(rows, cols, nr_frames));
        sampleFormatName = 'float';
    case Tiff.SampleFormat.UInt
        switch bitdepth
            case 8
                imgdata = uint8(zeros(rows, cols, nr_frames));
            case 16
                imgdata = uint16(zeros(rows, cols, nr_frames));
            case 32
                imgdata = uint32(zeros(rows, cols, nr_frames));
            case 64
                imgdata = uint64(zeros(rows, cols, nr_frames));
            otherwise
                error('Unsupported bitdepth. Supported are 8/16/32/64')
        end
        sampleFormatName = 'uint';
    case Tiff.SampleFormat.Int
        switch bitdepth
            case 8
                imgdata = int8(zeros(rows, cols, nr_frames));
            case 16
                imgdata = int16(zeros(rows, cols, nr_frames));
            case 32
                imgdata = int32(zeros(rows, cols, nr_frames));
            case 64
                imgdata = int64(zeros(rows, cols, nr_frames));
            otherwise
                error('Unsupported bitdepth. Supported are 8/16/32/64')
        end
        sampleFormatName = 'int';
    otherwise
        error('Unsupported datatype ''%s\''. Supported: ''float'', ''uint'', ''int'' ', datatype);
end


% -- File information --
spacing = '    ';
fprintf('%sRows: %i\n',spacing,rows)
fprintf('%sColumns: %i\n',spacing,cols)
fprintf('%sMax Frame: %i\n',spacing,nr_frames)
fprintf('%sBitdepth: %i\n',spacing,bitdepth)
fprintf('%s!!Assumed!! datatype: %s\n',spacing,sampleFormatName)

switch(iminfo.ByteOrder)
    case 'big-endian'
        machineformat = 'ieee-be';
    case 'little-endian'
        machineformat = 'ieee-le';
    otherwise
        error('Unknown ByteOrder:%s',iminfo.ByteOrder)
end


%  -- Rewind to the first frame that should be read --
bytesPerFrame = iminfo.StripByteCounts;

fID = fopen(filename,'r',machineformat);
fseek(fID, iminfo.StripOffsets, 'bof'); % Skip offset
fseek(fID, (minmax_frame(1)-1)*bytesPerFrame, 'cof'); % Skip to the first frame to read


msgAccumulator = ''; % needed for one-line printing

% -- read rest of the stack --
startTime = tic;
lastElapsedTime = 0;

for iFrame = 1:nr_frames
    elapsedTime = toc(startTime);
    if( (elapsedTime-lastElapsedTime)>0.25) % output every 0.25 seconds
        rewindMessages()
        rewPrintf('\nReading Frame %i/%i ..',iFrame, nr_frames);
        lastElapsedTime = elapsedTime;
    end
    imgdata(:,:,iFrame) = fread(fID, [rows,cols], '*uint16', 0, machineformat)';
end
fclose(fID);
rewindMessages()
rewPrintf('\nReading Frame %i/%i .. done\n',nr_frames, nr_frames);

if convert2double
    imgdata = double(imgdata);
end

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
warning('on','all'); % Restore warnings
status = fclose(fID);
%     fprintf('\nUSER EXIT');
end