function imgdata = read_tiff(filename, convert2double, allwarn)
% Usage: imgdata = read_tiff(filename, convert2double, allwarn)
%
% Loads an image or image stack from a Tiff file.
%
% filename: Tiff file to read
% convert2double: Datatype conversion to double after reading | default: true
% allwarn : Show warnings for all slices, not only first one | default: false
%
% NOTE: Deactivate double conversion if not neccessary, reading the same datatype used in the Tiff file.
%       This potentially saves a lot of memory (factor 4 for 16-bit Tiff movies) and thus increased performance.
%       You can still perform a cast to double (or other datatypes) later if desired.
%         Example: imgdata = double( read_tiff(filename) );

% By
% Simon Christoph Stein - August 2014
% E-Mail: scstein@phys.uni-goettingen.de
%

fprintf('Reading file ''%s''.. \n', filename)

if nargin < 2 || isempty(convert2double)
    convert2double = true;
end

if nargin < 3 || isempty(allwarn)
    allwarn = false;
end

% -- Create Tiff object --
% Check filename, append '.tif' if neccessary
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext, '.tif') == 0
    filename = [filename '.tif'];
end
t = Tiff(filename,'r');
cleanupTrigger = onCleanup(@() cleanupFunc(t));

% -- Extract header information --
rows = t.getTag('ImageLength'); % 1
cols = t.getTag('ImageWidth'); % 2
% t.getTag('ImageDepth'); % what is this for?

bitdepth = t.getTag('BitsPerSample');
sampleFormat = t.getTag('SampleFormat');
samplesPerPixel = t.getTag('SamplesPerPixel'); % TODO: For multichannel Tiff. Reading not supported here atm.
if samplesPerPixel > 1
    warning(' Multichannel Tiff not supported! Reading first channel..');
end



% -- Count number of frames for pre-allocation --
not(t.lastDirectory); % this forces the Tiff object to read the first frames tag data and check for errors.
if not(allwarn), warning('off','all');  end; % disable warnings for counting
while not(t.lastDirectory)
    t.nextDirectory;
end
nr_frames = t.currentDirectory;
warning('on','all');


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
fprintf('%sFrames: %i\n',spacing,nr_frames)
fprintf('%sBitdepth: %i\n',spacing,bitdepth)
fprintf('%sDatatype: %s\n',spacing,sampleFormatName)

%  -- Rewind to the first frame and read it --
t.setDirectory(1);
msgAccumulator = ''; % needed for one-line printing
rewPrintf('Reading Frame %i/%i ..',1, nr_frames);
imgdata(:,:,1) = t.read();

% -- read rest of the stack --
startTime = tic;
lastElapsedTime = 0;

if not(allwarn), warning('off','all');  end; % Disable warnings after first slice (we could miss something important here, but usually they just repeat)
for IFD = 2:nr_frames
    elapsedTime = toc(startTime);    
    if( (elapsedTime-lastElapsedTime)>0.25) % output every 0.25 seconds
        rewindMessages()
        rewPrintf('\nReading Frame %i/%i ..',IFD, nr_frames);
        lastElapsedTime = elapsedTime;
    end
    
    
    t.nextDirectory;
    imgdata(:,:,IFD) = t.read();
end
rewindMessages()
rewPrintf('\nReading Frame %i/%i .. done\n',nr_frames, nr_frames);
warning('on','all'); % Restore warnings

t.close;


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
%     fprintf('\nUSER EXIT');
end