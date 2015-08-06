function trackMovie(movieData,varargin)
% Track features in a movie which has been processed by a detection method
%
% Sebastien Besson, 5/2011
%
% Copyright (C) 2014 LCCB 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;



%Get the indices of any previous tracking processes from this function                                                                              
iProc = movieData.getProcessIndex('TrackingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(TrackingProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

trackProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(trackProc,paramsIn);

%% --------------- Initialization ---------------%%

% Check detection process first
if isempty(p.DetProcessIndex)
    p.DetProcessIndex =movieData.getProcessIndex('DetectionProcess',1,1);

    if isempty(p.DetProcessIndex)
        error(['Detection has not been run! '...
            'Please run detection prior to tracking!'])
    end
end
detProc=movieData.processes_{p.DetProcessIndex};

if ~detProc.checkChannelOutput(p.ChannelIndex)
    error(['Missing detection output ! Please apply detection before ' ...
        'running tracking!'])
end

% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = detProc.outFilePaths_{1,i};
end
trackProc.setInFilePaths(inFilePaths);
    
% Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilename= ['Channel_' num2str(i) '_tracking_result'];
    outFilePaths{1,i} = [p.OutputDirectory filesep outFilename '.mat'];
    if p.saveResults.export
        outFilePaths{2,i} = [p.OutputDirectory filesep outFilename '_mat.mat'];
    end
end
mkClrDir(p.OutputDirectory);
trackProc.setOutFilePaths(outFilePaths);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting tracking...')

for i = p.ChannelIndex
    movieInfo = detProc.loadChannelOutput(i);
    
    % Call function - return tracksFinal for reuse in the export
    % feature
    tracksFinal = trackCloseGapsKalmanSparse(movieInfo, p.costMatrices, p.gapCloseParam,...
        p.kalmanFunctions, p.probDim, 0, p.verbose);
    save(outFilePaths{1,i},'tracksFinal');
    
    % Optional export
    if p.saveResults.export
        if ~p.gapCloseParam.mergeSplit
            [M.trackedFeatureInfo M.trackedFeatureIndx]=...
                convStruct2MatNoMS(tracksFinal);
        else
            [M.trackedFeatureInfo M.trackedFeatureIndx,M.trackStartRow,M.numSegments]=...
                convStruct2MatIgnoreMS(tracksFinal);
        end
        save(outFilePaths{2,i},'-struct','M');
        clear M;
    end
end

disp('Finished tracking!')
