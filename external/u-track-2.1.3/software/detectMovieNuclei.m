function detectMovieNuclei(movieData,varargin)
% detectMovieNuclei detects nuclei objects in a movie
%
% detectMovieNuclei 
%
% SYNOPSIS detectMovieNuclei(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   
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

% Sebastien Besson, Nov 2012

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('NucleiDetectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(NucleiDetectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
nucDetProc = movieData.getProcess(iProc);
%Parse input, store in parameter structure
p = parseProcessParams(nucDetProc,paramsIn);

%% --------------- Initialization ---------------%%

if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',nucDetProc.getName());
else
    wtBar=-1;
end

% Set up the input directories
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    if isempty(p.ProcessIndex)
        inFilePaths{1,i} = movieData.getChannelPaths{i};
    else
       inFilePaths{1,i} = movieData.processes_{p.ProcessIndex}.outFilePaths_{1,i}; 
    end
end
nucDetProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'nuclei_channel' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory);
nucDetProc.setOutFilePaths(outFilePaths);

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting nuclei...')

logMsg = @(chan) ['Please wait, detecting nuclei for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;

nFrames=movieData.nFrames_;
nChan = length(p.ChannelIndex);
nDetFrames=p.lastFrame-p.firstFrame+1;
nTot = nChan*nDetFrames;

movieInfo(nFrames,1)=deal(struct('xCoord',[],'yCoord',[],'amp',[]));

for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
    disp(logMsg(iChan))
    disp(inFilePaths{1,iChan});
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    if ishandle(wtBar), waitbar((i-1)/numel(p.ChannelIndex),wtBar,logMsg(iChan)); end
    
    for j=p.firstFrame:p.lastFrame
        if isempty(p.ProcessIndex)
            I = double(movieData.channels_(iChan).loadImage(j));
        else
            I = double(movieData.getProcess(p.ProcessIndex).loadChanelOutput(iChan, j));
        end
        movieInfo(j,1) = detectNuclei(I, p.radius, p.confluent,...
            'edgeFilter', p.edgeFilter, 'sigma', p.sigma, 'p', p.p,...
            'useDblLog', p.useDblLog);
        
        % Update the waitbar
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = (i-1)*nDetFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
    
    % Save results
    save(outFilePaths{1,iChan}, 'movieInfo');
    
end

if ishandle(wtBar), close(wtBar); end
disp('Finished detecting nuclei...')
