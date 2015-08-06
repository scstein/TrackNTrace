classdef TrackingProcess < DataProcessingProcess
    % A class definition for a generic tracking process.
    %
    % Chuangang Ren, 11/2010
    % Sebastien Besson (last modified Dec 2011)
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
    
    
    methods(Access = public)
        
        function obj = TrackingProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = TrackingProcess.getName;
                super_args{3} = @trackMovie;
                if isempty(funParams)
                    funParams = TrackingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'tracksFinal', 'gapInfo', 'staticTracks'};
            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan', @(x) obj.checkChanNum(x));
            ip.addOptional('iFrame', [] ,@(x) obj.checkFrameNum(x));
            ip.addParamValue('output', outputList{1}, @(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
            output = ip.Results.output;
            iFrame = ip.Results.iFrame;
            if ischar(output),output={output}; end
            
            % Data loading
            s = load(obj.outFilePaths_{iChan}, 'tracksFinal');
            
            varargout = cell(numel(output), 1);
            for i = 1:numel(output)
                switch output{i}
                    case {'tracksFinal', 'staticTracks'}
                        varargout{i} = s.tracksFinal;
                    case 'gapInfo'
                        varargout{i} = findTrackGaps(s.tracksFinal);
                end
                if strcmp(output{i}, 'tracksFinal') && ~isempty(iFrame),
                    % Filter tracks existing in input frame
                    trackSEL=getTrackSEL(s.tracksFinal);
                    validTracks = (iFrame>=trackSEL(:,1) &iFrame<=trackSEL(:,2));
                    [varargout{i}(~validTracks).tracksCoordAmpCG]=deal([]);
                    
                    for j=find(validTracks)'
                        varargout{i}(j).tracksCoordAmpCG = varargout{i}(j).tracksCoordAmpCG(:,1:8*(iFrame-trackSEL(j,1)+1));
                    end
                end
            end
        end
        
        function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Tracks';
            output(1).var='tracksFinal';
            output(1).formatData=@TrackingProcess.formatTracks;
            output(1).type='overlay';
            output(1).defaultDisplayMethod = @(x) TracksDisplay(...
                'Color',colors(x,:));
            output(2).name='Gap length histogram';
            output(2).var='gapInfo';
            output(2).formatData=@(x) x(:,4);
            output(2).type='graph';
            output(2).defaultDisplayMethod=@(x)HistogramDisplay('XLabel','Gap length',...
                'YLabel','Counts');
            output(3).name='Static tracks';
            output(3).var='staticTracks';
            output(3).formatData=@TrackingProcess.formatTracks;
            output(3).type='overlay';
            output(3).defaultDisplayMethod = @(x) TracksDisplay(...
                'Color',colors(x,:), 'useDragtail', false);
        end
        
        
    end
    methods(Static)
        function name = getName()
            name = 'Tracking';
        end
        function h = GUI()
            h= @trackingProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            
            funParams.ChannelIndex =1:numel(owner.channels_);
            funParams.DetProcessIndex = [];
            funParams.OutputDirectory = [outputDir  filesep 'tracks'];
            
            % --------------- gapCloseParam ----------------
            funParams.gapCloseParam.timeWindow = 5; %IMPORTANT maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
            funParams.gapCloseParam.mergeSplit = 0; % (SORT OF FLAG: 4 options for user) 1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            funParams.gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.
            funParams.gapCloseParam.diagnostics = 1; %FLAG 1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
            
            
            % --------------- kalmanFunctions ----------------
            
            kalmanFunctions = TrackingProcess.getKalmanFunctions(1);
            fields = fieldnames(kalmanFunctions);
            validFields = {'reserveMem','initialize','calcGain','timeReverse'};
            kalmanFunctions = rmfield(kalmanFunctions,fields(~ismember(fields,validFields)));
            funParams.kalmanFunctions = kalmanFunctions;
            
            % --------------- saveResults ----------------
            funParams.saveResults.export = 0; %FLAG allow additional export of the tracking results into matrix
            
            % --------------- Others ----------------
            
            funParams.verbose = 1;
            funParams.probDim = 2;
            
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            
            
        end
        
        function kalmanFunctions = getKalmanFunctions(index)
            % Brownian + Directed motion models
            kalmanFunctions(1).name = 'Brownian + Directed motion models';
            kalmanFunctions(1).reserveMem = func2str(@kalmanResMemLM);
            kalmanFunctions(1).initialize  = func2str(@kalmanInitLinearMotion);
            kalmanFunctions(1).initializeGUI  = @kalmanInitializationGUI;
            kalmanFunctions(1).calcGain    = func2str(@kalmanGainLinearMotion);
            kalmanFunctions(1).timeReverse = func2str(@kalmanReverseLinearMotion);
            
            % Microtubule plus-end dynamics
            kalmanFunctions(2).name = 'Microtubule plus-end dynamics';
            kalmanFunctions(2).reserveMem = func2str(@kalmanResMemLM);
            kalmanFunctions(2).initialize  = func2str(@plusTipKalmanInitLinearMotion);
            kalmanFunctions(2).initializeGUI  = @kalmanInitializationGUI;
            kalmanFunctions(2).calcGain    = func2str(@plusTipKalmanGainLinearMotion);
            kalmanFunctions(2).timeReverse = func2str(@kalmanReverseLinearMotion);
            
            if nargin>0
                assert(all(ismember(index, 1:numel(kalmanFunctions))));
                kalmanFunctions=kalmanFunctions(index);
            end
        end
        
        
        function costMatrix = getDefaultLinkingCostMatrices(owner,timeWindow,varargin)
            
            % Brownian + Directed motion models
            costMatrices(1).name = 'Brownian + Directed motion models';
            costMatrices(1).funcName = func2str(@costMatRandomDirectedSwitchingMotionLink);
            costMatrices(1).GUI = @costMatRandomDirectedSwitchingMotionLinkGUI;
            costMatrices(1).parameters.linearMotion = 0; % use linear motion Kalman filter.
            costMatrices(1).parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            costMatrices(1).parameters.maxSearchRadius = 5; %IMPORTANT maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            costMatrices(1).parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
            costMatrices(1).parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            costMatrices(1).parameters.nnWindow = timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
            costMatrices(1).parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
            costMatrices(1).parameters.diagnostics = owner.nFrames_-1;
            
            % plusTip markers
            plusTipCostMatrix.name = 'Microtubule plus-end dynamics';
            plusTipCostMatrix.funcName = func2str(@plusTipCostMatLinearMotionLink);
            plusTipCostMatrix.GUI = @plusTipCostMatLinearMotionLinkGUI;
            plusTipCostMatrix.parameters.linearMotion = 1; % use linear motion Kalman filter.
            plusTipCostMatrix.parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            plusTipCostMatrix.parameters.maxSearchRadius = 10; %IMPORTANT maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            plusTipCostMatrix.parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
            plusTipCostMatrix.parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            plusTipCostMatrix.parameters.nnWindow = timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
            plusTipCostMatrix.parameters.kalmanInitParam.initVelocity = []; %Kalman filter initialization parameters.
            plusTipCostMatrix.parameters.kalmanInitParam.convergePoint = []; %Kalman filter initialization parameters.
            plusTipCostMatrix.parameters.kalmanInitParam.searchRadiusFirstIteration = 20; %Kalman filter initialization parameters.
            plusTipCostMatrix.parameters.diagnostics = [];
            costMatrices(2)=plusTipCostMatrix;
            
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addRequired('timeWindow',@isscalar);
            ip.addOptional('index',1:length(costMatrices),@isvector);
            ip.parse(owner,timeWindow,varargin{:});
            index = ip.Results.index;
            costMatrix=costMatrices(index);
        end
        
        function costMatrix = getDefaultGapClosingCostMatrices(owner,timeWindow,varargin)
            
            % Linear motion
            costMatrices(1).name = 'Brownian + Directed motion models';
            costMatrices(1).funcName = func2str(@costMatRandomDirectedSwitchingMotionCloseGaps);
            costMatrices(1).GUI = @costMatRandomDirectedSwitchingMotionCloseGapsGUI;
            costMatrices(1).parameters.linearMotion = 0; %use linear motion Kalman filter.
            
            costMatrices(1).parameters.minSearchRadius = 2; %minimum allowed search radius.
            costMatrices(1).parameters.maxSearchRadius = 5; %maximum allowed search radius.
            costMatrices(1).parameters.brownStdMult = 3*ones(timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
            
            costMatrices(1).parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
            costMatrices(1).parameters.nnWindow = timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
            costMatrices(1).parameters.brownScaling = [0.5 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
            costMatrices(1).parameters.timeReachConfB = timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
            costMatrices(1).parameters.ampRatioLimit = [0.5 2]; % (FLAG + VALUES small-big value) for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
            
            % If parameters.linearMotion = 1
            costMatrices(1).parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
            costMatrices(1).parameters.linStdMult = 3*ones(timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
            costMatrices(1).parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
            costMatrices(1).parameters.timeReachConfL = timeWindow; %similar to timeReachConfB, but for the linear part of the motion.
            costMatrices(1).parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
            % ---------------------------------
            
            costMatrices(1).parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
            costMatrices(1).parameters.resLimit = []; % text field resolution limit, which is generally equal to 3 * point spread function sigma.
            
            % Linear motion
            plusTipCostMatrix.name = 'Microtubule plus-end dynamics';
            plusTipCostMatrix.funcName = func2str(@plusTipCostMatCloseGaps);
            plusTipCostMatrix.GUI = @plusTipCostMatCloseGapsGUI;
            plusTipCostMatrix.parameters.maxFAngle = 30; %use linear motion Kalman filter.
            plusTipCostMatrix.parameters.maxBAngle = 10; %use linear motion Kalman filter.
            plusTipCostMatrix.parameters.backVelMultFactor = 1.5;
            plusTipCostMatrix.parameters.fluctRad = 1.0;
            plusTipCostMatrix.parameters.breakNonLinearTracks = false;
            costMatrices(2)=plusTipCostMatrix;
            
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addRequired('timeWindow',@isscalar);
            ip.addOptional('index',1:length(costMatrices),@isvector);
            ip.parse(owner,timeWindow,varargin{:});
            index = ip.Results.index;
            costMatrix=costMatrices(index);
        end
        
        function displayTracks = formatTracks(tracks)
            % Format tracks structure into compound tracks for display
            % purposes
            
            % Determine the number of compound tracks
            nCompoundTracks = cellfun('size',{tracks.tracksCoordAmpCG},1)';
            nTracksTot = [0 cumsum(nCompoundTracks(:))'];
            nTracks = nTracksTot(end);
            
            % Fail fast if no track
            if nTracks == 0
                displayTracks = struct.empty(1,0);
                return
            end
            
            displayTracks(sum(nCompoundTracks),1) = struct('xCoord', [], 'yCoord', [], 'events' ,[], 'number', []);
            for i = find(nCompoundTracks)'
                % Get the x and y coordinate of all compound tracks
                for  j = 1 : nCompoundTracks(i)
                    iTrack = nTracksTot(i) + j ;
                    displayTracks(iTrack).xCoord = tracks(i).tracksCoordAmpCG(j, 1:8:end);
                    displayTracks(iTrack).yCoord = tracks(i).tracksCoordAmpCG(j, 2:8:end);
                    displayTracks(iTrack).number = i;
                    if isfield(tracks, 'label'),
                        displayTracks(iTrack).label = tracks(i).label;
                    end
                end
                
                % Fill split events NaNs in compound tracks
                nTimes = numel(displayTracks(iTrack).xCoord);
                splitEvents = find(tracks(i).seqOfEvents(:,2)==1 & ~isnan(tracks(i).seqOfEvents(:,4)))';
                eventTimes = tracks(i).seqOfEvents(splitEvents,1)-1;
                for iEvent = splitEvents(eventTimes < nTimes)
                    iTrack1 = nTracksTot(i)+ tracks(i).seqOfEvents(iEvent,3);
                    iTrack2 = nTracksTot(i)+ tracks(i).seqOfEvents(iEvent,4);
                    t = tracks(i).seqOfEvents(iEvent,1)-1;
                    displayTracks(iTrack1).xCoord(t) = displayTracks(iTrack2).xCoord(t);
                    displayTracks(iTrack1).yCoord(t) = displayTracks(iTrack2).yCoord(t);
                    displayTracks(iTrack1).events(t) = 's';
                    displayTracks(iTrack2).events(t) = 's';
                end
                
                % Fill merge events NaNs in compound tracks
                mergeEvents = find(tracks(i).seqOfEvents(:,2)==2 & ~isnan(tracks(i).seqOfEvents(:,4)))';
                eventTimes = tracks(i).seqOfEvents(mergeEvents,1)-1;
                for iEvent = mergeEvents(eventTimes < nTimes)
                    iTrack1 = nTracksTot(i)+ tracks(i).seqOfEvents(iEvent,3);
                    iTrack2 = nTracksTot(i)+ tracks(i).seqOfEvents(iEvent,4);
                    t = tracks(i).seqOfEvents(iEvent,1);
                    displayTracks(iTrack1).xCoord(t) = displayTracks(iTrack2).xCoord(t);
                    displayTracks(iTrack1).yCoord(t) = displayTracks(iTrack2).yCoord(t);
                    displayTracks(iTrack1).events(t) = 'm';
                    displayTracks(iTrack2).events(t) = 'm';

                end
            end
        end
    end
end

