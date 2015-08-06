classdef MotionAnalysisProcess < PostTrackingProcess
    % A concrete class for analyzing tracks diffusion
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
    
    % Sebastien Besson, March 2012
    
    methods (Access = public)
        function obj = MotionAnalysisProcess(owner, varargin)
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
                
                super_args{1} = owner;
                super_args{2} = MotionAnalysisProcess.getName;
                super_args{3} = @analyzeMovieMotion;
                if isempty(funParams)  % Default funParams
                    funParams = MotionAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@PostTrackingProcess(super_args{:});
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'diffAnalysisRes', 'tracks'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            nOutput = numel(output);
            
            % Data loading
            s = load(obj.outFilePaths_{1,iChan},output{:});
            
            varargout = cell(nOutput);
            for i = 1:nOutput
                switch output{i}
                    case 'tracks'
                        tracksFinal = s.(output{i});
                        if ~isempty(iFrame),
                            % Filter tracks existing in input frame
                            trackSEL=getTrackSEL(tracksFinal);
                            validTracks = (iFrame>=trackSEL(:,1) &iFrame<=trackSEL(:,2));
                            [tracksFinal(~validTracks).tracksCoordAmpCG]=deal([]);
                            
                            for j=find(validTracks)'
                                tracksFinal(j).tracksCoordAmpCG = tracksFinal(j).tracksCoordAmpCG(:,1:8*(iFrame-trackSEL(j,1)+1));
                            end
                            varargout{i} = tracksFinal;
                        else
                            varargout{i} = tracksFinal;
                        end
                    case 'diffAnalysisRes'
                        varargout{i} = s.(output{i});
                end
            end
        end
        
        function output = getDrawableOutput(obj)
            types = MotionAnalysisProcess.getTrackTypes();
            colors = vertcat(types.color);
            output(1).name='Classified tracks';
            output(1).var='tracks';
            output(1).formatData=@MotionAnalysisProcess.formatTracks;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)TracksDisplay('Color', colors);
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'Motion Analysis';
        end
        function h = GUI()
            h = @motionAnalysisProcessGUI;
        end
        
        function alpha = getAlphaValues()
            alpha=[0.01 0.05 0.1 0.2];
        end
        
        function methods = getConfinementRadiusMethods()
            methods(1).type = 0;
            methods(1).name = 'Mean positional standard deviation';
            methods(2).type = 1;
            methods(2).name = 'Minimum positional standard deviation';
            methods(3).type = 2;
            methods(3).name = 'Rectangle approximation';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'MotionAnalysis'];
            funParams.probDim = 2;
            funParams.checkAsym = 0;
            funParams.alphaValues = [0.05 0.1];
            funParams.confRadMin=0;
        end
        
        function displayTracks = formatTracks(tracks)
            % Format classified tracks into structure for display
            
            % Read track types and classification matrix
            types = MotionAnalysisProcess.getTrackTypes();
            track_class = vertcat(tracks.classification);
            
            % Assign label to each track of known type
            for i = 1 : numel(types) - 1
                idx = types(i).f(track_class);
                if any(idx), [tracks(idx).label] = deal(i); end
            end
            
            % Assign last label to unlabelled tracks
            if ~isfield(tracks, 'label')
                 [tracks.label] = deal(numel(types));
            else
                idx = cellfun(@isempty, {tracks.label});
                [tracks(idx).label] = deal(numel(types));
            end
            
            % Format tracks using TrackingProcess utility function
            displayTracks = TrackingProcess.formatTracks(tracks);
        end
        
        function types = getTrackTypes()
            % Get the color map for classified tracks
            %
            % see also: plotTracksDiffAnalysis2D
            
            types(1).name = 'linear & 1D confined diffusion';
            types(1).f = @(x) x(:, 1) == 1 & x(:, 3) == 1;
            types(1).color = [1 0.7 0];
            types(2).name = 'linear & 1D normal diffusion';
            types(2).f = @(x) x(:, 1) == 1 & x(:, 3) == 2;
            types(2).color = [1 0 0];
            types(3).name = 'linear & 1D super diffusion';
            types(3).f = @(x) x(:, 1) == 1 & x(:, 3) == 3;
            types(3).color = [0 1 0];
            types(4).name = 'linear & too short to analyze 1D diffusion';
            types(4).f = @(x) x(:, 1) == 1 & isnan(x(:, 3));
            types(4).color = [1 1 0];
            types(5).name = 'random/unclassified & 2D confined diffusion';
            types(5).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 1;
            types(5).color = [0 0 1];
            types(6).name = 'random/unclassified & 2D normal diffusion';
            types(6).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 2;
            types(6).color = [0 1 1];
            types(7).name = 'random/unclassified & 2D super diffusion';
            types(7).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 3;
            types(7).color = [1 0 1];
            types(8).name = 'random & too short to analyze 2D diffusion';
            types(8).f = @(x) x(:, 1) == 0 & isnan(x(:, 2));
            types(8).color = [.6 0 1];
            types(9).name = 'too short for any analysis';
            types(9).f = @(x) 1;
            types(9).color = [.7 .7 .7];
        end
        
    end
end