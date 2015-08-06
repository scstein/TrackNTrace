classdef DetectionProcess < ImageAnalysisProcess
    % An abstract class for all detection processes with an output
    % structure compatible with the tracker
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
    
    % Chuangang Ren, 11/2010
    % Sebastien Besson (last modified May 2012)
    
    methods(Access = public)
        
        function obj = DetectionProcess(owner, name, funName, funParams )
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;
            end
            if nargin > 3
                obj.funParams_ = funParams;
            end
        end
        
        function status = checkChannelOutput(obj,iChan)
            
            %Checks if the selected channels have valid output files
            nChan = numel(obj.owner_.channels_);
            if nargin < 2 || isempty(iChan), iChan = 1:nChan; end
            
            status=  ismember(iChan,1:nChan) & ....
                arrayfun(@(x) exist(obj.outFilePaths_{1,x},'file'),iChan);
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'movieInfo'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            s = load(obj.outFilePaths_{1,iChan},output{:});
           
            if numel(ip.Results.iFrame)>1,
                varargout{1}=s.(output{1});
            else
                varargout{1}=s.(output{1})(iFrame);
            end
        end
        function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Objects';
            output(1).var='movieInfo';
            output(1).formatData=@DetectionProcess.formatOutput;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none','Color',colors(x,:));
        end  
        
    end
    methods(Static)
        function name = getName()
            name = 'Detection';
        end
        function h = GUI()
            h = @abstractProcessGUI;
        end
        function procClasses = getConcreteClasses()
            procClasses = ...
                {@SubResolutionProcess;
                @CometDetectionProcess;
                @AnisoGaussianDetectionProcess;
                @NucleiDetectionProcess;
                @PointSourceDetectionProcess;
                };
            procClasses = cellfun(@func2str, procClasses, 'Unif', 0);
        end
        
        function y =formatOutput(x)
            % Format output in xy coordinate system
            if isempty(x.xCoord)
                y = NaN(1,2);
            else
                y = horzcat(x.xCoord(:,1),x.yCoord(:,1));
            end
        end
    end
end