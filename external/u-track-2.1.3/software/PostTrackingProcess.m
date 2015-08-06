classdef PostTrackingProcess < DataProcessingProcess
    % An abstract class for all post-tracking processes 
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
    
    % Sebastien Besson Jul 2012
    
    methods(Access = public)
        
        function obj = PostTrackingProcess(owner, name, funName, funParams )
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;
            end
            if nargin > 3
                obj.funParams_ = funParams;
            end
        end
        
    end
    methods(Static)
        function name = getName()
            name = 'Track analysis';
        end
        function h = GUI()
            h = @abstractProcessGUI;
        end
        function procClasses = getConcreteClasses()
            procClasses = ...
                {'MotionAnalysisProcess';
                'CometPostTrackingProcess';                
                'TrackGroupingProcess'};
        end
    end
end