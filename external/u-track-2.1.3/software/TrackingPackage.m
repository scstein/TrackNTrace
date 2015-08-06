classdef TrackingPackage < Package
    % An abstract class for a geeneric Tracking Package
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
    
    methods
        function obj = TrackingPackage(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'TrackingPackage'];
            end
                 
            % Call the superclass constructor
            obj = obj@Package(super_args{:});        
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'U-Track';
        end 
        function m = getDependencyMatrix(i,j)   
            m = [0 0 0;  %1 DetectionProcess
                 1 0 0;  %2 TrackingProcess
                 1 1 0;];%3 PostTrackingProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = trackingPackageGUI(varargin{:});
        end
        function classes = getProcessClassNames(index)
            classes = {
                'DetectionProcess',...
                'TrackingProcess',...
                'PostTrackingProcess'};
            if nargin==0, index=1:numel(classes); end
            classes=classes(index);
        end
        
        function objects = getConcretePackages()
            objects(1).name = 'Single particles';
            objects(1).packageConstr = @UTrackPackage;
            objects(2).name = 'Microtubules plus-ends';
            objects(2).packageConstr = @PlusTipTrackerPackage;                        
            objects(3).name = 'Nuclei';
            objects(3).packageConstr = @NucleiTrackingPackage;
        end
        
    end 
end