classdef  Reader < handle
    % Concrete implementation of MovieObject for a single movie
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
    
    properties
        sizeX
        sizeY
        sizeZ
        sizeC
        sizeT
        bitDepth
    end
    
    methods(Abstract)
        getSizeX(obj)
        getSizeY(obj)
        getSizeZ(obj)
        getSizeC(obj)
        getSizeT(obj)
        getBitDepth(obj)
        getImageFileNames(obj)
        getChannelNames(obj)
        loadImage(obj)
    end
    
    methods ( Access = public )
        
        function I = loadStack(obj, c, t, varargin)
        % loadStack reads a Z-stack and returns a YXZ Matrix
        %
        %   loadStack(c,t,Z) returns a Z-stack for channel c and
        %   time frame t  with Z planes indicated by the vector Z
        %
        %   loadStack(c,t) returns loadStack(c,t, 1: getSizeZ())
        %
        % output: a 3D matrix with dimensions YXZ
        %
        % Note: Override only if need you need to overload or change how
        % the input is checked. Otherwise override loadStack_.
        %
        % Example:
        %   reader = movieData.getReader();
        %   % The following are all equivalent
        %   zStackMatrix = reader.loadStack(1,1,1:reader.getSizeZ());
        %   zStackMatrix = reader.loadStack(1,1);
        %   zStackMatrix = reader.loadStack(1);
        %
       
            % Input check
            ip = inputParser;
            ip.addRequired('c', ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC() ) );
            % Parse c first for validation check
            %  since getSizeT and getSizeZ may depend on a valid c.
            ip.parse(c);
            ip.addRequired('t', ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT() ) );
            ip.addOptional('z', 1 : obj.getSizeZ(), ...
                @(x) all(ismember(x, 1 : obj.getSizeZ() ) ) );
            ip.parse(c, t, varargin{:});
            t = ip.Results.t;
            z = ip.Results.z;
            
            I = obj.loadStack_(c,t,z);
        end
    end
    methods( Access = protected )
        function I = loadStack_(obj, c, t, z)
        % loadStack_ is a validation-free version called by loadStack
        %
        % Provides a generic implementation of loadStack.
        % Guarantees all readers will have a loadStack method if they have
        % a loadImage method.
        %
        % Override if there are backend optimizations.
        % This is likely slower because of repeated input validation.
        %
        % parameters are c, t, z with only c required
        % t defaults to 1
        % z defaults to 1:sizeZ
        %
        
            % Load first image to get class and dimensions.
            first = obj.loadImage_( c , t , z(1) );
            % "I" will be a YXZ matrix.
            I = zeros( [ size(first) length(z) ] , class(first));
            I(:,:,1) = first;
            
            for zi = 2:length(z)
                I(:,:,zi) = obj.loadImage_( c , t , z(zi) );
            end
                
        end
    %end
    %methods( Abstract, Access = protected )
        % I = loadImage_(obj, c, t, z )
        function I = loadImage_(obj, c, t, z )
            % loadImage_ is a validation-free implementation called by
            %   loadImage (or will be in the near future)
            %
            % Forward to validation version for now.
            % This WILL be an abstract method in the future.

            I = obj.loadImage(c,t,z);
        end
    end
    
end
