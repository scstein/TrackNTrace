%setErrorbarStyle(he, pos, de) modifies the width of the error bars
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

% Francois Aguet, Feb 22 2011 (last modif. 01/21/2012)

function setErrorbarStyle(he, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('he');
ip.addOptional('de', 0.2, @isscalar);
ip.addParamValue('Position', 'both', @(x) any(strcmpi(x, {'both', 'top', 'bottom'})));
ip.parse(he, varargin{:});
de = ip.Results.de;

he = get(he, 'Children');
xd = get(he(2), 'XData');
if strcmpi(ip.Results.Position, 'bottom')
    xd(4:9:end) = xd(1:9:end);
    xd(5:9:end) = xd(1:9:end);
else
    xd(4:9:end) = xd(1:9:end) - de;
    xd(5:9:end) = xd(1:9:end) + de;
end
if strcmpi(ip.Results.Position, 'top')
    xd(7:9:end) = xd(1:9:end);
    xd(8:9:end) = xd(1:9:end);
else
    xd(7:9:end) = xd(1:9:end) - de;
    xd(8:9:end) = xd(1:9:end) + de;
end
set(he(2), 'XData', xd);