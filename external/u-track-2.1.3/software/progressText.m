function progressText(fractionDone,text)
%PROGRESSTEXT shows progress of a loop as text on the screen
%
% SYNOPSIS: progressText(fractionDone,text)
%
% INPUT fractionDone: fraction of loop done (0-1)
%		text (opt): {yourTexthere} : XX% done xx:xx:xx remaining
%                   Note: text can be changed for every call
% OUTPUT none
%
% EXAMPLE
%   n = 1000;
%   progressText(0,'Test run') % Create text
%   for i = 1:n
%       pause(0.01) % Do something important
%       progressText(i/n) % Update text
%   end
%
% REMARKS progressText will set lastwarn to ''
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn based on progressbar.m
% DATE: 29-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

persistent starttime lastupdate clearText printText finalText warnText clearWarn

% constants
nCharsBase = 27; % change this if changing output format

% Test input
if nargin < 1 || isempty(fractionDone) 
    fractionDone = 0;
end
if nargin < 2 || isempty(text)
    text = '';
else
    text = [text,' : '];
end

if fractionDone == 0 || isempty(starttime)
    % set up everything
    
    % Set time of last update to ensure calculation
    lastupdate = clock - 1;
    
    % Task starting time reference
    starttime = clock;
    
    % create fprintf-expression
    
    printText = sprintf('%s%%2d%%%% done %%s remaining',text);
    initialText = sprintf('%s 0%%%% done xx:xx:xx remaining',text);
    finalText = sprintf('%s100%%%% done %%s elapsed',text);
    % get length of fprintf expression
    nChars = nCharsBase + length(text);
    clearText = repmat('\b',1,nChars);
    
    % print initialText and return
    fprintf(1,initialText);
    
%     %fprintfExpression removes old expression before overwriting
%     fprintfExpression = [clearText printText];
%     fprintfExpressionFinal = [clearText, finalText];

% empty warning
lastwarn('');
warnText = '';
clearWarn = '';
    
    return
elseif ~isempty(text)
    % text has been changed. Create fprintfExpressions first, then update
    % clearText
    printText = sprintf('%s%%2d%%%% done %%s remaining',text);
    finalText = sprintf('%s100%%%% done %%s elapsed',text);
    fprintfExpression = [clearText clearWarn printText, warnText];
    fprintfExpressionFinal = [clearText, clearWarn, finalText, warnText, '\n'];
    
    nChars = nCharsBase + length(text);
    clearText = repmat('\b',1,nChars);
elseif ~isempty(lastwarn)
    % add warnings to the end of the progressText
    % find warning
    w = lastwarn;
    lastwarn('')
    nw = length(w);
    % erase warning
    fprintf(1,repmat('\b',1,11+nw)); %--- is this correct???
    % create new warnText
    w = regexprep(w,['(',char(10),'\s+)'],' - ');
    warnTextNew = sprintf('\n  Warning @%7.3f%%%% done : %s ',100*fractionDone,w);
    warnLength = length(warnTextNew) - 1; % subtract 2 for %%-sign
    warnText = [warnText,warnTextNew];
    
    fprintfExpression = [clearText clearWarn printText, warnText];
    fprintfExpressionFinal = [clearText, clearWarn, finalText, warnText, '\n'];
    
    % prepare cleanWarn for next time
    clearWarn = [clearWarn,repmat('\b',1,warnLength)];
else
    % all is normal. Just generate output
    fprintfExpression = [clearText clearWarn printText, warnText];
    fprintfExpressionFinal = [clearText, clearWarn, finalText, warnText, '\n'];
end

% write progress
percentDone = floor(100*fractionDone);

% get elapsed time
runTime = etime(clock,starttime);

% check whether there has been a warning since last time
% if ~isempty(lastwarn)
%     lastwarn('');
%     fprintfExpression = regexprep(fprintfExpression,'(\\b)*','\\n');
%     fprintfExpressionFinal = regexprep(fprintfExpressionFinal,'(\\b)*','\\b\\n');
% end

if percentDone == 100 % Task completed
    fprintf(1,fprintfExpressionFinal,convertTime(runTime)); % finish up
    clear starttime lastupdate clearText printText finalText % Clear persistent vars
    return
end

% only update if significant time has passed
if etime(clock,lastupdate) < 0.3
    return
end

% update
timeLeft = runTime/fractionDone - runTime;
fprintf(1,fprintfExpression,percentDone,convertTime(timeLeft));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfcn
function timeString = convertTime(time)

timeStruct = sec2struct(time);
if timeStruct.hour > 99
    timeString = '99:59:59';
else
    timeString = timeStruct.str(1:end-4);
end
