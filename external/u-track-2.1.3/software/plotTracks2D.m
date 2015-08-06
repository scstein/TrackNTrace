function h=plotTracks2D(trackedFeatureInfo,timeRange,colorTime,markerType,...
    indicateSE,newFigure,image,flipXY,ask4sel,offset,minLength)
%PLOTTRACKS2D plots a group of tracks in 2D and allows user to click on them and extract track information
%
%SYNOPSIS plotTracks2D(trackedFeatureInfo,timeRange,colorTime,markerType,...
%    indicateSE,newFigure,image,flipXY,ask4sel,offset)
%
%INPUT  trackedFeatureInfo: -- EITHER --
%                           Output of trackWithGapClosing:
%                           Matrix indicating the positions and amplitudes
%                           of the tracked features to be plotted. Number
%                           of rows = number of tracks, while number of
%                           columns = 8*number of time points. Each row
%                           consists of
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points
%                           where the track does not exist.
%                           -- OR --
%                           Output of trackCloseGapsKalman:
%                           Structure array with number of entries equal to
%                           the number of tracks (or compound tracks when
%                           merging/splitting are considered). Contains the
%                           fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of
%                              frames the compound track spans. Each row
%                              consists of
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       timeRange         : 2-element row vector indicating time range to plot.
%                           Optional. Default: whole movie.
%       colorTime         : String with the following options:
%                           -'1' if time is to be color-coded (green in the
%                           beginning, blue in the middle, red in the end).
%                           -'k', 'b', 'r', etc. if all tracks are in black,
%                              blue, red, etc.
%                           -'2' if tracks are colored by cycling through
%                              the plot's default color order.
%                           -'3' as 2, but using extendedColors (23
%                              different colors).
%                           Optional. Default: 'k'.
%       markerType        : String indicating marker type for plotting.
%                           Only used if colorTime is not '1'.
%                           Optional. Default: 'none'.
%       indicateSE        : 1 if track starts and ends are to be indicated
%                           with circles and squares, respectively; 0
%                           otherwise. Optional. Default: 1.
%       newFigure         : 1 if plot should be made in a new figure
%                           window, 0 otherwise (in which case it will be
%                           plotted in an existing figure window).
%                           newFigure can also be a handle to the axes in
%                           which to plot.
%                           Optional. Default: 1.
%       image             : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image
%       flipXY            : 1 if x and y coord should be flipped for
%                           plotting. Optional. Default: 0.
%       ask4sel           : 1 if user should be asked to select tracks in
%                           plot in order to show track information.
%                           Optional. Default: 1.
%       offset            : [dx,dy] that is to be added to the coordinates.
%                           Optional. Default: [0,0]
%       minLength         : Minimum length of tracks to be ploted.
%                           Optional. Default: 1.
%
%OUTPUT The plot.
%
%Khuloud Jaqaman, August 2006
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--plotTracks2D: Incorrect number of input arguments!');
    return
end

%keep only tracks with minimum requested length
if nargin < 11 || isempty(minLength)
    minLength = 1;
end
if minLength > 1
    criteria.lifeTime.min = minLength;
    indx = chooseTracks(trackedFeatureInfo,criteria);
    trackedFeatureInfo = trackedFeatureInfo(indx,:);
end

%get number of tracks and number of time points
if isstruct(trackedFeatureInfo) %if tracks are in structure format
    numTracks = length(trackedFeatureInfo);
    tmp = vertcat(trackedFeatureInfo.seqOfEvents);
    numTimePoints = max(tmp(:,1));
    clear tmp
else %if tracks are in matrix format
    [numTracks,numTimePoints] = size(trackedFeatureInfo);
    numTimePoints = numTimePoints/8;
end

errFlag = 0;

%check whether a time range for plotting was input
if nargin < 2 || isempty(timeRange)
    timeRange = [1 numTimePoints];
else
    if timeRange(1) < 1 || timeRange(2) > numTimePoints
        disp('--plotTracks2D: Wrong time range for plotting!');
        errFlag = 1;
    end
end

%check whether colorTime was input
if nargin < 3 || isempty(colorTime)
    colorTime = 'k';
end
% make sure colorTime 1,2 are strings
if isnumeric(colorTime) && isscalar(colorTime)
    colorTime = num2str(colorTime);
end

%check whether markerType was input
if nargin < 4 || isempty(markerType)
    markerType = 'none';
end

%check whether indicateSE was input
if nargin < 5 || isempty(indicateSE)
    indicateSE = 1;
else
    if indicateSE ~= 0 && indicateSE ~= 1
        disp('plotTracks2D: indicateSE should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether newFigure was input
if nargin < 6 || isempty(newFigure)
    newFigure = 1;
else
    if newFigure ~= 0 && newFigure ~= 1 && ~(ishandle(newFigure) && strmatch(get(newFigure,'Type'),'axes'))
        disp('--plotTracks2D: newFigure should be 0 or 1 or an axes handle!');
        errFlag = 1;
    end
end
% check for axes handle
if ishandle(newFigure) && strcmp(get(newFigure,'Type'),'axes')
    axH = newFigure;
    newFigure = 0;
elseif newFigure == 0
    axH = gca;
end

%check whether user supplied an image
if nargin < 7 || isempty(image)
    image = [];
end

if nargin < 8 || isempty(flipXY)
    flipXY = false;
end
if nargin < 9 || isempty(ask4sel)
    ask4sel = true;
end
if nargin < 10 || isempty(offset)
    offset = [0,0];
end

%exit if there are problems in input variables
if errFlag
    disp('--plotTracks2D: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(trackedFeatureInfo) %if tracks are input in structure format

    %store the input structure as a variable with a different name
    inputStructure = trackedFeatureInfo;
    clear trackedFeatureInfo;

    %get number of segments making each track
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(inputStructure(i).tracksCoordAmpCG,1);
    end

    %if all tracks have only one segment ...
    if max(numSegments) == 1

        %indicate that there are no compound tracks with merging and splitting branches
        mergeSplit = 0;

        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        %in this case of course every compound track is simply one track
        %without branches
        trackStartRow = (1:numTracks)';

        %store tracks in a matrix
        trackedFeatureInfo = NaN*ones(numTracks,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(i,8*(startTime-1)+1:8*endTime) = inputStructure(i).tracksCoordAmpCG;
        end

    else %if some tracks have merging/splitting branches

        %indicate that in the variable mergeSplit
        mergeSplit = 1;

        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        trackStartRow = ones(numTracks,1);
        for iTrack = 2 : numTracks
            trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
        end

        %put all tracks together in a matrix
        trackedFeatureInfo = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(trackStartRow(i):trackStartRow(i)+...
                numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
                inputStructure(i).tracksCoordAmpCG;
        end

    end

else %if tracks are not input in structure format

    %indicate that there are no compound tracks with merging and splitting branches
    mergeSplit = 0;

    %indicate that each track consists of one segment
    numSegments = ones(numTracks,1);

    %locate the row of the first track of each compound track in the
    %big matrix of all tracks
    %in this case of course every compound track is simply one track
    %without branches
    trackStartRow = (1:numTracks)';

end

%get the x,y-coordinates of features in all tracks
if flipXY
    tracksY = trackedFeatureInfo(:,1:8:end)' + offset(1);
    tracksX = trackedFeatureInfo(:,2:8:end)' + offset(2);
else
    tracksX = trackedFeatureInfo(:,1:8:end)' + offset(1);
    tracksY = trackedFeatureInfo(:,2:8:end)' + offset(2);
end

% % add offset
% tracksX = tracksX ;
% tracksY = tracksY ;

%find x-coordinate limits
minXCoord = min(floor(min(tracksX(:))),0);
maxXCoord =  ceil(max(tracksX(:)));

%find y-coordinate limits
minYCoord = min(floor(min(tracksY(:))),0);
maxYCoord =  ceil(max(tracksY(:)));

%calculate the number of time points to be plotted
numTimePlot = timeRange(2) - timeRange(1) + 1;

%define colors to loop through in case colorTime = '2'
if isempty(image)
    colorLoop = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]; %colors: k,r,g,b,y,m,c
else
    colorLoop = [0.7 0.7 0.7; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]; %colors: 'ligh pink',r,g,b,y,m,c
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if the user wants to plot in a new figure window
if newFigure

    %open new figure window
    figure

    if ~isempty(image) %if user supplied an image
        tmpImage = image(image~=0);
        minImage = min(tmpImage);
        maxImage = max(tmpImage);
        imshow(image,[minImage maxImage]); %plot the image
    else %if user did not supply an image
        imshow(ones(maxYCoord,maxXCoord),[]); %plot an empty image
    end

    %set figure axes limits
    axis([minXCoord maxXCoord minYCoord maxYCoord]);

    %show coordinates on axes
    axH = gca;
    set(axH,'visible','on');

    %label axes
    if flipXY
        xlabel('y-coordinate (pixels)');
        ylabel('x-coordinate (pixels)');
    else
        xlabel('x-coordinate (pixels)');
        ylabel('y-coordinate (pixels)');
    end

end

%hold on axes
set(axH,'NextPlot','add');

%extract the portion of tracksX and tracksY that is of interest
tracksXP = tracksX(timeRange(1):timeRange(2),:);
tracksYP = tracksY(timeRange(1):timeRange(2),:);

switch colorTime

    case '1' %if user wants to color-code time

        xData=arrayfun(@(x)[tracksXP(~isnan(tracksXP(:,x)),x); NaN],1:size(tracksXP,2),'Unif',false);
        yData=arrayfun(@(x)[tracksYP(~isnan(tracksYP(:,x)),x); NaN],1:size(tracksYP,2),'Unif',false);
        plot(axH,vertcat(xData{:}),vertcat(yData{:}),'k:');
        
        %get the overall color per time interval
        colorOverTime = timeColormap(numTimePlot);

        %overlay tracks with color coding wherever a feature has been detected
        for i=1:numTimePlot-1
            validData=~all(isnan(tracksXP(i:i+1,:)),1);
            xData=vertcat(tracksXP(i:i+1,validData),NaN(1,sum(validData)));
            yData=vertcat(tracksYP(i:i+1,validData),NaN(1,sum(validData)));
            plot(axH,xData(:),yData(:),'color',colorOverTime(i,:));
        end

    case '2' %no time color-coding, loop through series of colors to color tracks

        %plot tracks by looping through colors
        %missing intervals are indicated by a dotted line
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot(axH,tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:');
            plot(axH,tracksXP(:,i),tracksYP(:,i),'color',colorLoop(mod(i-1,7)+1,:),...
                'marker',markerType);
        end
        
    case '3' % no time color-coding, use extendedColors
        
        %plot tracks by looping through colors
        %missing intervals are indicated by a dotted line
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot(axH,tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:');
            plot(axH,tracksXP(:,i),tracksYP(:,i),'color',extendedColors(mod(i-1,23)+1),...
                'marker',markerType);
        end

    otherwise %no time color-coding, all tracks same color

        %plot tracks with the line color indicated
        %missing intervals are indicated by a dotted line
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot(axH,tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:');
            plot(axH,tracksXP(:,i),tracksYP(:,i),colorTime,'marker',markerType);
            %             plot(axH,tracksXP(obsAvail,i),tracksYP(obsAvail,i),':','Color',[0.7 0.7 0.7]);
            %             plot(axH,tracksXP(:,i),tracksYP(:,i),'marker',markerType,'Color',[0.7 0.7 0.7]);
        end

end %(switch colorTime)

%show merges and splits
if mergeSplit

    %go over all tracks
    for iTrack = 1 : numTracks

        %parse sequence of events of this compound track and find merges and
        %splits
        seqOfEvents = inputStructure(iTrack).seqOfEvents;
        indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
        indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';

        %go over all splits
        for iSplit = indxSplit

            %get time of splitting
            timeSplit = seqOfEvents(iSplit,1);

            %determine row where starting track is located
            rowS = trackStartRow(iTrack) + seqOfEvents(iSplit,3) - 1;

            %determine row where splitting track is located
            rowSp = trackStartRow(iTrack) + seqOfEvents(iSplit,4) - 1;

            %plot split as a dash-dotted line
            plot(axH,[tracksX(timeSplit,rowS) tracksX(timeSplit-1,rowSp)], ...
                [tracksY(timeSplit,rowS) tracksY(timeSplit-1,rowSp)],'k-.');
            %             plot(axH,[tracksX(timeSplit,rowS) tracksX(timeSplit-1,rowSp)], ...
            %                 [tracksY(timeSplit,rowS) tracksY(timeSplit-1,rowSp)],'-.','Color',[0.7 0.7 0.7]);

        end

        %go over all merges
        for iMerge = indxMerge

            %get time of merging
            timeMerge = seqOfEvents(iMerge,1);

            %determine row where ending track is located
            rowE = trackStartRow(iTrack) + seqOfEvents(iMerge,3) - 1;

            %determine row where merging track is located
            rowM = trackStartRow(iTrack) + seqOfEvents(iMerge,4) - 1;

            %plot merge as a dashed line
            plot(axH,[tracksX(timeMerge-1,rowE) tracksX(timeMerge,rowM)], ...
                [tracksY(timeMerge-1,rowE) tracksY(timeMerge,rowM)],'k--');
            %             plot(axH,[tracksX(timeMerge-1,rowE) tracksX(timeMerge,rowM)], ...
            %                 [tracksY(timeMerge-1,rowE) tracksY(timeMerge,rowM)],'--','Color',[0.7 0.7 0.7]);

        end

    end %(for iTrack = 1 : numTracks)

end %(if mergeSplit)

if indicateSE %if user wants to indicate starts and ends

    %if there are merges and splits
    if mergeSplit

        %go over all tracks
        for iTrack = 1 : numTracks

            %parse sequence of events of this compound track and find starts and
            %ends
            seqOfEvents = inputStructure(iTrack).seqOfEvents;
            indxStart = (find(seqOfEvents(:,2) == 1 & isnan(seqOfEvents(:,4)) ...
                & seqOfEvents(:,1) >= timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
            indxEnd = (find(seqOfEvents(:,2) == 2 & isnan(seqOfEvents(:,4)) ...
                & seqOfEvents(:,1) >= timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';

            %get the information of the starts
            startInfo = [];
            for i = 1 : length(indxStart)
                iStart = indxStart(i);

                %get start time
                timeStart = seqOfEvents(iStart,1);

                %determine row where starting track is located in big matrix
                %of tracks
                rowS = trackStartRow(iTrack) + seqOfEvents(iStart,3) - 1;

                %get coordinates at the start
                startInfo(i,:) = [tracksX(timeStart,rowS) tracksY(timeStart,rowS) timeStart];

            end

            %get the information of the ends
            endInfo = [];
            for i = 1 : length(indxEnd)
                iEnd = indxEnd(i);

                %get end time
                timeEnd = seqOfEvents(iEnd,1);

                %determine row where ending track is located in big matrix
                %of tracks
                rowE = trackStartRow(iTrack) + seqOfEvents(iEnd,3) - 1;

                %get coordinates at the end
                endInfo(i,:) = [tracksX(timeEnd,rowE) tracksY(timeEnd,rowE) timeEnd];

            end

            %place circles at track starts and squares at track ends
            switch colorTime
                case {'1','2','3'}
                    if ~isempty(startInfo)
                        plot(axH,startInfo(:,1),startInfo(:,2),'k',...
                            'LineStyle','none','marker','o');
                    end
                    if ~isempty(endInfo)
                        plot(axH,endInfo(:,1),endInfo(:,2),'k',...
                            'LineStyle','none','marker','square');
                    end
                    % JD 2/09: case 2 is identical to case 1
%                 case '2'
%                     if ~isempty(startInfo)
%                         plot(axH,startInfo(:,1),startInfo(:,2),'k',...
%                             'LineStyle','none','marker','o');
%                     end
%                     if ~isempty(endInfo)
%                         plot(axH,endInfo(:,1),endInfo(:,2),'k',...
%                             'LineStyle','none','marker','square');
%                     end
                otherwise
                    if ~isempty(startInfo)
                        plot(axH,startInfo(:,1),startInfo(:,2),colorTime,...
                            'LineStyle','none','marker','o');
                    end
                    if ~isempty(endInfo)
                        plot(axH,endInfo(:,1),endInfo(:,2),colorTime,...
                            'LineStyle','none','marker','square');
                    end
            end

        end %(for iTrack = 1 : numTracks)

    else %if there are no merges and splits

        %find the beginning and end of each track
        for i=numTracks:-1:1
            timePoint = find(~isnan(tracksX(:,i)));
            startInfo(i,:) = [tracksX(timePoint(1),i) ...
                tracksY(timePoint(1),i) timePoint(1)];
            endInfo(i,:) = [tracksX(timePoint(end),i) ...
                tracksY(timePoint(end),i) timePoint(end)];
        end

        %place circles at track starts and squares at track ends if they happen to
        %be in the plotting region of interest
        switch colorTime
            case {'1','2','3'}
                indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
                plot(axH,startInfo(indx,1),startInfo(indx,2),'k','LineStyle','none','marker','o');
                indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
                plot(axH,endInfo(indx,1),endInfo(indx,2),'k','LineStyle','none','marker','square');
                %             case '2'
                %                 indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
                %                 plot(axH,startInfo(indx,1),startInfo(indx,2),'k','LineStyle','none','marker','o');
                %                 indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
                %                 plot(axH,endInfo(indx,1),endInfo(indx,2),'k','LineStyle','none','marker','square');
            otherwise
                indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
                plot(axH,startInfo(indx,1),startInfo(indx,2),colorTime,...
                    'LineStyle','none','marker','o');
                indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
                plot(axH,endInfo(indx,1),endInfo(indx,2),colorTime,...
                    'LineStyle','none','marker','square');
        end

    end %(if mergeSplit)

end %(if indicateSE)

if ask4sel

    %ask the user whether to click on figure and get frame information
    userEntry = input('select points in figure? y/n ','s');

    while strcmp(userEntry,'y')

        %let the user choose the points of interest
        [x,y] = getpts;

        %find the time points of the indicated points
        for i=1:length(x)

            %find the distances between those points and the tracks
            distTrack2Point = (tracksXP-x(i)).^2+(tracksYP-y(i)).^2;

            %determine the minimum distance for each chosen point
            [frameChosen,rowChosen] = find(distTrack2Point==min(distTrack2Point(:)));

            %go over all chosen rows
            for j = 1 : length(rowChosen)

                %find the track corresponding to each minimum distance
                trackChosen = find(trackStartRow <= rowChosen(j),1,'last');
                segmentChosen = rowChosen(j) - trackStartRow(trackChosen) + 1;

                disp(['Track: ' num2str(trackChosen) ...
                    '   Segment: ' num2str(segmentChosen) ...
                    '   Frame: ' num2str(frameChosen(j)+timeRange(1)-1) ...
                    '   Coordinates: ' num2str(tracksXP(frameChosen(j),rowChosen(j))) ...
                    ' ' num2str(tracksYP(frameChosen(j),rowChosen(j)))  ]);

            end

        end

        %ask the user again whether to click on figure and get frame information
        userEntry = input('select points again? y/n ','s');

    end

end

%%%%% ~~ the end ~~ %%%%%

