classdef TNTscalebar < matlab.mixin.SetGet & handle
    %TNTscalebar Class implementing a scalebar object for images
    %   A TNTscalebar consists out of line and label underneath.
    %   Accepted syntax:
    %    TNTscalebar
    %    TNTscalebar pixelsize
    %    TNTscalebar pixelsize unit
    %    TNTscalebar off
    %    TNTscalebar(Name,Value)
    %    TNTscalebar(ax,__)
    %    hdl = TNTscalebar(__)
    %   Examples:
    %    TNTscalebar 0.1 ?m
    %    sb = TNTscalebar(ax, 'Pixelsize', 0.2, 'Color', 'blue', 'Location', 'SouthWest');
    %
    %   The scalebar length is calculated based on the y dimension. It is
    %   assumed that the the axis values are in pixels. To use it with a
    %   correctly scaled axis set pixelsize to 1.
    %
    %   Extra:
    %    TNTscalebar.addToogle
    %    TNTscalebar.addToogle(fig)
    %   Adds a scalebar toogle button to the current/specified figure's 
    %   toolbar.
    %    TNTscalebar.removeToogle
    %    TNTscalebar.removeToogle(fig)
    %   Removes toogle button from the current/specified figure's toolbar.
    %
    % (C) Christoph Thiele, 2020.
    
    properties
        % The location of the scalebar in the axes.
        Location char {mustBeMemberChari(Location,{'North','South','East', 'West','NorthEast','SouthEast','NorthWest','SouthWest','manual'})} = 'SouthEast';
        % VerticalAlignment of the scalebar. For auto it is set by the Location.
        VerticalAlignment {mustBeMemberChari(VerticalAlignment,{'center','top','bottom','auto'})} = 'auto';
        % HorizontalAlignment of the scalebar. For auto it is set by the Location.
        HorizontalAlignment {mustBeMemberChari(HorizontalAlignment,{'center','left','right','auto'})} = 'auto';
        % Switch the Length of scalebar is set/updated automatically
        LengthMode {mustBeMemberChari(LengthMode,{'manual','auto'})} = 'auto';
        % Pixelsize in the unit of the scalebar
        Pixelsize(1,1) double {mustBePositive} = 1;
        % Accepted scalebar length values
        Numerals(:,1) double {mustBePositive} = [1 2 5];
        % Minimal length of the scalebar relative to the y axis
        MinLength(1,1) double {mustBePositive} = 0.2;
        % Unit of the scalebar
        Unit char = '';
        % Formating of the scalebar label (formatSpec for sprintf). 
        NumberFormat char = '%g';
        % (true/false) Ensures that the scalebar is the topmost element
        KeepOnTop matlab.lang.OnOffSwitchState = matlab.lang.OnOffSwitchState.on;
    end
    properties (Dependent)
        % FontSize of the scalebar label
        FontSize(1,1) double {mustBePositive}
        % FontName of the scalebar label
        FontName
        % FontWeight of the scalebar label
        FontWeight
        % LineWidth of the scalebar line
        LineWidth(1,1) double {mustBePositive}
        % Color of the scalebar line and label
        Color
        % Scalebar position in normalised coordinates
        Position
        % Visiblity of the scalebar
        Visible matlab.lang.OnOffSwitchState
        % Length of the scalebar in the scalebar unit
        Length(1,1) double
        % (true/false) Automatically update the scalebar length and position on pan/zoom
        AutoUpdate matlab.lang.OnOffSwitchState
    end % properties (Dependent)
    properties ( SetAccess = private )
        % Parent axes
        Parent
        Line
        Text
    end % properties ( Access = private )
    properties ( SetAccess = private, GetAccess = private )
        % Parent axes
        Length_
        Position_
        UpdateListener
        DeleteListener
    end % properties ( Access = private, GetAccess = private )
    
    methods
        function obj = TNTscalebar(varargin)
            %TNTscalebar Construct an instance of this class
            %   Detailed explanation goes here
            
            % Get the parent axis
            cax = gobjects(0);
            if nargin>0
                if nargin == 1 && strcmpi(varargin{1},'off')
                    % Parse: TNTscalebar off
                    varargin = {'Visible','off'};
                elseif nargin>0 && nargin<3 && iscellstr(varargin) %#ok<ISCLSTR>
                    % Parese: TNTscalebar pixelsize unit
                    if ~isnan(str2double(varargin{1}))
                        if nargin==2
                            varargin = {'Pixelsize',str2double(varargin{1}),'Unit',[' ' varargin{2}]};
                        else
                            varargin = {'Pixelsize',str2double(varargin{1})};
                        end
                        varargin = [varargin,{'Visible','on'}];
                    end
                end
                if isgraphics(varargin{1})
                    % use ancestor to make sure it also works for images as
                    % input.
                    cax = ancestor(varargin{1},'axes');
                    varargin = varargin(2:end);
                else
                    axind = strcmpi('Parent',varargin);
                    if any(axind) && find(axind,1,'last')<numel(varargin)
                        cax = ancestor(varargin{find(axind,1,'last')+1},'axes');
                        varargin = varargin(~axind & circshift(~axind,1));
                    end
                end
            else
                varargin = {'Visible','on'};
            end
            if isempty(cax)
                cax = gca;
            end
            
            
            % Check for existing scalebar
            oldobj = getappdata(cax,'TNTscalebar');
            if isempty(oldobj)
                obj.Parent = cax;
            
                % Create objects
                obj.Line = line(nan,nan,'Parent',obj.Parent,'LineStyle','-','LineWidth',5,'Color',[1 1 1],'Tag','TNTscalebar_Line','PickableParts','none','HandleVisibility','off');
                obj.Line.Annotation.LegendInformation.IconDisplayStyle;
                obj.Text = text(nan,nan,'','Parent',obj.Parent,'FontSize',14,'Color',obj.Line.Color,'Units','data','HorizontalAlignment','center','VerticalAlignment','top','Tag','TNTscalebar_Text','HandleVisibility','off');
                
                % Save in appdata
                setappdata(obj.Parent,'TNTscalebar',obj);
                % Delete the object if the scalebar is deleted (e.g. by a new plot)
                obj.DeleteListener = listener(obj.Line,'ObjectBeingDestroyed', @obj.delete_callback);
                %listener(obj.Text,'ObjectBeingDestroyed', @obj.delete);
                obj.UpdateListener = listener(obj.Parent,{'XLim','YLim','Position'},'PostSet',@obj.autoupdate_callback);
                % Note: This listener does not cover all cases. Setting
                % xlim/ylim to automatic or restoring the initial view does
                % not trigger the PostSet callback. Here a listener to the
                % XRuler/YRuler MarkedClean event could help.
            else
                % This destroys the just created new object and uses the
                % existing one.
                obj = oldobj;
            end
            % Set properties
            obj.setProperties(varargin{:});
            % 
            obj.updateLength();
            obj.updateText();
            obj.updatePosition();
            obj.updatePosition();
            syncToogle(ancestor(obj.Parent,'figure'));
            
            if nargout == 0
                % Do not return the object if no output is requested.
                clear obj;
            end
        end
        function delete(obj)
            %DELETE Destructor method.
            delete(obj.UpdateListener);
            delete(obj.DeleteListener);
            
            delete(obj.Line);
            delete(obj.Text);
            % Remove from appdata if the axes still exist
            if ~isempty(obj.Parent) && isgraphics(obj.Parent) && isappdata(obj.Parent,'TNTscalebar')
                rmappdata(obj.Parent,'TNTscalebar');
            end
            
        end % destructor
        function update(obj)
            if ~isempty(obj.Parent)
                obj.updateLength();
                obj.updateText();
                obj.updatePosition();
            end
        end
        %% Set/Get
        % Properties (get only)
        function value = get.Parent( obj)
            value = obj.Parent;
        end % get.Parent
        % Properties (set only)
        function set.Location( obj, value)
            obj.Location = value;
            obj.updatePosition();
        end % set.Location
        function set.VerticalAlignment( obj, value)
            obj.VerticalAlignment = value;
            obj.updatePosition();
        end % set.VerticalAlignment
        function set.HorizontalAlignment( obj, value)
            obj.HorizontalAlignment = value;
            obj.updatePosition();
        end % set.HorizontalAlignment
        function set.LengthMode( obj, value)
            obj.LengthMode = value;
            if strcmpi(obj.LengthMode,'auto')
                obj.update();
            end
        end % set.LengthMode
        function set.Pixelsize( obj, value)
            obj.Pixelsize = value;
            obj.update();
        end % set.Pixelsize
        function set.Numerals( obj, value)
            obj.Numerals = value;
            obj.update();
        end % set.Numerals
        function set.MinLength( obj, value)
            obj.MinLength = value;
            obj.update();
        end % set.MinLength
        function set.Unit( obj, value)
            obj.Unit = value;
            obj.update();
        end % set.Unit
        function set.NumberFormat( obj, value)
            obj.NumberFormat = value;
            obj.update();
        end % set.NumberFormat
        % Dependent properties
        function value = get.FontSize( obj )
            value = obj.Text.FontSize;
        end % get.FontSize
        function set.FontSize( obj, value)
            obj.Text.FontSize = value;
            obj.updatePosition();
        end % set.FontSize
        function value = get.FontName( obj )
            value = obj.Text.FontName;
        end % get.FontName
        function set.FontName( obj, value)
            obj.Text.FontName = value;
        end % set.FontName
        function value = get.FontWeight( obj )
            value = obj.Text.FontWeight;
        end % get.FontWeight
        function set.FontWeight( obj, value)
            obj.Text.FontWeight = value;
        end % set.FontWeight
        function value = get.LineWidth( obj )
            value = obj.Line.LineWidth;
        end % get.LineWidth
        function set.LineWidth( obj, value)
            obj.Line.LineWidth = value;
            obj.updatePosition();
        end % set.LineWidth
        function value = get.Color( obj )
            value = obj.Line.Color;
        end % get.Color
        function set.Color( obj, value)
            obj.Line.Color = value;
            obj.Text.Color = value;
        end % set.Color
        function value = get.Position( obj )
            value = obj.calculatePosition();
        end % get.Position
        function set.Position( obj, value)
            obj.Position_ = value;
            % Setting the Location property triggers a position update
            obj.Location = 'manual';
        end % set.Position
        function value = get.Visible( obj )
            value = obj.Line.Visible;
        end % get.Visible
        function set.Visible( obj, value)
            obj.Line.Visible = value;
            obj.Text.Visible = value;
            obj.syncToogle();
        end % set.Visible
        function set.Length( obj, value)
            if isempty(value)
                obj.LengthMode = 'auto';
                obj.updateLength();
            else
                obj.Length_ = value;
                obj.LengthMode = 'manual';
            end
        end % set.Length
        function value = get.Length( obj )
            value = obj.Length_;
        end % get.Length
        function set.AutoUpdate( obj, value)
            obj.UpdateListener.Enabled = logical(value);
        end % set.AutoUpdate
        function value = get.AutoUpdate( obj )
            value = matlab.lang.OnOffSwitchState(obj.UpdateListener.Enabled);
        end % get.AutoUpdate
    end % methods
    methods ( Access = protected )
        function setProperties( obj, varargin )
            %SETPROPERTIES Set user-specified chart properties safely.
            
            if ~isempty( varargin )
                try
                    set( obj, varargin{:} )
                catch e
                    % Destroy the chart object and throw the exception
                    % from the calling function.
                    obj.delete()
                    e.throwAsCaller()
                end % try/catch
            end % if
            
        end % setProperties
        
        function updatePosition(obj)
            if ~isempty(obj.Parent) %obj.Parent can be empty when loading a saved figure
                xrange = xlim(obj.Parent);
                yrange = ylim(obj.Parent);
                
                pos_norm = obj.calculatePosition();
                sb_position = [...
                    pos_norm(1)*diff(xrange)+xrange(1),...
                    pos_norm(2)*diff(yrange)+yrange(1)];
                
                [height,lineheight] = obj.getHeight();
                length = obj.Length_/obj.Pixelsize;
                width = max(length,obj.Text.Extent(3));
                
                line_y = [1 1].*sb_position(2)-height*pos_norm(4)+lineheight/2;
                line_x = sb_position(1)+length*[-0.5 0.5]-width*(pos_norm(3)-0.5);
                text_y = sb_position(2)-height*pos_norm(4)+lineheight;
                text_x = mean(line_x);
                
                set(obj.Line,'XData',line_x,'YData',line_y);
                set(obj.Text,'Position',[text_x,text_y,0]);
                
                if obj.KeepOnTop
                    uistack([obj.Line, obj.Text],'top');
                end
            end
        end % updatePosition
        
        function updateText(obj)
            obj.Text.String = sprintf([obj.NumberFormat obj.Unit],obj.Length_);
        end % updateText
        function updateLength(obj)
            if strcmpi(obj.LengthMode,'auto')
                xrange = xlim(obj.Parent);
                numerals = [obj.Numerals;10*obj.Numerals(1)];
                xlength = abs(xrange(2)-xrange(1))*obj.Pixelsize;
                sb_order = 10^floor(log10(xlength*obj.MinLength));
                obj.Length_ = sb_order*numerals(find(ceil(xlength*obj.MinLength/sb_order)<=numerals,1));
            end
        end % updateLength
        
        function pos = calculatePosition(obj)
            switch lower(obj.Location)
                case 'southeast'
                    pos = [0.95 0.95 1.0 1.0];
                case 'north'
                    pos = [0.50 0.05 0.5 0.0];
                case 'south'
                    pos = [0.50 0.95 0.5 1.0];
                case 'east'
                    pos = [0.95 0.50 1.0 0.5];
                case 'west'
                    pos = [0.05 0.50 0.0 0.5];
                case 'northwest'
                    pos = [0.05 0.05 0.0 0.0];
                case 'northeast'
                    pos = [0.95 0.05 1.0 0.0];
                case 'southwest'
                    pos = [0.05 0.95 0.0 1.0];
                case 'manual'
                    pos = obj.Position_;
                    pos(end+1:4) = 0.5; % Make sure that we have 4 elements
                otherwise
                    error('Unknown location property.');
            end
            switch lower(obj.HorizontalAlignment)
                case 'center'
                    pos(3) = 0.5;
                case 'left'
                    pos(3) = 0;
                case 'right'
                    pos(3) = 1;
            end
            switch lower(obj.VerticalAlignment)
                case 'center'
                    pos(4) = 0.5;
                case 'top'
                    pos(4) = 0;
                case 'bottom'
                    pos(4) = 1;
            end
        end % getLocation
        function [height,lineheight] = getHeight(obj)
            ppd = obj.getPointPerData();
            lineheight = ppd(2)*obj.Line.LineWidth;
            upper = min(obj.Line.YData)-lineheight/2;
            lower = obj.Text.Position(2)+obj.Text.Extent(4);
            height = lower-upper;
            if isnan(height)
                % Guess based on font size and linewidth
                height = (obj.Text.FontSize+obj.Line.LineWidth)*ppd(2);
            end
        end % getHeight
        
        function ppd = getPointPerData(obj)
            unit = obj.Parent.Units;
            obj.Parent.Units = 'Points';
            pos = obj.Parent.Position;
            obj.Parent.Units = unit;
            
            xrange = diff(xlim(obj.Parent));
            yrange = diff(ylim(obj.Parent));
            ppd = [xrange,yrange]./pos(3:4);
        end
        
        function autoupdate_callback(obj,~,src)
            if isobject(obj)
                if isequal(src.AffectedObject,obj.Parent)
                    obj.update();
                else
                    obj.delete();
                end
            else
                warning('Callback should be deleted.');
            end
        end
        function delete_callback(obj,~,~)
            if isobject(obj)
                obj.delete();
            else
                warning('Callback should be deleted.');
            end
        end
        
        function syncToogle(obj)
            % Sync with toogle in figure toolbar if exists.
            fig = ancestor(obj.Parent,'figure');
            toogle = findall(fig,'Tag','Annotation.InsertScalebar');
            if ~isempty(toogle) && isgraphics(toogle)
                if isequal(obj.Parent,fig.CurrentAxes)
                    toogle.State = obj.Line.Visible;
                end
            end
        end
    end % methods ( Access = protected )
    
    methods(Static)
        function toogle = addToogle(fig)
            if nargin < 1 || isempty(fig) || ~isgraphics(fig)
                fig = gcf;
            else
                fig = ancestor(fig,'figure');
            end
            toogle = findall(fig,'Tag','Annotation.InsertScalebar');
            if isempty(toogle)
                tbar = findall( fig, 'tag', 'FigureToolBar');
                if ~isempty(tbar) && isgraphics(tbar)
                    toogle = uitoggletool(tbar,'CData',scalebar_icon(),'Tooltip','Insert scalebar','Tag','Annotation.InsertScalebar','ClickedCallback',@callback_scalebartoggle);
                    addlistener(fig,'CurrentAxes','PostSet',@callback_changeCurrentAxes);
                    syncToogle(fig);
                end
            end
            if nargout == 0
                clear toogle;
            end
        end
        function removeToogle(fig)
            if nargin < 1 || isempty(fig) || ~isgraphics(fig)
                fig = gcf;
            else
                fig = ancestor(fig,'figure');
            end
            toogle = findall(fig,'Tag','Annotation.InsertScalebar');
            if ~isempty(toogle)
                delete(toogle);
            end
        end
    end
end

function mustBeMemberChari(value,list)
    mustBeMember(lower(value),lower(list));
end

function syncToogle(fig)
    % Synchronise the toogle with state of the scalebar in the current axes
    cax = get(fig,'CurrentAxes');
    toogle = findall(fig,'Tag','Annotation.InsertScalebar');    
    if ~isempty(toogle)
        % Find scalebar, if any.
        if ~isempty(cax) && isgraphics(cax) && isappdata(cax,'TNTscalebar')
            sb = getappdata(cax,'TNTscalebar');
            toogle.State = sb.Visible;
        else
            toogle.State = 'off';
        end
    end
end

function callback_scalebartoggle(obj,evt)
    cax = get(ancestor(obj,'figure'),'CurrentAxes');
    if isgraphics(cax)
        TNTscalebar(cax,'Visible',obj.State);
    else
        obj.State = 'off';
    end
end

function callback_changeCurrentAxes(obj,evt)
    syncToogle(evt.AffectedObject);
end

function icon = scalebar_icon()
    icon = [...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 NaN;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0.5;...
        0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0.5;...
        0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0.5;...
        0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0.5;...
        0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0.5;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5;...
        NaN 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
    icon = repmat(icon,1,1,3);
end

