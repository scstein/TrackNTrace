% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
%     Copyright (C) 2016  Simon Christoph Stein, scstein@phys.uni-goettingen.de
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [filename_movies, globalOptions,importOptions, candidateOptions,refinementOptions,trackingOptions,postprocOptions, candidateData_loaded, refinementData_loaded, trackingData_loaded, postprocData_loaded, movieSize_loaded, firstFrame_lastFrame_loaded, metadata_loaded, outputPath_loaded, GUIreturns] = startupGUI(formats,filename_movies)
% TrackNTrace startup GUI. Here the input (movies/TNT data file) can be
% selected when RunTrackNTrace.m is called.
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2016
%
% Extended CT, 2019:
%	- support for any filetype covered by an importPlugin (in argument formats)
%	- filename_movies can be passed as input argument

GUIreturns.userExit = false;

globalOptions = [];
importOptions = [];
candidateOptions = [];
refinementOptions = [];
trackingOptions = [];
postprocOptions = [];

candidateData_loaded = [];
refinementData_loaded = [];
trackingData_loaded = [];
postprocData_loaded = [];
movieSize_loaded = [];
firstFrame_lastFrame_loaded = [];
metadata_loaded = [];
outputPath_loaded = -1;

if nargin<2 || ~iscell(filename_movies)
    filename_movies ={};
else % Take existing files from provided list
    %TODO check extension against format list
    ind = logical(cellfun(@(f)exist(f,'file'),filename_movies));
    if ~all(ind)
        cellfun(@(f)fprintf('TNT: File not found: %s\n',f),filename_movies(~ind),'UniformOutput',false);
    end
    filename_movies = filename_movies(ind);
end

% -- Preparing the GUI --
h_main = openfig('startupGUI_Layout.fig');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % For cleanup
movegui(h_main,'center');

h_all = guihandles(h_main);

% Setup GUI elements
set(h_all.button_addMovies, 'Callback', @callback_addMovies);
set(h_all.button_removeMovie, 'Callback', @callback_removeMovie);
set(h_all.button_load, 'Callback', @callback_loadSettings);
set(h_all.button_start, 'Callback',@callback_exit);


if numel(filename_movies) > 0 % show passed movies
    set(h_all.button_start,'Enable','on');
    set(h_all.listbox_movieList,'String',filename_movies);
end

% % GUI main
uiwait(h_main);
drawnow; % makes figure disappear instantly (otherwise it looks like it is existing until script finishes)

% Load settings from a file
    function callback_loadSettings(hObj,event)        
        [infile, path] = uigetfile({'*.mat','TNT settings'});
        if isfloat(infile);
            return;
        end; % User clicked cancel
        
        outputPath_loaded = path;
        
        % Note: Loading has to be done this way, as variables "can not be
        % added to a static workspace" (e.g. the one of this GUI).
        warning off
        allOptions = load([path,infile],'globalOptions','importOptions','candidateOptions','refinementOptions','trackingOptions','postprocOptions','filename_movie','candidateData','refinementData','trackingData','postprocData','movieSize','firstFrame_lastFrame');
        warning on
        globalOptions   = allOptions.globalOptions;
        if isfield(allOptions,'importOptions')
            importOptions = allOptions.importOptions;
        end
        if isfield(allOptions,'candidateOptions')
            candidateOptions = allOptions.candidateOptions;
        end
        if isfield(allOptions,'refinementOptions')
            refinementOptions   = allOptions.refinementOptions;
        end
        if isfield(allOptions,'trackingOptions')
            trackingOptions  = allOptions.trackingOptions;
        end
        if isfield(allOptions,'postprocOptions')
            postprocOptions  = allOptions.postprocOptions;
        end
        if isfield(allOptions,'candidateData')
            candidateData_loaded  = allOptions.candidateData;
        end
        if isfield(allOptions,'refinementData')
            refinementData_loaded  = allOptions.refinementData;
        end
        if isfield(allOptions,'trackingData')
            trackingData_loaded  = allOptions.trackingData;
        end
        if isfield(allOptions,'postprocData')
            postprocData_loaded  = allOptions.postprocData;
        end
        if isfield(allOptions,'movieSize')
            movieSize_loaded  = allOptions.movieSize;
        end
        if isfield(allOptions,'firstFrame_lastFrame')
            firstFrame_lastFrame_loaded  = allOptions.firstFrame_lastFrame;
        end
        if isfield(allOptions,'metadata')
            metadata_loaded  = allOptions.metadata;
        end
        
        filename_movies = {allOptions.filename_movie};
        
        callback_exit();
    end

% Opens a file chooser dialog to choose multiple input (movie) files
% for processing. Note: filenames will be seperated by ';'
    function callback_addMovies(hObj,event)
        % Get current text field to set starting path of uigetfile
        filename_movies = get(h_all.listbox_movieList,'String');
        
        % We have to do this, since repopulating a listbox does not
        % automatically reset its value..
        if numel(filename_movies) == 0
            set(h_all.listbox_movieList,'Value',1);
        end
        
        path = pwd;
        if ~isempty(filename_movies)
            [path,~,~] = fileparts(filename_movies{end});
        end
        [movieList, path] = uigetfile(formats,'Select files to load',[path,filesep],'MultiSelect','on');
        if( isfloat(movieList) ); return; end; % User pressed cancel.
        
        % atach path to every entry in the list then add to listbox
        if iscell(movieList)
            for i=1:length(movieList)
                movieList{i} =  [path,movieList{i}];
            end
            % Add to listbox
            filename_movies = [filename_movies; movieList.'];
        elseif ischar(movieList)
            movieList = [path,movieList];
            % Add to listbox
            filename_movies = [filename_movies; {movieList}];
        end
        
        set(h_all.listbox_movieList,'String',filename_movies);
        
        % Enable the Start button
        if numel(filename_movies) > 0
            set(h_all.button_start,'Enable','on');
        end
    end

% Remove a movie from the movie list
    function callback_removeMovie(hObj,event)
        selected_entry = get(h_all.listbox_movieList,'Value');
        filename_movies = get(h_all.listbox_movieList,'String');
        
        % When listbox is empty, do nothing
        if numel(filename_movies) == 0
           return; 
        end
        
        % When last selected item is deleted, select the one before it
        if selected_entry == numel(filename_movies)
            set(h_all.listbox_movieList,'Value',selected_entry-1);
        end
        
        filename_movies(selected_entry) = [];
        set(h_all.listbox_movieList,'String',filename_movies);
        
        % Disable the Start button
        if numel(filename_movies) == 0
            set(h_all.button_start,'Enable','off');
        end
    end



    function callback_exit(hObj, event)
        delete(h_main);
    end

% Called when closing the application via the 'X' button (or via close)
    function onAppClose(hObj, event)
        GUIreturns.userExit = true;
        delete(h_main);
    end

end
