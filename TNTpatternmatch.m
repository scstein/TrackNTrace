% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
%     Copyright (C) 2024
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
function TNTpatternmatch()

addPathsPM()

% State variables
pm_align_peaks = false(1);
pm_tg_start = 0;
pm_tg_end = inf;
pm_tcspc = [];
pm_ids = [];

n_pattern = 4;
p_enable = true(n_pattern,1);
p_file = cell(n_pattern,1);
p_filter = cell(n_pattern,1);
p_events = true(n_pattern,1);

p_tcspc = cell(n_pattern,1);
p_change = true(n_pattern,1);

data_path = pwd;
settings_path = pwd;

% -- Preparing the GUI --
h_main = openfig('TNTpatternmatch.fig');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % For cleanup
movegui(h_main,'center');

h_all = guihandles(h_main);


% Setup GUI elements
set(h_all.button_load, 'Callback', @callback_loadSettings);
set(h_all.button_save, 'Callback', @callback_saveSettings);
set(h_all.button_preview, 'Callback', @callback_preview);
set(h_all.button_apply, 'Callback', @callback_apply);
set(h_all.button_p1_file, 'Callback', @callback_selectFile);
set(h_all.button_p2_file, 'Callback', @callback_selectFile);
set(h_all.button_p3_file, 'Callback', @callback_selectFile);
set(h_all.button_p4_file, 'Callback', @callback_selectFile);
set(h_all.edit_tg_start, 'Callback', {@callback_IntEdit,0,inf});
set(h_all.edit_tg_end, 'Callback', {@callback_IntEdit,0,inf});
initgui()


% Functions to keeep gui and variables in sync
    function initgui()
        xlabel(h_all.axes_preview, 'time bin');
        ylabel(h_all.axes_preview, 'counts');
    end

    function gui2state()
        pm_align_peaks = logical(h_all.cb_align_peaks.Value);
        pm_tg_start = str2num(h_all.edit_tg_start.String);
        pm_tg_end = str2num(h_all.edit_tg_end.String);
        
        for pid = 1:n_pattern
            p_enable(pid) = logical(h_all.(['cb_p' num2str(pid) '_enable']).Value);
            % check for change
            new_file = h_all.(['edit_p' num2str(pid) '_file']).UserData;
            new_events = logical(h_all.(['radio_p' num2str(pid) '_filt_events']).Value);
            new_filter = h_all.(['edit_p' num2str(pid) '_filter']).String;
            p_change(pid) = ~isequaln({new_file,new_events,new_filter},{p_file{pid},p_events(pid),p_filter{pid}});
            if p_change(pid)
                p_file{pid} = new_file;
                p_events(pid) = new_events;
                p_filter{pid} = new_filter;
            end
        end
    end

    function state2gui()
        h_all.cb_align_peaks.Value = pm_align_peaks;
        h_all.edit_tg_start.String = num2str(pm_tg_start);
        h_all.edit_tg_end.String = num2str(pm_tg_end);

        for pid = 1:n_pattern
            h_all.(['cb_p' num2str(pid) '_enable']).Value = p_enable(pid);
            h_all.(['radio_p' num2str(pid) '_filt_events']).Value = p_events(pid);
            h_all.(['radio_p' num2str(pid) '_all_photon']).Value = ~p_events(pid);
            h_all.(['edit_p' num2str(pid) '_filter']).String = p_filter{pid};
            callback_selectFile(h_all.(['edit_p' num2str(pid) '_file']),[],p_file{pid})
        end
    end

% Load settings from a file
    function callback_loadSettings(hObj,event)        
        [infile, path] = uigetfile({'*.mat','TNT patterns'},settings_path);
        if isfloat(infile)
            return;
        end % User clicked cancel
        
        settings_path = path;
        
        % Note: Loading has to be done this way, as variables "can not be
        % added to a static workspace" (e.g. the one of this GUI).
        warning off
        allSettings = load([path,infile],'pm_align_peaks','pm_tg_start','pm_tg_end','pm_tcspc','pm_ids','n_pattern','p_enable','p_file','p_filter','p_tcspc','p_events','p_change','data_path');
        warning on
             
        pm_align_peaks = allSettings.pm_align_peaks;
        pm_tg_start = allSettings.pm_tg_start;
        pm_tg_end = allSettings.pm_tg_end;
        pm_tcspc = allSettings.pm_tcspc;
        pm_ids = allSettings.pm_ids;

        n_pattern = allSettings.n_pattern;
        p_enable = allSettings.p_enable;
        p_file = allSettings.p_file;
        p_filter = allSettings.p_filter;
        p_events = allSettings.p_events;
        p_tcspc = allSettings.p_tcspc;
        p_change = allSettings.p_change;

        data_path = allSettings.data_path;
        
        state2gui()
        if ~isempty(pm_tcspc)
            plotPatterns()
        end
    end

    function callback_saveSettings(hObj,event)        
        [infile, path] = uiputfile({'*.mat','TNT patterns'}, settings_path);
        if isfloat(infile)
            return;
        end % User clicked cancel
        
        settings_path = path;
        
        gui2state()
        save([path,infile],'pm_align_peaks','pm_tg_start','pm_tg_end','pm_tcspc','pm_ids','n_pattern','p_enable','p_file','p_filter','p_tcspc','p_events','p_change','data_path','-v7.3')
    end
    
    function pid = obj2pid(hObj)
        pid = 0;
        if isgraphics(hObj) && ~isempty(hObj.Tag)
            tag_id = regexpi(hObj.Tag,'_p(\d+)_','tokens','once');
            if ~isempty(tag_id)
                pid = str2num(tag_id{1});
            end
        end
    end

    function callback_selectFile(hObj,~,file)
        pid = obj2pid(hObj);
        
        if nargin<3
            newPath = h_all.(['edit_p' num2str(pid) '_file']).UserData;
            if isempty(newPath)
                newPath = data_path;
            else
                [newPath,~,~] = fileparts(data_path);
            end
            [newMovie, newPath] = uigetfile({'*.mat','Visualizer files';'*.*','All files'},'Select movie or TNT file as pattern',[newPath,filesep]);

            if isfloat(newMovie)
                return
            end
            data_path = newPath;
            tntFile = [newPath filesep newMovie];
        elseif isempty(file)
            % reset
            empty_str = 'Click select to choose a file.';
            set(h_all.(['edit_p' num2str(pid) '_file']),'UserData',[]);
            set(h_all.(['edit_p' num2str(pid) '_file']),'String',empty_str);
            set(h_all.(['edit_p' num2str(pid) '_file']),'Tooltip','');
            return
        else
            tntFile = file;
        end

        [~,filename,ext] = fileparts(tntFile);
        filename = [filename,ext];
        set(h_all.(['edit_p' num2str(pid) '_file']),'UserData',tntFile);
        set(h_all.(['edit_p' num2str(pid) '_file']),'String',filename);
        set(h_all.(['edit_p' num2str(pid) '_file']),'Tooltip',tntFile);

    end

    function callback_preview(~,~)
        enabledUI = findall(h_main,'Enable','on','Type','UIControl','-not','Style','text');
        set(enabledUI,'Enable','off');
        drawnow();
        try
            gui2state();
            generatePattern();   
            plotPatterns();
        catch err
            disp( getReport( err, 'extended', 'hyperlinks', 'on' ) )
        end
        set(enabledUI,'Enable','on');
    end
    
    function callback_apply(~,~)
        [infile, path] = uigetfile({'*.mat','TNT results'},data_path,'MultiSelect','on');
        if isfloat(infile)
            return;
        end % User clicked cancel
        data_path = path;
        infile = cellstr(infile);
        
        % Generate the patterns
        enabledUI = findall(h_main,'Enable','on','Type','UIControl','-not','Style','text');
        set(enabledUI,'Enable','off');
        drawnow();
        gui2state();
        generatePattern();   
        plotPatterns();
        drawnow();
        % Apply to each file
        outfiles = {};
        for fn = infile(:)
            try
                outfiles{end+1} = applyPM([path filesep fn{1}]);
            catch
                fprintf('TNT: Error proccessing file: %s.\n', fn{1});
            end
        end
        set(enabledUI,'Enable','on');
        fprintf('######\nTNT: Processed %i files:\n',numel(outfiles));
        for iMovie=1:numel(outfiles)
            [~,tnt_fname] = fileparts(outfiles{iMovie});
            fprintf('     <a href="matlab:TNTvisualizer(''%s'')">%s</a>\n',outfiles{iMovie},tnt_fname);
        end
    end

    function callback_exit(hObj, event)
        delete(h_main);
    end

% Called when closing the application via the 'X' button (or via close)
    function onAppClose(hObj, event)
        delete(h_main);
    end

    function generatePattern()
        tcspc_maxpos = [];
        tcspc_len = [];
        pn = 1;
        p_valid = false(size(p_enable));
        for pid = 1:n_pattern
            if p_enable(pid) && ~isempty(p_file{pid})
                if p_change(pid)
                    % Only update on change
                    tcspc = getPattern(p_file{pid}, p_events(pid), p_filter{pid});
                    if ~isempty(tcspc)
                        % using previously loaded pattern if reading failed
                        p_tcspc{pid} = tcspc;
                        p_change(pid) = false;
                    end
                end
                p_valid(pid) = ~isempty(p_tcspc{pid});
                if p_valid(pid)
                    [~,tcspc_maxpos(pn)] = max(p_tcspc{pid});
                    tcspc_len(pn) = length(p_tcspc{pid});
                    pn = pn+1;
                end
            end
        end
        pm_ids = find(p_valid);
        if pn==1
            fprintf('TNT: No patterns found.\n');            
            return
        end
        if pm_align_peaks
            maxpos = floor(mean(tcspc_maxpos));
            for pid = 1:(pn-1)
                p_tcspc{pid} = circshift(p_tcspc{pid}, maxpos-tcspc_maxpos(pid));
            end
        end
        tcspc_end = min(pm_tg_end,min(tcspc_len));
        tcspc_start = max(1,pm_tg_start);
        pm_tcspc = zeros(tcspc_end-tcspc_start+1,pn-1);
        for pid = 1:(pn-1)
            pm_tcspc(:,pid) = p_tcspc{pm_ids(pid)}(1,tcspc_start:tcspc_end);
        end
        pm_tcspc(pm_tcspc==0|isnan(pm_tcspc)) = 0.01;
        pm_tcspc = pm_tcspc./sum(pm_tcspc,1);
    end

    function plotPatterns()
        if ~isempty(pm_tcspc)
            ax = h_all.axes_preview;
            semilogy(ax, (1:size(pm_tcspc,1))+(-1+max(pm_tg_start,1)), pm_tcspc);
            legend(ax, arrayfun(@(n)sprintf('Pattern %i',n),pm_ids,'UniformOutput',false));
            xlabel(ax, 'time bin');
            ylabel(ax, 'probability');
            xlim(ax,[pm_tg_start,pm_tg_start+size(pm_tcspc,1)]);
        else
            cla(h_all.axes_preview);
        end
    end

    function newfile = applyPM(file)
        newfile = [];
        if ~endsWith(file,'.mat')
            fprintf('TNT: File not a TNT results file: %s.\n', file);
            return
        end
        tntres = load(file);
        if ~isfield(tntres,'postprocData')
            fprintf('TNT: No post-processing data found in file %s.\n', file);
            return
        elseif ~isfield(tntres.postprocOptions,'TCSPC')
            fprintf('TNT: No TCSPC data in file %s.\n', file);
            return
        end
        % align peak
        if pm_align_peaks
            [~,maxpos_mix] = max(sum(tntres.postprocOptions.TCSPC,1));
            [~,maxpos_ref] = max(sum(pm_tcspc,2));
            tcspc_mix = circshift(tntres.postprocOptions.TCSPC,maxpos_ref-maxpos_mix,2);
        else
            tcspc_mix = tntres.postprocOptions.TCSPC;
        end
        % calculate matrix
        tcspc_start = max(1,pm_tg_start);
        tcspc_len = min([size(pm_tcspc,1) (size(tcspc_mix,2)-tcspc_start+1) (pm_tg_end-pm_tg_start+1)]);
        tcspc_end = tcspc_start+tcspc_len-1;
        
        ln_p_ref = log(pm_tcspc(1:tcspc_len,:));
        QMLE = tcspc_mix(:,tcspc_start:tcspc_end) * ln_p_ref;
        [QMLE_max,QMLE_ind] = max(QMLE,[],2);
        % Calculate the posterior proberbility for the nth species
        % fQMLE = exp(QMLE)./(sum(exp(QMLE),2)); % P(S1)/(P(S1)+P(S2))
        fQMLE = exp(QMLE-mean(QMLE,2))./(sum(exp(QMLE-mean(QMLE,2)),2)); % Avoids nans due to overflowing floats
        
        newParamDescription = [{'pm-id','pm-pmax'} arrayfun(@(n)sprintf('pm-p%i',n),pm_ids,'UniformOutput',false)];
        newCols = [pm_ids(QMLE_ind) max(fQMLE,[],2) fQMLE];
        if max(tntres.postprocData(:,1))==size(tcspc_mix,1)
            % one TCSPC per track
            label2struct = @(fnames)cell2struct(num2cell(1:numel(fnames))',matlab.lang.makeValidName(fnames(:)));
            postind = label2struct(tntres.postprocOptions.outParamDescription);
            newCols = newCols(tntres.postprocData(:,postind.Track_ID),:);
        end
        if size(tntres.postprocData,1)~=size(newCols,1)
            fprintf('TNT: Cannot match TCSPC with loc data in file: %s.\n', file);
        end
        tntres.postprocData = [tntres.postprocData newCols];
        tntres.postprocOptions.outParamDescription = [tntres.postprocOptions.outParamDescription; newParamDescription'];
        % save results
        newfile = [file(1:end-4) 'pm.mat'];
        try
            save(newfile,'-struct','tntres','-v7.3');
        catch
            fprintf('TNT: Error saving file: %s.\n', newfile);
            
        end
    end
end

%% --- PM functions ---

function tcspc = getPattern(file, use_events, filter_str)
    tcspc = [];
    if use_events
        if ~endsWith(file,'.mat')
            fprintf('TNT: File not a TNT results file: %s.\n', file);
            return
        end
        tntres = load(file,'postprocOptions','postprocData');
        if ~isfield(tntres,'postprocData')
            fprintf('TNT: No post-processing data found in file %s.\n', file);
            return
        elseif ~isfield(tntres.postprocOptions,'TCSPC')
            fprintf('TNT: No TCSPC data in file %s.\n', file);
            return
        end
        tcspc = getTCSPCfromTNT(tntres,filter_str);
    else
        if endsWith(file,'.mat')
            tntres = load(file,'filename_movie');
            file = tntres.filename_movie;
        end
        tcspc = getTCSPCfromRaw(file);
    end
end

function tcspcdata = getTCSPCfromTNT(tntres, filter_str)
    tcspcdata = [];
    try
        getFilterFun = @(paramsNames, filter_str)str2func(['@(posData)true(size(posData,1),1)&' regexprep(vectorize(['( ' filter_str ' )']),...
            strcat('(?<!\w)',matlab.lang.makeValidName(paramsNames),'(?!\w)'),... % Replace spaces with underscore and makes sure the match is not within a name.
            cellfun(@(n)sprintf('posData(:,%i)',n),num2cell(1:numel(paramsNames)),'UniformOutput',false),...
            'ignorecase')]);
        label2struct = @(fnames)cell2struct(num2cell(1:numel(fnames))',matlab.lang.makeValidName(fnames(:)));

        postind = label2struct(tntres.postprocOptions.outParamDescription);
        
        if isempty(filter_str)
            filter_TrackIDs = unique(tntres.postprocData(:,postind.Track_ID));
        else
            filter_fun = getFilterFun(tntres.postprocOptions.outParamDescription,filter_str);
            filter_ind = filter_fun(tntres.postprocData);
            filter_TrackIDs = unique(tntres.postprocData(filter_ind,postind.Track_ID));
        end
        tcspcdata = sum(tntres.postprocOptions.TCSPC(filter_TrackIDs,:),1);
    catch
            fprintf('TNT: Error applying filter (%s).\n', filter_str);        
    end
end

function tcspcdata = getTCSPCfromRaw(file)
    tcspcdata = [];
    % fund a suitable import plugin
    [importPlugins,importFormats] = loadImportPlugins();
    warnstate = warning('backtrace');
    warning off backtrace;
    importPlugin = selectImportPlugin(file,importPlugins);
    warning(warnstate);
    if importPlugin>0
        importPlugin = importPlugins(importPlugin);
    else
        fprintf('TNT: No import plugin for file %s.\n', file);
        return
    end
    
    % Use plugin to read TCSPC with the default options
    importOptions = importPlugin.getOptions();
    if isfield(importOptions.info,'getTCSPC') && isa(importOptions.info.getTCSPC,'function_handle')
        maxPhotons = 1e8;
        if nargin(importOptions.info.getTCSPC) == 3
            [tcspcdata,resolution] = importOptions.info.getTCSPC(file,maxPhotons,importOptions); % Resolution in s
        else
            [tcspcdata,resolution] = importOptions.info.getTCSPC(file,maxPhotons); % Resolution in s
        end
        tcspcdata = accumarray(tcspcdata,1)';
    else
        fprintf('TNT: The import plugin %s does not support preview of the TCSPC.\n', importOptions.plugin_name)
    end
end

%% --- General functions ---
% Adds pathes needed for the visualizer.
function addPathsPM()
    fullPathToThisFile = mfilename('fullpath');
    [path,~,~] = fileparts(fullPathToThisFile);
    addpath(genpath([path,filesep,'subfun']));
    addpath(genpath([path,filesep,'external']));
    addpath(genpath([path,filesep,'helper']));
    addpath(genpath([path,filesep,'plugins']));
end

% Callback for edit fields containing integer values. Checks if a correct
% number was entered and restricts it to the given bounds.
function callback_IntEdit(hObj,~, minVal,maxVal)
    if nargin<3 || isempty(minVal)
        minVal=0;
    end
    if nargin<4 || isempty(maxVal)
        maxVal=inf;
    end

    % Accept end as inf, as MATLAB users are used to end as the last element
    if(strcmp(get(hObj,'String'), 'end'))
        value = inf;
    else
        value = round(str2num(get(hObj,'String'))); %#ok<ST2NM> str2num distinguishes NaN and invalid input, str2double does not.
    end

    if isempty(value)
        set(hObj,'ForegroundColor','r');
        set(hObj,'String','INVALID');
        uicontrol(hObj);
    else
        value = max(minVal,value);
        value = min(maxVal,value);
        set(hObj,'ForegroundColor','k');
        set(hObj,'String',sprintf('%i',value));
    end
end

