h = openfig('TestGUI');
handles = guihandles(h);

% set(handles.foofoo, 'Parent', handles.panel_candidates);
% panel_pos = get(handles.panel_candidates,'Position');

% name = plugin_crossCorrelation(handles.panel_candidates, options);
name = plugin_intensityFiltering(handles.panel_candidates);
