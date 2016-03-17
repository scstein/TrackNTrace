function [plugin] = plugin_default()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Default name';

% Type of plugin (integer)
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = INSERT_CORRECT_NUMBER;

% The functions this plugin implements
mainFunc =  @arbitraryName_main;

% Description of output parameters
outParamDescription = {'OutputVariableName1';'OutputVariableName2'};

% Create the plugin
plugin = TNTplugin(name,type, mainFunc,outParamDescription);

% Add init or postprocessing function
plugin.initFunc = @arbitraryName_init;
plugin.postFunc = @arbitraryName_post;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Simple description comes here.';

% Activate/Deactivate this plugin's parallel processing
% deactivation necessary if plugin cannot operate on frame-by-frame basis
% will be true by default if not provided
plugin.useParallelProcessing = true;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('anIntegerVarname',...
    'int',...
    {1, 0, inf},...
    'Tooltip description.');
plugin.add_param('aDoubleVarname',...
    'float',...
    {0, 6.5, inf},...
    'Tooltip description.');
plugin.add_param('aBooleanVarname',...
    'bool',...
    false,...
    'Tooltip description.');
plugin.add_param('aListVarname',...
    'list',...
    {'Default entry','Another entry'},...
    'Tooltip description.');
plugin.add_param('aStringVarname',...
    'list',...
    ['aChar','Array'],...
    'Tooltip description.');
plugin.add_param('aFileChooserVarname',...
    'filechooser',...
    {'\default\Path\To\File\', 'TypeExtension'},...
    'Tooltip description.');
end


function [outputData] = arbitraryName_main(img,inputData,inputOptions,currentFrame)

% This is the main function which is called inside the main frame loop (for
% candidate search and fitting) or as the main tracking function (for
% tracking).

% img is the 2D image of the current movie (accessed by the global variable
% 'movie') in frame currentFrame. Frames start at 1.

% In case of candidate search, inputData is not given and inputOptions is
% the candidateOptions struct defined by this plugin above. outputData is
% expected to be a row array of [x,y] positions (one column each), with an
% arbitrary number of additional columns. x and y correspond to column and
% row pixel, respectively.

% In case of fitting, inputData is the row array of [x,y] positions
% (additional columns are possible but x and y are the first two) returned
% by the candidate plugin, inputOptions is the fittingOptions struct, and
% outputData is a 2D row array of particle data which MUST contain at least
% columns for x,y,z positions (put z=0 if 3D fitting does not happen), and
% possibly amplitude (or intensity), and background.

% In case of tracking, img and currentFrame are not available. inputData is
% a cell array of the outputData from the fitting stage with one cell per
% frame, inputOptions is the trackingOptions struct. outputData is expected
% to be a 2D row array with columns [Track-ID,Frame,x,y,z,[additional]].
% Both the id and the frame start at 1 as usual.
% 
% Note that all options structs, the movie, and the correction image (see
% correctMovie.m) are always accessible through the "global" key:
% 
% global globalOptions 
% global candidateOptions 
% global fittingOptions 
% global trackingOptions 
% global movie 
% global imgCorrection
% 
% That means, if your main fitting function, for example, which is called
% for every frame, would need adjacent frames for better precision, you can
% access the movie (the current frame is given by currentFrame), correct it
% (movie_corrected = correctMovie(movie(:,:,[frames_to_process])) ), and
% extract information from it. "Correcting" here means dark image
% correction and photon conversion, if enabled.

end


function [optionsStruct] = arbitraryName_init(optionsStruct)

% This function is called before the main loop starts, i.e. fittingOptions
% = arbitraryName_init(fittingOptions). You can process all options from
% the GUI here and store additional derived ones here since you can't 
% access anything inside the GUI.
% 
% For example, a filechooser in the GUI could have given you a path to a
% calibration file in inputOptions. You can load the file, process it, and
% store necessary variables here:
% 
% load(optionsStruct.someFilename); 
% calibrationVariable = someFunction(some,number,of,variables); 
% optionsStruct.anotherVarible = calibrationVariable;
% 
% Do take care not to overwrite optionsStruct with something different
% unless you want TNT to crash in a spectacular way.

end


function [outputData,optionsStruct] = arbitraryName_post(outputData,optionsStruct)

% This function is meant to post-process your data which is why both the
% data and the options struct are available for in- and output. A natural
% example would be filtering the fitted positions from the localization
% stage and saving a STORM-like image in the options struct.

end