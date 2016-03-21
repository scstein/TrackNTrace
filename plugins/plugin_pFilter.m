function [plugin] = plugin_pFilter()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'p-value filtering';

% Type of plugin.
% 1: Candidate detection
% 2: Spot fitting
% 3: Tracking
type = 1;

% The functions this plugin implements
mainFunc =  @findCandidates_pFiltering;

% Description of output parameters
outParamDescription = {'x';'y'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = 'Candidate detection based on a p-value hypothesis test against local background. \n\nIt is assumed that image pixels belonging to actual emitters cannot be purely described by flat background noise. Therefore, all pixels are compared against a Gaussian background. If a cluster of pixels fails this test, it is assumed to belong to an emitter and will count as a candidate.';

% Deactivate this plugin's parallel processing, as plugin relys on persistent
% variables depending on previous loop iterations
plugin.useParallelProcessing = false;

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% types are int, float, bool, list, string, filechooser
plugin.add_param('pThreshold',...
    'float',...
    {0.92, 0, 1},...
    'p-Value threshold above which pixels can be identified as candidates. \nHigher means less but higher quality candidates.');
plugin.add_param('iterationCount',...
    'int',...
    {10, 0, inf},...
    'Background is estimated every N times in the main loop. \nSet to higher number to speed up function.');
plugin.add_param('clusterSize',...
    'int',...
    {4,2,inf},...
    'Number of pixels belonging to a local maximum candidate (locally connected cluster). \nShould be similar to the particle size.');
end


%   -------------- User functions --------------

function [candidateData] = findCandidates_pFiltering(img,options,currentFrame)
% Find particle candidates in fluorescence images based on intensity, local
% background and size.
% 
% INPUT:
%     img: 2D matrix of pixel intensities, data type and normalization
%     arbitrary.
%
%     candidateOptions: Struct of input parameters provided by GUI.
% 
%     currentFrame: Integer, current movie frame in main loop.
%
% OUTPUT:
%     candidatePos - 2D array of Nx2 matrix of particle candidate positions
%     [column pixel, row pixel] without subpixel position. Middle of upper
%     left pixel would be [1,1].

persistent img_bck_mean img_bck_std

if mod(currentFrame,options.iterationCount)==0 || isempty(img_bck_mean)
    [img_bck_mean,img_bck_std] = calcBG(img,size(img,1),size(img,2),options.clusterSize);
end

% obtain threshold binary image
pMap = normcdf(img,img_bck_mean,img_bck_std)>options.pThreshold;

% find clustes
connected_clusters = bwconncomp(pMap,4);
connected_clusters.PixelIdxList = connected_clusters.PixelIdxList(find(cellfun(@(var) numel(var)>options.clusterSize-1,connected_clusters.PixelIdxList)));
connected_clusters.NumObjects = numel(connected_clusters.PixelIdxList);

% get their centroids
centroids = regionprops(connected_clusters,'Centroid');
candidateData = cat(1,centroids.Centroid);

end %END pFiltering
