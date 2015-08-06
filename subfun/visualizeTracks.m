function visualizeTracks( movie, trajectoryData, FPS, traj_lifetime, n_colors, use_bw, manContrast )
% USAGE: visualizeTracks(movie, trajectoryData)
% [ Full USAGE: visualizeTracks(movie, trajectoryData, FPS, traj_lifetime, n_colors, use_bw) ]
% 
% Visualizer for tracks computed by the tracker.
%
% Input:
%   movie: 3D matrix (rows,cols,frames) of the analyzed movie
%   trajectoryData: 2D matrix with columns particleID|frame|x|y|...
%                   This is the output of the tracker
%   FPS: frames per second to play movie with | default: 30
%   traj_lifetime: trajectories are kept for #traj_lifetime frames after
%                  the particle has vanished. | default: 0
%   n_colors: number of colors used to display trajectories | default: 20
%             Colors are generated using distinguishable_colors.m by 
%             Timothy E. Holy (Matlab File Exchange).
%   use_bw: black/white image, otherwise colormap hot | default false
%   manContrast: Manual contrast adjustment for first frame
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2015
%

% input parsing
if nargin <3 || isempty(FPS)
    FPS = 30;
end

if nargin <4 || isempty(traj_lifetime)
   traj_lifetime = 0;
end

if nargin < 5 || isempty(n_colors)
    n_colors = 20;
end

if nargin <6 || isempty(use_bw)
   use_bw = false;
end

if nargin <7 || isempty(manContrast)
   manContrast = false;
end

% Convert data into cell array which is better for plotting
id_tracks = unique(trajectoryData(:,1));
n_tracks = numel(id_tracks);

cell_traj = cell(n_tracks ,1);
cnt = 1;
for iTrack = 1:n_tracks
    cell_traj{iTrack} = trajectoryData( trajectoryData(:,1)== id_tracks(cnt) , 2:4);
    cnt = cnt+1;
end

% create colors
if use_bw % background color
    bg = 'k';
else
    bg = 'r';
end

track_colors = repmat( distinguishable_colors(n_colors, bg), ceil(n_tracks/n_colors) ,1);
track_colors = track_colors(1:n_tracks,:);



% Play the movie
timePerFrame = 1/FPS;
elapsed_time = 0;

h = figure;
for iF = 1:size(movie,3)
    if(iF ==1)
        xl = [0.5,size(movie,2)+0.5];
        yl = [0.5,size(movie,1)+0.5];
        frame = movie(:,:,1);
        zl = [min(frame(:)), max(frame(:))];
    else
       xl = xlim;   
       yl = ylim;
    end
    
    % Skip frame if computer is too slow
    if elapsed_time >timePerFrame
       elapsed_time = elapsed_time - timePerFrame;
       continue;
    end
    tic_start = tic;
   
    % Plot movie frame
    imagesc(movie(:,:,iF)); axis image; colormap gray; 
    if use_bw
       colormap gray; 
    else
       colormap hot; 
    end     
    title(sprintf('Frame %i/%i',iF,size(movie,3)));
   
   % Draw the tracks of currently visible particles
   hold on;
   for iTrack = 1:n_tracks
       if iF < cell_traj{iTrack}(1,1) || iF > cell_traj{iTrack}(end,1) + traj_lifetime % skip all particles not visible (any longer)
           continue
       end
       mask_toPlot = cell_traj{iTrack}(:,1)<=iF;
       plot(cell_traj{iTrack}(mask_toPlot, 2), cell_traj{iTrack}(mask_toPlot, 3), '.--','Color',track_colors(iTrack,:));
   end
   hold off;
   
   % allow manual tweaking of contrast on first frame
   if(iF == 1 && manContrast)
       imcontrast;
       pause;
       zl = caxis;
       tic_start = tic;
   end
   
   % save the axis limits in case the user zoomed
   xlim(xl);
   ylim(yl);
   caxis(zl);
   
   drawnow;
   pause(1/FPS); 
    
   elapsed_time = elapsed_time + toc(tic_start)- timePerFrame;
end

end

