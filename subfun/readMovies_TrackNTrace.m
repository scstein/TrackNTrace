function [movie_list,dark_stack] = readMovies_TrackNTrace(filename)
% function [movieList,darkStack] = readMovies_TrackNTrace(filename)
% Find filenames of all relevant TIF images/movies from directories
% provided by user to be used for particle tracking
%
% INPUT:
%     filename: struct of two strings, one for movie folder and one for
%     dark image stack folder. See RunTrackNTrace for details.
%
% OUTPUT
%     movieList: cell array of full paths of all movies
%
%     darkStack: 3D array of correction images. Multiple images get stacked
%     along the third dimension.


% Parse input
dir_mov = filename.movies;
dir_dark = filename.dark_movies;

dark_stack = [];
if ~isempty(dir_dark)
    files_dark = dir(dir_dark); %get files
    files_dark = files_dark(find(cellfun(@(var) ~isempty(strfind(var,'dark') & strfind(var,'tif')),{files_dark(:).name}))); %#ok<FNDSB>
    dark_stack = [];
    [dir_dark,~,~] = fileparts(dir_dark);
    dir_dark = [dir_dark,filesep]; %little trick to return filename if dir_dark is not a directory but a single path string to a stack
    
    for i=1:numel(files_dark)
        dark_movie = read_tiff([dir_dark,files_dark(1).name],false);
        dark_img = CalculateDark(dark_movie);
        dark_stack = cat(3,dark_stack,dark_img);
    end
end


files_movie = dir(dir_mov);
files_movie = files_movie(find(cellfun(@(var) (isempty(~strfind(var,'dark')) & ~isempty(strfind(var,'tif'))),{files_movie(:).name}))); %#ok<FNDSB>
[dir_mov,~,~] = fileparts(dir_mov);
dir_mov = [dir_mov,filesep]; %see above

movie_list = cell(numel(files_movie),1);
for i=1:numel(files_movie)
    movie_list(i) = {[dir_mov,files_movie(i).name]};
end