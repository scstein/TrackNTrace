function [intMap,zMap] = reconstructSMLM(xyz,precision,weight,magnification,orgSize,mode,modeOptions)
%RECONSTRUCTSMLM generates a super-resolved SMLM image for the intensity and the value z (e.g. a lifetime).
% The inputs xyz, precision, and orgSize have to be in the same unit e.g. pixel.
% [intMap,zMap] = reconstructSMLM(xyz,precision,weight,magnification,orgSize,mode,modeOptions)
% Inputs:
%  xyz              [nx2 OR nx3] array with x,y position. If a third column is present
%                   it is averaged across all localisations at each position.
%  precision        [nx1 OR 1x1] localisation precision; default: 0
%  weight           [nx1 OR 1x1 OR 0x0] weight for each localisation; default: 1
%  magnifcation     [1x1] number of pixel each unit is divided into; default: 10
%  orgSize          [2x1] size of the original image; default: max(xyz(:,[2 1]))
%  mode             [char] One of the following options:
%                   'Gaussian': Draws an pixel integrated gaussian with 
%                       sigma=precision at each localisation's position.
%                   'Jitter': Draws a point with a normal distributed random
%                       offset (sigma=precision) for each localisation. This is
%                       repeated n (default: 5) times. n is set by modeOptions.
%                   Otherwise: 2D histogram of localisation positions.
%
% Code for Gaussian histogram and Jitter adapted form first TNT release.
%
% Copyright (C) 2020, Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de

% Check inputs and set defaults
if nargin == 0
    return;
end
if size(xyz,2)>3 && size(xyz,1)<=3
    xyz = xyz';
end
assert(any(size(xyz,2)==[2 3]), 'Input argument xyz musst have two or three columns.');
zFlag = size(xyz,2)==3;
if ~zFlag && nargout ==2
    warning('First input needs 3 columns to generate a zMap. zMap will be empty.');
end

if exist('precision','var') && ~isempty(precision)
    precision = precision(:);
    assert(isscalar(precision) || numel(precision)==size(xyz,1), 'Input argument precision musst be scalar or same length as xyz.');
else
    precision = 0;
end

if exist('weight','var') && ~isempty(weight)
    weight = weight(:);
    assert(isscalar(weight) || numel(weight)==size(xyz,1), 'Input argument weight musst be scalar or same length as xyz.');
else
    weight = 1;
end

if ~exist('magnification','var') || isempty(magnification)
    magnification = 10;
end
assert(isscalar(magnification)&&isnumeric(magnification), 'Input argument magnification musst be a numeric scalar.');

if ~exist('orgSize','var') || isempty(orgSize)
    orgSize = ceil(max(xyz(:,[2 1]),[],1));
elseif numel(orgSize) == 1
    % assume quaratic image
    orgSize = [orgSize, orgSize];
end
assert(isvector(orgSize)&&isnumeric(orgSize), 'Input argument orgSize musst be a numeric vector.');

if all(precision==0)
    mode = 'Histogram';
elseif ~exist('mode','var') || isempty(mode)
    mode = 'Gaussian';
end
if ~exist('modeOptions','var') || isempty(modeOptions)
    modeOptions = 5;
end

% Calculate reconstruction
intMap = zeros(orgSize(1:2)*magnification);
zMap = intMap;

if numel(weight)==1
    weight = weight.*ones(size(xyz,1),1);
end
switch lower(mode)
    case 'gaussian'
        % add gaussian functions with sigma = precision
        
        % expand singletons, if necessary
        if numel(precision)==1
            precision = precision.*ones(size(xyz,1),1);
        end
        % filter localisations
        valid_loc = ~isinf(precision)&~isnan(precision)& ~isnan(weight);
        if zFlag
            valid_loc = valid_loc & ~isnan(xyz(:,3));
        end
        precision = precision(valid_loc,:);
        weight = weight(valid_loc,:);
        if zFlag
            zData = xyz(valid_loc,3);
        end
        
        pos = xyz(valid_loc,1:2)*magnification;
        pos_idx = round(pos);
        for iPos = 1:size(pos,1)
            sigma = precision(iPos)*magnification;
            halfw = ceil(3*sigma);
            if pos_idx(iPos,1)-halfw<1 || pos_idx(iPos,1)+halfw>size(intMap,2) || pos_idx(iPos,2)-halfw<1 || pos_idx(iPos,2)+halfw>size(intMap,1)
                % out of bounds
                continue;
            end
            img = weight(iPos)*gaussianMask(rem(pos(iPos,1),1),rem(pos(iPos,2),1),sigma,sigma,halfw);
            intMap(pos_idx(iPos,2)-halfw:pos_idx(iPos,2)+halfw,pos_idx(iPos,1)-halfw:pos_idx(iPos,1)+halfw) = intMap(pos_idx(iPos,2)-halfw:pos_idx(iPos,2)+halfw,pos_idx(iPos,1)-halfw:pos_idx(iPos,1)+halfw)...
                +img;
            if zFlag
                zMap(pos_idx(iPos,2)-halfw:pos_idx(iPos,2)+halfw,pos_idx(iPos,1)-halfw:pos_idx(iPos,1)+halfw) = zMap(pos_idx(iPos,2)-halfw:pos_idx(iPos,2)+halfw,pos_idx(iPos,1)-halfw:pos_idx(iPos,1)+halfw)...
                    +zData(iPos)*img;
            end
        end
        if zFlag
            zMap(intMap>0) = zMap(intMap>0)./intMap(intMap>0);
        end
        
    case 'jitter'
        % create jittered histogram where an average histogram of noisy positions is calculated
        % expand singletons, if necessary
        if numel(precision)==1
            precision = precision.*ones(size(xyz,1),1);
        end
        % filter localisations
        valid_loc = ~isinf(precision)&~isnan(precision)& ~isnan(weight);
        if zFlag
            valid_loc = valid_loc & ~isnan(xyz(:,3));
        end
        precision = precision(valid_loc,:);
        weight = weight(valid_loc,:);
        if zFlag
            zData = xyz(valid_loc,3);
        end
        
        pos = xyz(valid_loc,[2,1]);
        
        for i=1:max(1,modeOptions)
            pos_jitter = ceil(magnification*(pos+randn(size(pos,1),2).*[precision,precision]));
            valid_loc = 0 < pos_jitter(:,1) & pos_jitter(:,1) <= magnification*orgSize(2) ...
                      & 0 < pos_jitter(:,2) & pos_jitter(:,2) <= magnification*orgSize(1) ;

            intMap = intMap + accumarray(pos_jitter(valid_loc,:),weight(valid_loc),orgSize(1:2)*magnification);
            if zFlag
                zMap = zMap + accumarray(pos_jitter(valid_loc,:),weight(valid_loc).*zData(valid_loc),orgSize(1:2)*magnification);
            end
        end
        if zFlag
            zMap(intMap>0) = zMap(intMap>0)./intMap(intMap>0);            
        end
    otherwise % Normal histogram
        
        valid_loc = 0 < xyz(:,1) & xyz(:,1) <= orgSize(2) ...
                  & 0 < xyz(:,2) & xyz(:,2) <= orgSize(1) ...
                  & ~isnan(weight);
        if zFlag
            valid_loc = valid_loc & ~isnan(xyz(:,3));
        end
        
        pixpos = ceil(magnification*xyz(valid_loc,[2 1]));
        intMap = accumarray(pixpos,weight(valid_loc),orgSize(1:2)*magnification);
        
        if zFlag
            zMap = accumarray(pixpos,weight(valid_loc).*xyz(valid_loc,3),orgSize(1:2)*magnification);
            zMap(intMap>0) = zMap(intMap>0)./intMap(intMap>0);            
        end
end %SWITCH
end

function [img] = gaussianMask(x,y,sigma_x,sigma_y,halfw)
    % normalised, integrated 2D gaussian
    [x_grid,y_grid] = meshgrid (-halfw:halfw,-halfw:halfw);

    img = 0.5.*(erfc(-(x_grid-x+0.5)/(sqrt(2)*sigma_x))-erfc(-(x_grid-x-0.5)/(sqrt(2)*sigma_x))).*...
          0.5.*(erfc(-(y_grid-y+0.5)/(sqrt(2)*sigma_y))-erfc(-(y_grid-y-0.5)/(sqrt(2)*sigma_y)));
    img = img./sum(img(:));
end