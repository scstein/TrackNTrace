function [ im ] = BinLocalizations( positions, imsize, zoomfactor )
    
% Offset and rescale positions
positions = single(positions);
positions = positions*zoomfactor;
%positions = positions + 0.5;

% Filter localizations outside the FOV
binimsize = ceil(imsize*zoomfactor);
keep = positions(:,1)>=0;
keep = keep & positions(:,2)>=0;
keep = keep & positions(:,1)<=binimsize;
keep = keep & positions(:,2)<=binimsize;
positions = positions(keep,:);

% Bin localizations
im = zeros(binimsize,binimsize);
for i = 1:size(positions,1)
    x = round(positions(i,1));
    y = round(positions(i,2));
    if x==0
        x=x+1;
    end
    if y==0
        y=y+1;
    end
    im(y,x) = im(y,x)+1;
end

end

