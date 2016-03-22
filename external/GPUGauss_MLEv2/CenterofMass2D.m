% \file CenterofMass2D.m
% \author Peter Relich
% \date May 24, 2011
% \brief 2D Center of Mass function invoked in the gaussMLEv2 files
%

%% 2D Center of Mass
function [x,y] = CenterofMass2D( sz, data )
    % \brief Compute the 2D center of mass of a subregion
    %
    % Inputs:
    % \sz nxn size of the subregion
    % \data subregion to search
    % 
    % Outputs:
    % \x x coordinate to return
    % \y y coordinate to return
    
    Nfits = length(data(1,1,:));
    
    tmpx = zeros(1,Nfits);
    tmpy = zeros(1,Nfits);
    tmpsum = zeros(1,Nfits);
    
    for ii = 1:sz
        for jj = 1:sz
            tmpx = tmpx + permute(data( ii,jj,: ),[1,3,2])*(ii-1);
            tmpy = tmpy + permute(data( ii,jj,: ),[1,3,2])*(jj-1);
            tmpsum = tmpsum + permute(data( ii,jj,: ),[1,3,2]);
        end
    end
    
    x = tmpx./tmpsum;
    y = tmpy./tmpsum;
end

