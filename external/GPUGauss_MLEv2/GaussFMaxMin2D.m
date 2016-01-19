% \file GaussFMaxMin2D.m
% \author Peter Relich
% \date June 02, 2011
% \brief Max/Min estimation function invoked in the gaussMLEv2 files
% 

%% Pixel Max Min filtering
function [MaxN, MinBG] = GaussFMaxMin2D( sz, sigma, data )
    % \brief returns filtered max and min pixels of a given subregion
    %
    % Inputs:
    % \sz nxn size of the subregion
    % \sigma used in filter calculation
    % \data the subregion to search
    %
    % Outputs:
    % \MaxN maximum pixel value
    % \MingBG minimum background value
    Nfits = length(data(1,1,:));
    MaxN = zeros(1,Nfits); MinBG = 10e10*ones(1,Nfits); %starts out big
    
    norm = 0.5/sigma^2;
    
    % loop over all pixels
    for kk = 1:sz
        for ll = 1:sz
            filteredpixel = zeros(1,Nfits);
            sum = 0.0;
            for ii = 1:sz
                for jj = 1:sz
                    filteredpixel = filteredpixel + exp( -(ii-kk-2)^2*norm )...
                        *exp( -(ll-jj-2)^2*norm )*permute(data( ii,jj,: ),[1,3,2]);
                    sum = sum + exp( -(ii-kk-2)^2*norm ) * exp( -(ll-jj-2)^2*norm );
                end
            end
            filteredpixel = filteredpixel/sum;
            
            for nn = 1:Nfits
                MaxN(nn) = max(MaxN(nn), filteredpixel(nn));
                MinBG(nn) = min(MinBG(nn), filteredpixel(nn));
            end
        end
    end
end