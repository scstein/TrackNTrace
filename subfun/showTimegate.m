function [fighdl] = showTimegate(inputfile,importOptions,timegate,nsFlag)
%SHOWTIMEGATE Plots a TCPSC histogram indicating the choosen timegate
%   [fighdl] = showTimegate(inputfile,importOptions,timegate,nsFlag)
%
%   Copyright (C) 2020  Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de
    maxPhotons = 1e6;
    
    if nargin < 3 || isempty(timegate)
        timegate = 0;
    end
    if nargin < 4 || isempty(nsFlag)
        nsFlag = false;
    end
    
    if nargin(importOptions.info.getTCSPC) == 3
        [tcspcdata,Resolution] = importOptions.info.getTCSPC(inputfile,maxPhotons,importOptions); % Resolution in s
    else
        [tcspcdata,Resolution] = importOptions.info.getTCSPC(inputfile,maxPhotons); % Resolution in s
    end
        
    if nsFlag
        timegate = timegate./(Resolution*1e9);
    end
    timegate(end+1:2) = inf; % Ensure that second entry exists
        
    tcspc = accumarray(tcspcdata,1);
    tcspc_tg = accumarray(tcspcdata,tcspcdata>timegate(1)&tcspcdata<timegate(2),size(tcspc));
    
    semilogy([tcspc tcspc_tg]);
    
    xlabel('Time bin');
    ylabel('Counts');
    if ~isempty(tcspc)
        legend({'TCSPC',sprintf('Within the timegate: %.1f%%',100*sum(tcspc_tg)/sum(tcspc))});
    end
end

