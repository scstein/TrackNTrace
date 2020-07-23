function [bidirect_shift, head] = PTU_getBidirectShift(name)
% [bidirect_shift, head] = PTU_getBidirectShift(name)
% Reads the first frames of a bidirectional scanned measurement and determines
% the shift between forward and reverse scan with subpixel accuraccy.
% This is done with an adaptive grid search which minimise the difference 
% between adjacent lines.
%
% Copyright (C) 2020, Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de


photons = 5e6;
if strcmp(name(end-2:end),'ptu')
    
    head = PTU_Read_Head(name);
    if ~isempty(head)
        
        if ~isfield(head,'ImgHdr_BiDirect') || head.ImgHdr_BiDirect == 0
            bidirect_shift = 0;
            return;
        end
        
        nphot = head.TTResult_NumberOfRecords;        
        nx    = head.ImgHdr_PixX; %
        ny    = head.ImgHdr_PixY;
        if isfield(head,'ImgHdr_MaxFrames')
            nz = head.ImgHdr_MaxFrames;
        else
            nz = 1;
        end

%         class_t = 'uint16';         % class for the tcspc channel
%         class_c = 'uint8';          % class for the channel
        class_x = 'single';         % class for the x position, non-interger position for later binning
        class_y = intminclass(ny);  % class for the y position
%         class_f = intminclass(nz);  % class for the frame number
        
        if head.ImgHdr_Ident ~= 9
            error('Not supported scan hardware');
        elseif head.ImgHdr_Ident == 9 % FLIMBee
            
          if isfield(head,'ImgHdr_MaxFrames')
              nz = head.ImgHdr_MaxFrames;
          else
              tim_p_frame =  1./head.ImgHdr_LineFrequency*ny;
              tot_time = head.TTResult_StopAfter*1e-3;
              nz = ceil(tot_time/tim_p_frame);
          end
          class_f = intminclass(nz);
          
          
          anzch      = 32;
          Resolution = max(1e9*head.MeasDesc_Resolution);
          chDiv      = 1e-9*Resolution/head.MeasDesc_Resolution;
          Ngate      = ceil(1e9*head.MeasDesc_GlobalResolution./Resolution); % Last bin was always empty
          
          LineStart = 4;
          LineStop  = 2;
          Frame     = 3;
          
          if isfield(head,'ImgHdr_LineStart')
              LineStart = 2^(head.ImgHdr_LineStart-1);
          end
          if isfield(head,'ImgHdr_LineStop')
              LineStop = 2^(head.ImgHdr_LineStop-1);
          end
          if isfield(head,'ImgHdr_Frame')
              Frame = 2^(head.ImgHdr_Frame-1);
          end
          y        = [];
%           tmpx     = [];
%           chan     = [];
          markers  = [];
          dt       = zeros(ny,1);
          
          
%           im_sync  = zeros(photons,1);
%           im_tcspc = zeros(photons,1,class_t);
%           im_chan  = zeros(photons,1,class_c);
          im_line  = zeros(photons,1,class_y);
          im_col   = zeros(photons,1,class_x);
          im_frame = zeros(photons,1,class_f);
          
          im_frame_index = zeros(1,nz);
          im_frame_index(1) = 1;
          
          cn_phot = 0;
          Turns1   = [];
          Turns2   = [];
          Framechange = [];
          
          cnt      = 0;
          tend     = 0;
          line     = 1;
          nframe    = 1;
          
          
          h1 = waitbar(0,['Frame ' num2str(nframe) ' out of ' num2str(nz)]);
          
          % bidirectional scan
          
          [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
          %                 Framechange = [];
          while (num>0 && nframe < 2) % process PHOTONS entries but at least the first full frame 
              
              cnt = cnt + num;
              if ~isempty(y)
                  tmpy = tmpy+tend;
              end
              tend  = tmpy(end)+loc;
              
              ind = ((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));
              
              y       = [y; tmpy(ind)];                         %#ok<AGROW>
%               tmpx    = uint16([tmpx; floor(tmptcspc(ind)./chDiv)+1]);
%               chan    = uint8([chan; tmpchan(ind)+1]);
              markers = uint8([markers; tmpmarkers(ind)]);
              
              if LineStart==LineStop
                  tmpturns = y(markers==LineStart);
                  if numel(Turns1)>numel(Turns2)           % first turn is a LineStop
                      Turns1 = [Turns1; tmpturns(2:2:end)];      %#ok<AGROW>
                      Turns2 = [Turns2; tmpturns(1:2:end)];      %#ok<AGROW>
                  else
                      Turns1 = [Turns1; tmpturns(1:2:end)];      %#ok<AGROW>
                      Turns2 = [Turns2; tmpturns(2:2:end)];      %#ok<AGROW>
                  end
              else
                  Turns1 = [Turns1; y(markers==LineStart)]; %#ok<AGROW>
                  Turns2 = [Turns2; y(markers==LineStop)];  %#ok<AGROW>
              end
              Framechange = [Framechange; y(markers==uint8(Frame))];
              ind     = ~markers; %Remove markers
              y       = y(ind);
%               tmpx    = tmpx(ind);
%               chan    = chan(ind);
              markers = markers(ind);
              
              
              if numel(Turns2)>1 % at least two lines
                  % Get scanning times (all in syncs)
                  Turn_centertime = median(Turns2-Turns1(1:(numel(Turns2))));  % Time between line start and stop marker
                  Turn_linetime   = median(Turns1(2:end)-Turns1(1:end-1));     % Time between start of one and the next line
                  Turn_switchtime = Turn_linetime-Turn_centertime;             % Time spend inbetween two lines
                  
                  % Find the photons between the markers. For empty lines
                  % Turns2_pos will be smaller than Turns1_pos which will
                  % generate an empty index.
                  % Turns1/2 are shifted by 0.1 syncs in order to include
                  % photons at the sync of Turn1 but to exclude photon at
                  % the sync of Turn2. This avoids potential out of range
                  % errors in im_col.
                  if any(diff(y)>0)                                   % More than one photon in one sync. This can happen when measureing with more than one detector. interp1 can not deal with that so we need to sort it out.
                      [yu,ifyu] = unique(y,'first');                  % There are more efficent ways to implement this, taking advatage of the fact that y is sorted. E.g. taking diff(y)==0 and then some index shifting.
                      [~,ilyu] = unique(y,'last');
                      Turns1_pos = interp1(yu,ifyu,Turns1-0.1,'next','extrap');
                      Turns2_pos = interp1(yu,ilyu,Turns2-0.1,'previous','extrap');
                      
                      Turns1_pos_extended = interp1(yu,ifyu,Turns1-Turn_switchtime/2-0.1,'next','extrap');
                      Turns2_pos_extended = interp1(yu,ilyu,Turns2+Turn_switchtime/2-0.1,'previous','extrap');
                  else
                      Turns1_pos = interp1(y,1:numel(y),Turns1-0.1,'next','extrap');      % Position of the next photon after the marker. For empty lines it will be in the next line or return.
                      Turns2_pos = interp1(y,1:numel(y),Turns2-0.1,'previous','extrap');  % Position of the previous photon before the marker. For empty lines it will be in the previous line.
                      
                      Turns1_pos_extended = interp1(y,1:numel(y),Turns1-Turn_switchtime/2-0.1,'next','extrap');
                      Turns2_pos_extended = interp1(y,1:numel(y),Turns2+Turn_switchtime/2-0.1,'previous','extrap');
                  end
              
                  for j=2:2:numel(Turns2)
                      
                      t1 = Turns1(j-1);
                      t2 = Turns2(j-1);
                      
                      if any(Framechange <= t1)
                          ind = Framechange <= t1;
                          Framechange(ind) = [];
                          nframe = nframe + sum(ind);
                          line = 1;
                          im_frame_index(nframe) = cn_phot+1;
                      end
                      
                      ind = Turns1_pos_extended(j-1):Turns2_pos_extended(j-1); %Linear indexing is faster than logical indexing
                      
                      if ~isempty(ind) && ~isnan(ind(1))      % ind becomes NaN if Turns1_pos and/or Turns1_pos are NaN. This can happen if the line is at the begining/end and has no photons.
                          ind_sum = numel(ind);
                          
                          im_frame(cn_phot+1:cn_phot+ind_sum)  = nframe*ones(ind_sum,1);
                          %                       im_sync(cn_phot+1:cn_phot+ind_sum)   = y(ind);
                          %                       im_tcspc(cn_phot+1:cn_phot+ind_sum)  = tmpx(ind);
                          %                       im_chan(cn_phot+1:cn_phot+ind_sum)   = chan(ind);
                          im_line(cn_phot+1:cn_phot+ind_sum)   = line.*ones(ind_sum,1);
                          im_col(cn_phot+1:cn_phot+ind_sum)    = 1 + (nx.*(y(ind)-t1)./(t2-t1+eps(t2))); % Runs from 1 to nx+1
                          
                          cn_phot = cn_phot+ind_sum;
                      end
                      
                      dt(line)  = t2-t1;
                      line = line +1;
                      
                      % Reverse line
                      t1 = Turns1(j);
                      t2 = Turns2(j);
                      
                      ind = Turns1_pos_extended(j):Turns2_pos_extended(j); %Linear indexing is faster than logical indexing
                      
                      if ~isempty(ind) && ~isnan(ind(1))      % ind becomes NaN if Turns1_pos and/or Turns1_pos are NaN. This can happen if the line is at the begining/end and has no photons.
                          ind_sum = numel(ind);
                          
                          im_frame(cn_phot+1:cn_phot+ind_sum)  = nframe*ones(ind_sum,1);
                          %                       im_sync(cn_phot+1:cn_phot+ind_sum)   = y(ind);
                          %                       im_tcspc(cn_phot+1:cn_phot+ind_sum)  = tmpx(ind);
                          %                       im_chan(cn_phot+1:cn_phot+ind_sum)   = chan(ind);
                          im_line(cn_phot+1:cn_phot+ind_sum)   = line.*ones(ind_sum,1);
                          im_col(cn_phot+1:cn_phot+ind_sum)    = nx + 1 - (nx.*(y(ind)-t1)./(t2-t1+eps(t2))); % Runs from nx+1 to 1
                          
                          cn_phot = cn_phot+ind_sum;
                      end
                      
                      dt(line)  = t2-t1;
                      line = line +1;
                  end
                  
                  ind     = Turns2_pos(floor(numel(Turns2_pos)/2)*2):numel(y);
                  y       = y(ind);
%                   tmpx    = tmpx(ind);
%                   chan    = chan(ind);
                  markers = markers(ind);
                  Turns1  = Turns1(Turns1>Turns2(end));
                  Turns2  = [];
              end
              
              waitbar(cnt/head.TTResult_NumberOfRecords,h1,sprintf(['Frame ' num2str(nframe,'%d') ' out of ' num2str(nz)]));
              %                                     waitbar(line/ny,h2,sprintf(['Frame ' num2str(frame)]));
              drawnow
              
              [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
          end
          
          % replacing is faster than deleting
          im_frame = im_frame(1:cn_phot);
%           im_tcspc = im_tcspc(1:cn_phot);
%           im_chan  = im_chan(1:cn_phot);
          im_line  = im_line(1:cn_phot);
          im_col   = im_col(1:cn_phot);
%           im_sync  = im_sync(1:cn_phot);
          
          waitbar(1,h1,'Aligning lines...');
             
          %%
          % The two accumarrays need a lot of free memory for the
          % construction of the subs matrix.
          padding = ceil((Turn_switchtime*nx/Turn_centertime)/2*1.2); % Take 20% more to avoid out of bound errors in accumarry
          subpixelsize = 5; % nm
          nsubpixel = ceil(head.ImgHdr_PixResol/subpixelsize*1000);
          % Make a corse image of all read frames. Determin the frame
          % number for which we have accumulated 100 photon per corse pixel
          % or the gradient decreased to half.
          % This helps to supress background for STORM like measurements.
          coarsefactor = ceil(0.5/head.ImgHdr_PixResol);% roughly 500 nm
          if im_frame(end)>1
              img = accumarray([ceil(double(im_line)/coarsefactor) ceil((padding+ceil(im_col))/coarsefactor) im_frame],1,[ceil([ny (nx+2*padding)]./coarsefactor) im_frame(end)]);
              maxframe = find(cumsum(squeeze(sum(img,1:2)/nx/ny))>100,1);
              
              [gradx,grady] = gradient(img);
              grad = movmean(squeeze(mean(mean(sqrt(gradx.^2+grady.^2),1),2)),5);
              
              maxframe = max([maxframe find(grad<mean(grad([1 end])),1)])+1;
          else
              maxframe = [];
          end
          if isempty(maxframe) || maxframe<2 % if there is only one frame or find is returning [] or 1;
              maxframe = 2;
          end
          % Bin image with high resolution accumulated over the selected
          % frames. Apply movemean with odd window side to avoid shifting.
          % This is not the fastest way to do it but the recasting and
          % clearing helps to save some memory.
          ind = im_frame<maxframe;
          clear im_frame;
          im_line = cast(im_line(ind),class(im_col));
          im_line(:,2) = nsubpixel*padding+ceil(nsubpixel*im_col(ind)); % im_line now becomes subs
          clear ind im_col;
          img = movmean(accumarray(im_line,1,[ny nsubpixel*(nx+2*padding)]),floor(nsubpixel/2)*2+1,2);
          clear im_line;
          isinrange = @(x, vrange) x >= vrange(1) & x <= vrange(2);
          mask = isinrange(1:size(img,2),[nsubpixel*padding size(img,2)-nsubpixel*padding]);
          imgshiftfun = @(img,shift)cat(3,circshift(img(1:2:end,:),-floor(shift/2),2),circshift(img(2:2:end,:),ceil(shift/2),2));
          imgdifffun = @(img_shifted,mask)sum((img_shifted(:,mask,1)-img_shifted(:,mask,2)).^2,'all')+sum((img_shifted(2:end,mask,1)-img_shifted(1:end-1,mask,2)).^2,'all');
          % to reassamble shifted image: @(img)reshape(permute(img,[3 1 2]),[],size(img,2))
%           [~,shiftmin] = min(movmean(arrayfun(@(sh)(sum((img(1:2:end,:)-circshift(img(2:2:end,:),sh,2)).^2,'all')),shiftrange),1));
          % Search 10th padding resolution, at least 5 px resolution,
          % maximum range 25px
          shiftrange = ceil(linspace(-min((padding),25)*nsubpixel,min(padding,25)*nsubpixel,10));
          [~,shiftmin] = min(arrayfun(@(sh)imgdifffun(imgshiftfun(img,sh),mask),shiftrange));
          bidirect_shift = 1/nsubpixel*shiftrange(shiftmin);          
          
          % Search with pixel precission
          shiftrange = ceil(bidirect_shift*nsubpixel-(padding/10*nsubpixel)):nsubpixel:floor(bidirect_shift*nsubpixel+(padding/10*nsubpixel));
          [~,shiftmin] = min(arrayfun(@(sh)imgdifffun(imgshiftfun(img,sh),mask),shiftrange));
          bidirect_shift = 1/nsubpixel*shiftrange(shiftmin);
          
          % Refine subpixel position
          shiftrange = (-nsubpixel:nsubpixel)+round(bidirect_shift*nsubpixel);
          [~,shiftmin] = min(arrayfun(@(sh)imgdifffun(imgshiftfun(img,sh),mask),shiftrange));
          bidirect_shift = 1/nsubpixel*shiftrange(shiftmin);

          %%
          head.ImgHdr_BiDirectFittedOffset = bidirect_shift;
          head.ImgHdr_PixelTime = 1e9.*mean(dt)/nx/head.TTResult_SyncRate;
          head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
          close(h1);

        end

    else
        warning('Invalid PTU file. No header found.');
    end
end
end
