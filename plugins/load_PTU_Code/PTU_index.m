function [head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame] = PTU_index(name,BidirectShift)
% Creates a matfile with a list of all photons with sync, tcspc and
% channel. For scans also with x, y, frame. Can be partially read from the
% drive using matfile(name) and indexing.
%
% name              file name of PTU file
% BidirectShift 
%  - false/empty    does nothing
%  - true(logical)  (Default) determines the shift with PTU_getBidirectShift
%  - double         shift in pixels
%
% Copyright (C) 2020, Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de

photons = 5e6;
if strcmp(name(end-2:end),'ptu')
    
    head = PTU_Read_Head(name);
    if ~isempty(head)
        nphot = head.TTResult_NumberOfRecords;
        freeMem = getFreeMem();
         % We need at least 2 doubles (8 bytes) per photon
        if freeMem< nphot*(8*3)
            bufferFlag = true;
        else
            bufferFlag = false;
        end
        [inputpath,inputname] = fileparts(name);
        mfile = fullfile(inputpath, makeValidMatfileName([inputname '_index.mat']));
        if bufferFlag && exist(mfile,'file')
            delete(mfile);
        end
        
        [~, ~, tmpchan, tmpmarkers] = PTU_Read(name, [1 1e4], head);
        dind = unique(tmpchan(~tmpmarkers));
        
        anzch      = 32;
        Resolution = max(1e9*head.MeasDesc_Resolution);
        chDiv      = 1e-9*Resolution/head.MeasDesc_Resolution;
        Ngate      = ceil(1e9*head.MeasDesc_GlobalResolution./Resolution); % Last bin was always empty
        
        class_t = intminclass(Ngate);         % class for the tcspc channel
        class_c = intminclass(dind);          % class for the channel
        
        if head.ImgHdr_Dimensions == 1 % Point measurement
            h1 = waitbar(0,'Reading photons');
            
            cnt      = 0; % Read photons (including markers)
            cn_ind  = 0; % Read photons (excluding markers)
            tend     = 0;
            [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
            
            while (num>0)
                
                cnt = cnt + num;
                tmpy = tmpy+tend;
                tend  = tmpy(end)+loc;
                
                ind = (tmpmarkers==0)&((tmpchan<anzch)&(tmptcspc<Ngate*chDiv)); %Remove all marker and invalid photons
                cn_num = sum(ind);
                            
                im_sync(cn_ind+1:cn_ind+cn_num,1)   = tmpy(ind);
                im_tcspc(cn_ind+1:cn_ind+cn_num,1)  = tmptcspc(ind);
                im_chan(cn_ind+1:cn_ind+cn_num,1)   = tmpchan(ind)+1;
                
                cn_ind = cn_ind+cn_num;
                waitbar(cnt/head.TTResult_NumberOfRecords,h1,'Reading photons...');
                drawnow
                
                if bufferFlag
                    matpush(mfile,'im_sync',im_sync);
                    im_sync = cast([],'like',im_sync);
                    matpush(mfile,'im_tcspc',im_tcspc);
                    im_tcspc = cast([],'like',im_tcspc);
                    matpush(mfile,'im_chan',im_chan);
                    im_chan = cast([],'like',im_chan);
                    cn_ind = 0;
                end
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
                
            end
            
            close(h1);
            
            im_col = [];
            im_line = [];
            im_frame = [];
            im_frame_index = [];
        elseif head.ImgHdr_Ident == 9 % MultiFrame Scan
            
            nx    = head.ImgHdr_PixX;
            ny    = head.ImgHdr_PixY;
            if isfield(head,'ImgHdr_MaxFrames')
                nz = head.ImgHdr_MaxFrames;
            else % Not defined at the start of the measurement
                tim_p_frame =  1./head.ImgHdr_LineFrequency*ny;
                tot_time = head.TTResult_StopAfter*1e-3;
                nz = ceil(tot_time/tim_p_frame);
            end
            class_x = intminclass(nx);  % class for the x position
            class_y = intminclass(ny);  % class for the y position
            class_f = intminclass(nz);  % class for the frame number
            

            
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
            tmpx     = [];
            chan     = [];
            markers  = [];
            dt       = zeros(ny,1);
            
            if bufferFlag % They will be cleared after each (sub)frame
                % Initialise with defined class. This casts values which are
                % appended later automatically to that class.
                im_sync  = [];
                im_tcspc = cast([],class_t);
                im_chan  = cast([],class_c);
                im_line  = cast([],class_y);
                im_col   = cast([],class_x);
                im_frame = cast([],class_f);
            else %Preallocate to avoid size changes
                im_sync  = zeros(nphot,1);
                im_tcspc = zeros(nphot,1,class_t);
                im_chan  = zeros(nphot,1,class_c);
                im_line  = zeros(nphot,1,class_y);
                im_col   = zeros(nphot,1,class_x);
                im_frame = zeros(nphot,1,class_f);
            end

            im_frame_index = zeros(nz,1);
            im_frame_index(1) = 1;
            frame_start = nan(nz,1); % Save the sync of all the frame start marker to calculate an average frame rate
            
            cn_ind   = 0;   % Count index: relative position in im_sync, im_tcspc and so on. Is reseted by writing to the matfile.
            cn_phot  = 0;   % Count photon: absolute position in im_sync, im_tcspc and so on. Is NOT reseted by writing to the matfile.
            Turns1   = [];
            Turns2   = [];
            Framechange = [];
            
            cnt      = 0; % Counts: position in PTU file (includes marker)
            tend     = 0;
            line     = 1;
            nframe   = 1;

            
            h1 = waitbar(0,['Generating index: Frame ' num2str(nframe) ' out of ' num2str(nz)]);
                        
            if head.ImgHdr_BiDirect == 0
                    
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
                
                while (num>0)
                    
                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;
                    end
                    tend  = tmpy(end)+loc;
                    
                    ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));
                    
                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = uint16([tmpx; floor(tmptcspc(ind)./chDiv)+1]);
                    chan    = uint8([chan; tmpchan(ind)+1]);         
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
                    tmpx    = tmpx(ind);
                    chan    = chan(ind);
                    markers = markers(ind);
                    
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
                        Turns2_pos = interp1(yu,ilyu,Turns2-0.1,'previous');
                    else
                        Turns1_pos = interp1(y,1:numel(y),Turns1-0.1,'next','extrap');      % Position of the next photon after the marker. For empty lines it will be in the next line or return.
                        Turns2_pos = interp1(y,1:numel(y),Turns2-0.1,'previous');  % Position of the previous photon before the marker. For empty lines it will be in the previous line.
                    end
                    
                    if numel(Turns2)>0
                        for j=1:numel(Turns2)
                            
                            t1 = Turns1(j);
                            t2 = Turns2(j);
                            
                            if any(Framechange <= t1)
                                ind = Framechange <= t1;
                                frame_start(nframe+(1:sum(ind))-1)=Framechange(ind);
                                Framechange(ind) = [];
                                nframe = nframe + sum(ind);
                                line = 1;
                                im_frame_index(nframe,1) = cn_phot+1;
                            end
                            
                            ind = Turns1_pos(j):Turns2_pos(j);      % Linear indexing is faster than logical indexing
                            if ~isempty(ind) && ~isnan(ind(1))      % ind becomes NaN if Turns1_pos and/or Turns1_pos are NaN. This can happen if the line is at the begining/end and has no photons.
                                ind_sum = numel(ind);
                                ind_out = cn_ind+1:cn_ind+ind_sum; 
                                
                                im_frame(ind_out,1)  = nframe*ones(ind_sum,1);
                                im_sync(ind_out,1)   = y(ind);
                                im_tcspc(ind_out,1)  = tmpx(ind);
                                im_chan(ind_out,1)   = chan(ind);
                                im_line(ind_out,1)   = line.*ones(ind_sum,1);
                                im_col(ind_out,1)    = 1 + floor(nx.*(y(ind)-t1)./(t2-t1+eps(t2)));
                                
                                cn_phot = cn_phot+ind_sum;
                                cn_ind = cn_ind+ind_sum;
                            end
                            
                            dt(line)  = t2-t1;
                            line = line +1;
                        end
                        
                        ind     = max(Turns2_pos):numel(y);% Max to get the last non NaN
                        y       = y(ind);
                        tmpx    = tmpx(ind);
                        chan    = chan(ind);
                        markers = markers(ind);
                        Turns1  = Turns1(Turns1>Turns2(end));
                        Turns2  = [];
                    end
                    
                    waitbar(cnt/head.TTResult_NumberOfRecords,h1,sprintf(['Generating index: Frame %d out of %d'],nframe,nz));
                    drawnow
                    
                    if bufferFlag
                        matpush(mfile,'im_frame',im_frame);
                        im_frame = cast([],'like',im_frame);
                        matpush(mfile,'im_sync',im_sync);
                        im_sync = cast([],'like',im_sync);
                        matpush(mfile,'im_tcspc',im_tcspc);
                        im_tcspc = cast([],'like',im_tcspc);
                        matpush(mfile,'im_chan',im_chan);
                        im_chan = cast([],'like',im_chan);
                        matpush(mfile,'im_line',im_line)
                        im_line = cast([],'like',im_line);
                        matpush(mfile,'im_col',im_col);
                        im_col = cast([],'like',im_col);                        
                        cn_ind = 0;
                    end
                    
                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
                    
                end
                         
            else  % bidirectional scan
                if nargin<2 || isempty(BidirectShift)
                    BidirectShift = true;
                end
                if islogical(BidirectShift) && BidirectShift
                    [BidirectShift, head] = PTU_getBidirectShift(name);
                end
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
                %                 Framechange = [];
                while (num>0)
                    
                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;
                    end
                    tend  = tmpy(end)+loc;
                    
                    ind = ((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));
                    
                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = uint16([tmpx; floor(tmptcspc(ind)./chDiv)+1]);
                    chan    = uint8([chan; tmpchan(ind)+1]);
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
                    tmpx    = tmpx(ind);
                    chan    = chan(ind);
                    markers = markers(ind);                    
                    
                    if numel(Turns2)>1 % at least two lines
                        jmax = floor(numel(Turns2)/2)*2;
                        
                        if BidirectShift
                            syncShift = BidirectShift*mean(Turns2(1:jmax)-Turns1(1:jmax))/nx;
                            Turns1(1:2:jmax) = Turns1(1:2:jmax) + floor(syncShift/2);
                            Turns2(1:2:jmax) = Turns2(1:2:jmax) + floor(syncShift/2);
                            Turns1(2:2:jmax) = Turns1(2:2:jmax) + ceil(syncShift/2);
                            Turns2(2:2:jmax) = Turns2(2:2:jmax) + ceil(syncShift/2);
                        end

                        if any(diff(y)>0)                                   % More than one photon in one sync. This can happen when measureing with more than one detector. interp1 can not deal with that so we need to sort it out.
                            [yu,ifyu] = unique(y,'first');                  % There are more efficent ways to implement this, taking advatage of the fact that y is sorted. E.g. taking diff(y)==0 and then some index shifting.
                            [~,ilyu] = unique(y,'last');
                            Turns1_pos = interp1(yu,ifyu,Turns1-0.1,'next','extrap');
                            Turns2_pos = interp1(yu,ilyu,Turns2-0.1,'previous');
                        else
                            Turns1_pos = interp1(y,1:numel(y),Turns1-0.1,'next','extrap');  % Position of the next photon after the marker. For empty lines it will be in the next line or return.
                            Turns2_pos = interp1(y,1:numel(y),Turns2-0.1,'previous');       % Position of the previous photon before the marker. For empty lines it will be in the previous line.
                        end
                        
                        
                        for j=2:2:jmax
                            
                            t1 = Turns1(j-1);
                            t2 = Turns2(j-1);
                            
                            if any(Framechange <= t1)
                                ind = Framechange <= t1;
                                frame_start(nframe+(1:sum(ind))-1)=Framechange(ind);
                                Framechange(ind) = [];
                                nframe = nframe + sum(ind);
                                line = 1;
                                im_frame_index(nframe,1) = cn_phot+1;
                            end
                            
                            ind = Turns1_pos(j-1):Turns2_pos(j-1); % Linear indexing is faster than logical indexing
                            if ~isempty(ind) && ~isnan(ind(1))       % ind becomes NaN if Turns1_pos and/or Turns1_pos are NaN. This can happen if the line is at the begining/end and has no photons.
                                ind_sum = numel(ind);
                                ind_out = cn_ind+1:cn_ind+ind_sum; 
                                
                                im_frame(ind_out,1)  = nframe*ones(ind_sum,1);
                                im_sync(ind_out,1)   = y(ind);
                                im_tcspc(ind_out,1)  = tmpx(ind);
                                im_chan(ind_out,1)   = chan(ind);
                                im_line(ind_out,1)   = line.*ones(ind_sum,1);
                                im_col(ind_out,1)    = 1 + floor(nx.*(y(ind)-t1)./(t2-t1+eps(t2)));
                                
                                cn_phot = cn_phot+ind_sum;
                                cn_ind = cn_ind+ind_sum;
                            end
                            
                            dt(line)  = t2-t1;
                            line = line +1;
                            
                            % Reverse line
                            t1 = Turns1(j);
                            t2 = Turns2(j);
                            
                            ind = Turns1_pos(j):Turns2_pos(j); %Linear indexing is faster than logical indexing
                            if ~isempty(ind) && ~isnan(ind(1))       % ind becomes NaN if Turns1_pos and/or Turns1_pos are NaN. This can happen if the line is at the begining/end and has no photons.
                                ind_sum = numel(ind);
                                ind_out = cn_ind+1:cn_ind+ind_sum; 
                                
                                im_frame(ind_out,1)  = nframe*ones(ind_sum,1);
                                im_sync(ind_out,1)   = y(ind);
                                im_tcspc(ind_out,1)  = tmpx(ind);
                                im_chan(ind_out,1)   = chan(ind);
                                im_line(ind_out,1)   = line.*ones(ind_sum,1);
                                im_col(ind_out,1)    = nx - floor(nx.*(y(ind)-t1)./(t2-t1+eps(t2)));
                                
                                cn_phot = cn_phot+ind_sum;
                                cn_ind = cn_ind+ind_sum;
                            end
                            
                            dt(line)  = t2-t1;
                            line = line +1;
                        end
                        
                        ind     = Turns2_pos(jmax):numel(y);
                        
                        if ~isempty(ind) && ~isnan(ind(1))       % ind becomes NaN if Turns1_pos and/or Turns1_pos are NaN. This can happen if the line is at the begining/end and has no photons.
                            y       = y(ind);
                            tmpx    = tmpx(ind);
                            chan    = chan(ind);
                            markers = markers(ind);
                        end
                        Turns1  = Turns1(Turns1>Turns2(jmax));
                        Turns2  = Turns2(jmax+1:end);
                    end
                    
                    waitbar(cnt/head.TTResult_NumberOfRecords,h1,sprintf(['Generating index: Frame ' num2str(nframe,'%d') ' out of ' num2str(nz)]));
                    drawnow
                    
                    if bufferFlag
                        matpush(mfile,'im_frame',im_frame);
                        im_frame = cast([],'like',im_frame);
                        matpush(mfile,'im_sync',im_sync);
                        im_sync = cast([],'like',im_sync);
                        matpush(mfile,'im_tcspc',im_tcspc);
                        im_tcspc = cast([],'like',im_tcspc);
                        matpush(mfile,'im_chan',im_chan);
                        im_chan = cast([],'like',im_chan);
                        matpush(mfile,'im_line',im_line)
                        im_line = cast([],'like',im_line);
                        matpush(mfile,'im_col',im_col);
                        im_col = cast([],'like',im_col);
                        cn_ind = 0;
                    end
                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 photons], head);
                end
            end
            
            head.ImgHdr_PixelTime = 1e9.*mean(dt)/nx/head.TTResult_SyncRate; % in ns!
            head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
            frame_start = frame_start(~isnan(frame_start));
            if numel(frame_start)>1
                head.ImgHdr_FrameFrequency = head.TTResult_SyncRate/mean(diff(frame_start));
            end
            
            im_frame_index = im_frame_index(im_frame_index>0); %#ok<*NASGU>
            
            close(h1);
            
            if ~bufferFlag
                % replacing is faster than deleting
                im_frame = im_frame(1:cn_ind,1);
                im_tcspc = im_tcspc(1:cn_ind,1);
                im_chan  = im_chan(1:cn_ind,1);
                im_line  = im_line(1:cn_ind,1);
                im_col   = im_col(1:cn_ind,1);
                im_sync  = im_sync(1:cn_ind,1);
            end
        end
        
        head.MeasDesc_Ngate = Ngate;
        head.MeasDesc_Nchan = numel(dind);
            
        % Save and clear variables which are not required anymore, save
        % variables which will not change
        if bufferFlag
            save(mfile,'head','im_frame_index','-append');
        else
            save(mfile,'head','im_sync','im_tcspc','im_line','im_col','im_chan','im_frame','im_frame_index','-v7.3');
        end
        if nargout>1
            if bufferFlag
                load(mfile,'im_sync', 'im_tcspc', 'im_chan', 'im_line', 'im_col', 'im_frame');
            end
        end
    else
        warning('Invalid PTU file. No header found.');
    end
else
    warning('Unknown file type.');
end
end

function freeMem = getFreeMem()
try % Works only on Windows
    [~,freeMem] =  memory();
    freeMem = freeMem.PhysicalMemory.Available;
catch
    freeMem =  4e9;% Use 4 GB as default.
end
end
