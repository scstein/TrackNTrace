function [plugin] = plugin_fitLT()

%    -------------- Definition of plugin --------------
% Name of the component these options are for
name = 'fit lifetime';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
type = 4;

% The functions this plugin implements
mainFunc =  @fitLT;

% Description of output parameters
outParamDescription = {'Track-ID';'Frame';'x';'y';'z';'Amplitude';'Background';'sigma';'tau-fast';'nphoton'};

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Add initFunction
plugin.initFunc = @init_fitLT;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = ['Fit the lifetimes of localisations. \n\n',...
    'Extractes the TCSPC for each localisation or track \n'...
    'and determines the lifetime.'];

% Add parameters
% read comments of function TNTplugin/add_param for HOWTO
% plugin.add_text('TCSPC options');
plugin.add_param('position',...
    'list',...
    {'Preserve','Average','Refit'},...
    ['Chose whether the positions from the tracking are kept, averaged over the whole track or \n'...
     'refitted with an isotropic, integrated MLE Gaussian on a sum image of all frames of the track.']);
plugin.add_param('maskRadius',...
    'float',...
    {2, 0, inf},...
    'Mask radius in units of the PSF.');
plugin.add_param('sumTCSPC',...
    'bool',...
    true,...
    'If checked the photons for the TCSPC are accumulated over all frames of the track and fitted once.');
plugin.add_param('exportTCSPC',...
    'bool',...
    false,...
    'If checked the TCSPCs are exported as part of the postprocOptions.');
plugin.add_text('Lifetime fitting options');
plugin.add_param('fitType',...
    'list',...
    {'fast lifetime', 'monoexponential MLE', 'monoexponential PM', 'average DistFit', 'maximum correlation'},...
    ['Type of fit:\n'...
     '  fast lifetime: standard deviation of the arivial times. No fit. Ignores the following parameters.\n'...
     '  monoexponential MLE: fitting a monoexponential decay using a Nealder-Mead algorithm with a maximum likelihood estimator.\n'...
     '  monoexponential PM: pattern matching (MLE) with a library of monoexponential decays. Limited to 1000 lifetimes between minLT and maxLT but faster than fitting.\n'...
     '  average DistFit: fits the decay with a distribution of decays between minLT and maxLT and takes the amplitude-weighted rate average.\n'...
     '  maximum correlation: correlates the decay with a distribution of decays between minLT and maxLT and takes the highest correlation.\n']);
plugin.add_param('cutoff',...
    'float',...
    {0.3, 0, inf},...
    'Cut off after the maxium of the peak in ns.');
plugin.add_param('minAmp',...
    'int',...
    {100, 1, inf},...
    'Minimum number of photons to do a fit.');
plugin.add_param('minLT',...
    'float',...
    {0.5, 0, inf},...
    'Lower bound for fitting the lifetime in ns.');
plugin.add_param('maxLT',...
    'float',...
    {5, 0, inf},...
    'Upper bound for fitting the lifetime in ns.');
plugin.add_param('attempts',...
    'int',...
    {1, 1, inf},...
    'Number of fit attempts. The best fit is returned. (Only for MLE)');
plugin.add_param('tolerance',...
    'float',...
    {1e-8, 1e-30, 1e-1},...
    'Fit tolerance. (Only for MLE)');
end


%   -------------- User functions --------------
function options = init_fitLT(options,globalOptions,importOptions)
% We need access to the globalOptions to get the filename, frame range and
% binning. This is implimented as init function to avoid blocking of
% parallel execution.
    global filename_movie;
    options.filename_movie = filename_movie;
    options.hasTCSPC = importOptions.info.hasTCSPC;
    if globalOptions.previewMode
        options.framebinning = [globalOptions.binFrame, globalOptions.firstFrameTesting, globalOptions.lastFrameTesting];
    else
        options.framebinning = [globalOptions.binFrame, globalOptions.firstFrame, globalOptions.lastFrame];
    end
    if options.hasTCSPC
        if isfield(importOptions,'initFunc')
            options.cacheFunc = importOptions.initFunc(importOptions);
        else
            options.cacheFunc = [];
        end
        options.accumFunc = importOptions.postFunc(importOptions);
        if globalOptions.useTimegate
            options.timegate = [globalOptions.tgStart globalOptions.tgEnd];
        else
            options.timegate = [];
        end
    end
    
    % Set name of output variables
    global trackingOptions % Needed as we drag along names specified here
    options.outParamDescription = [trackingOptions.outParamDescription(:);'tau-fast';'nphoton'];

end

function [postprocData,options] = fitLT(trackingData,options)
	% Extracts the TCSPC from an area with with size maskRadius around the position
	% for each track/localisation. The lifetime can determined from the TCSPCs using
	% a tailfit with a distribution of decays or with a monoexponential MLE fit.
	% Copyright (C) 2020, Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de

    postprocData = [];
        
    msgAccumulator = ''; % Needed for rewindable command line printing (rewPrintf subfunction)
    mainTime_start = tic;
    mainTime = []; %#ok<NASGU>
    lastElapsedTime = 0;
    
    if isempty(trackingData)
        postprocData = trackingData;    
    elseif options.hasTCSPC
        if isempty(options.cacheFunc)
            [head,img]=options.accumFunc(options.filename_movie,{'head','tag'},[1 1 1]);
            cacheORfile = options.filename_movie;
        else
            [head,img,cacheORfile]=options.cacheFunc(options.filename_movie,{'head','tag','index'},[1 1 1]);
        end
        imgSZ = size(img);
        imgSZ = imgSZ(1:2);
        if ~issorted(trackingData(:,1)) %Tracking data has to be sorted by id (TNTnearestNeighbor does this by default).
            trackingData = sortrows(trackingData,1);
        end
        if options.sumTCSPC
            tcspc_mol = zeros(1,trackingData(end,1),head.MeasDesc_Ngate,head.MeasDesc_Nchan);
        else
            tcspc_mol = zeros(1,size(trackingData,1),head.MeasDesc_Ngate,head.MeasDesc_Nchan);            
        end
        switch options.position
            case 'Refit'
                rewPrintf('TNT: Refitting positions of tracks: constructing image\n');
                img = options.accumFunc(cacheORfile,{'tag'},options.framebinning,[],options.timegate); %TODO use cache if exists
                img = sum(img,3); % sum all channels
                halfw = max(round(3*quantile(trackingData(:,8),0.99)),4); % for dense data this window might be too large
                imgLoc = zeros(2*halfw+1,2*halfw+1,trackingData(end,1));
                imgInd = 0:2*halfw; % due to padding the coordinates are shifted by halfw
                % Pad image to avoid out of bounds
                img = padarray(img,[halfw,halfw]);

                assert(issorted(trackingData(:,1))); %Tracks need to be sorted but might have gaps in the numbering
                [track_ids,track_start,track_ind] = unique(trackingData(:,1),'first');
                track_len = diff([track_start; numel(track_ind)+1]);
                oldLoc = nan(numel(track_ids),6);
                offLoc = nan(numel(track_ids),2);% offset between imgLoc
                for ctrack = track_ids(:)'% loop over tracks
                    % Output process every 0.5 seconds
                    mainTime = toc(mainTime_start);
                    if( (mainTime-lastElapsedTime) > 0.5)
                        rewindMessages();
                        rewPrintf('TNT: Refitting positions of tracks: constructing image (%i/%i)\n',ctrack,track_ids(end));
                        lastElapsedTime = mainTime;
                    end
                    c_ind = track_start(ctrack)-1+(1:track_len(ctrack));
                    tempLoc = trackingData(c_ind,3:8);                           % x,y,z, Amp, BG, sigma
                    if track_len(ctrack)>1
                        tempLoc = [mean(tempLoc(:,1:3),1),sum(tempLoc(:,4:5),1),mean(tempLoc(:,6:end),1)];     % mean for pos/sigma, sum for Amp,BG
                    end
                    offLoc(ctrack,:) = round(tempLoc(1:2));
                    % crop fit window
                    imgLoc(:,:,ctrack) = sum(img(imgInd+offLoc(ctrack,2),imgInd+offLoc(ctrack,1),:,trackingData(c_ind,2)),4); % sum all frames. Gaps in the track are not filled.
                    % shift position to fit window
                    oldLoc(ctrack,:) = [(tempLoc(1:2) +1 -offLoc(ctrack,:) +halfw), tempLoc(3:end)];
                end
                clear img;
                % flatten to 2D image
                oldLoc(:,1) = oldLoc(:,1) + size(imgLoc,2) * (0:size(oldLoc,1)-1)';
                imgLoc = reshape(imgLoc,size(imgLoc,1),[]);
                % refit positions
                rewindMessages();
                rewPrintf('TNT: Refitting positions of tracks: performing fit\n');                
                newLoc = psfFit_Image(imgLoc,oldLoc(:,[1,2,4,5,6])',[1,1,1,1,1,0,0],true,true,halfw);          % 'optimize x,y,A,BG,sigma' (isoptric gaussian),integrated=true,MLE=true,halfw
                % Calculate positions back
                % newLoc = [xpos,ypos,A,BG,sigma_x,sigma_y,angle; exitflag]
                newLoc(1,:) = newLoc(1,:) - size(imgLoc,1) * (0:size(newLoc,2)-1);
                newLoc(1:2,:) = newLoc(1:2,:) -(1 -offLoc' +halfw);
                newLoc_valid = newLoc(end,:)>0;
                newLoc = [newLoc(1:2,:);zeros(1,size(newLoc,2));newLoc(3:4,:)./track_len';newLoc(5,:)]; % Add z, normalise amplitude and background to number of frames
                trackingData = [trackingData(:,1:2),newLoc(:,track_ind)',trackingData(:,9:end)];
                % Remove failed fits
                trackingData = trackingData(newLoc_valid(track_ind),:);
                clear imgLoc oldLoc newLoc offLoc track_ids track_start track_ind track_len;
            case 'Average'
                rewPrintf('TNT: Averageing positions of tracks\n');
                %for ctrack = trackingData(1,1):trackingData(end,1)% loop over tracks
                %    c_ind = trackingData(:,1)==ctrack;
                %    trackingData(c_ind,:) = [trackingData(c_ind,1:2),ones(sum(c_ind),1).*mean(trackingData(c_ind,3:5),1),trackingData(c_ind,6:end)]; % Take the mean of the position.
                %end
                % Equivalient to the loop but much faster
                [~,~,ind_accum] = unique(trackingData(:,1));
                posAVG = trackingData(:,3:5);
                [sub1, sub2] = ndgrid(ind_accum,1:size(posAVG,2));
                posAVG = accumarray([sub1(:) sub2(:)],posAVG(:))./accumarray(ind_accum,1);
                trackingData(:,3:5) = posAVG(ind_accum,:);
                clear sub1 sub2 ind_accum posAVG;
            case 'Preserve'
                % Do nothing
            otherwise
                warning off backtrace
                warning('Unrecognized position option. Fallback to ''Preserve''.');
                warning on backtrace
        end
        clear img;
        if isempty(trackingData)
            warning('trackingData empty.');
            return;
        end
        % get sigma
        if size(trackingData,2)<8
            warning('TNT: The plugin ''%s'' expects a localization sigma in the 8th column. Setting sigma to 1.3. \n',options.plugin_name);
            c_simga = 1.3;
        else
            c_simga = median(trackingData(:,8)); % Take the median sigma. If we are looking at something other than single molecules it would be necessary to take the fitted sigma for each.
        end
        for cframe = trackingData(1,2):trackingData(end,2)% loop over frames
            % Output process every 0.5 seconds
            mainTime = toc(mainTime_start);
            if( (mainTime-lastElapsedTime) > 0.5)
                rewindMessages();
                rewPrintf('TNT: Extracting TCSPC in frame %i/%i\n',cframe,trackingData(end,2));
                lastElapsedTime = mainTime;
            end
            
            % generate mask
            % the mask is padded by img_pad to avoid out of bound errors
            % for localisations at the edges. The padding is cut in the
            % last step.
            img_pad = 2;
            c_ind = trackingData(:,2)==cframe;
            % remove localisations out of the mask+padding (happens rarely when fitting with free sigma)
            c_pos = img_pad+round(trackingData(c_ind,3:4));
            c_pos_valid = all(c_pos>0,2) & all(c_pos<(imgSZ([2 1])+2*img_pad),2);
            c_pos = c_pos(c_pos_valid,:);
            c_ind(c_ind) = c_pos_valid;
            if ~any(c_ind)
                % empty frame
                continue;
            end
            masksz = ceil(c_simga*options.maskRadius);
            [x,y] = meshgrid(-masksz:masksz,-masksz:masksz);
            mask = ceil(sqrt(x.^2/masksz^2+y.^2/masksz^2))<2;
            
            % get TrackIDs
            if options.sumTCSPC
                track_ids = trackingData(c_ind,1);
            else
                track_ids = find(c_ind);
            end
            % TrackIDs are allways sorted, we avoid gaps in tcpsc_frame by using
            % 1:numel(track_ids) as temporary index and move it to position
            % track_ids in tcspc_mol.
            % construct mask (overlapping areas are removed from the mask)
            img_ind = zeros(imgSZ+2.*img_pad);
            img_ind(sub2ind(imgSZ+2.*img_pad,c_pos(:,2),c_pos(:,1))) = 1:numel(track_ids);
            % exclude pixels with overlaping molecules
            img_mask = conv2(img_ind>0,mask,'same')==1;
            img_ind = conv2(img_ind,mask,'same');
            img_ind(~img_mask) = 0;
            img_ind = img_ind(img_pad+(1:imgSZ(1)),img_pad+(1:imgSZ(2)));
            
            tcpsc_frame = options.accumFunc(cacheORfile,{'tcspc_pix'},...
                [options.framebinning(1) min(options.framebinning(2)+(cframe-[1 0])*min(realmax,options.framebinning(1))-[0 1],[inf options.framebinning(3)])],... %Calculate the first and last unbinned frame of cframe
                img_ind,false,true);
            
            % The last localisation(s) can be missing in tcpsc_frame if they
            % do not contain any photons. This can occur when there are
            % overlapping localisation from all directions.
            track_ids = track_ids(1:size(tcpsc_frame,2));
            tcspc_mol(1,track_ids,:,:) = tcspc_mol(1,track_ids,:,:) + cast(tcpsc_frame,class(tcspc_mol));
        end
        
        cacheORfile = options.filename_movie; %#ok<NASGU> cache no longer needed.
        % calculate fast LT
        resolution = head.MeasDesc_Resolution*1e9; % in ns
        
%         t = shiftdim((1:size(tcspc_mol,3))',-2)*resolution; % 1D time axis, but in the 3rd dimension.
%         tau_fast = real(sqrt(((squeeze(sum(t.^2.*double(tcspc_mol),3)))-((squeeze(sum(t.*double(tcspc_mol),3))).^2)./nphot)./(nphot-1))); % Normalise with (N-1) to make it consitent with (my)var
        
        tcspc_mol = permute(tcspc_mol(:,:,:,1),[2 3 4 1]);%[track gate chan], just use the first channel for now
        
        nphot = sum(tcspc_mol,2);
        t = (1:size(tcspc_mol,2))*resolution; % 1D time axis in the 2nd dimension.
        tau_fast = real(sqrt(((squeeze(sum(t.^2.*double(tcspc_mol),2)))-((squeeze(sum(t.*double(tcspc_mol),2))).^2)./nphot)./(nphot-1))); % Normalise with (N-1) to make it consitent with (my)var
        
        % fit
        
        switch options.fitType
            case {'monoexponential MLE','monoexponential PM'}
                
                [tcspc_cut,t,tail_ind] = tcspc_apply_cutoff(tcspc_mol,options.cutoff,resolution);
                
                if strcmp(options.fitType,'monoexponential PM')
                    n_lt = 500; % number of possible lifetime values
                    n_bg = 50;  % number of possible background values
                else
                    n_lt = 100;
                    n_bg = 20;
                end
                    
                taus_pm = linspace(max(options.minLT,eps),min(options.maxLT,t(end)),n_lt);
                bgs_pm = linspace(0,1-1/n_bg,n_bg-1).^2; % square linspace to sample low bg values more densely
                
                [taus_pm,bgs_pm] = meshgrid(taus_pm(:),bgs_pm(:)');
                taus_pm = taus_pm(:);
                bgs_pm = bgs_pm(:);
                taus_pm(end+1) = nan;
                bgs_pm(end+1) = 1;
                
                rewindMessages();
                rewPrintf('TNT: PM: Calculating decay patterns\n');
                t_delta = mean(diff(t));
                pfun_monoexp = @(tau,T,t) 1./tau.*exp(-t./tau)./(1-exp(-T./tau)) * t_delta; % Normalised monoexponetial decay
                pfun_monoexpBG = @(tau,b)b./numel(t)+(1-b).*pfun_monoexp(tau,t(end),t(:)');
                monoexp_logdecays = log(pfun_monoexpBG(taus_pm,bgs_pm));
                
                bunch_sz = 1e4;
                PM_ind = nan(size(tcspc_cut,1),1);
                for idx = 1:ceil(numel(PM_ind)/bunch_sz)
                    ind = ((idx-1)*bunch_sz+1):min(idx*bunch_sz,numel(PM_ind));
                    [~,PM_ind(ind)] = max(tcspc_cut(ind,:) * monoexp_logdecays',[],2);
                    rewindMessages();
                    rewPrintf('TNT: PM: Matching decay patterns %i/%i\n',ind(end),size(tcspc_cut,1));
                end
                
                clear monoexp_logdecays;
                
                PM_tau = taus_pm(PM_ind);
                PM_bg = bgs_pm(PM_ind);
                PM_np = sum(tcspc_cut,2);
                
                PM_tau(PM_np==0) = nan;
                PM_bg(PM_np==0) = nan;
                
                if strcmp(options.fitType,'monoexponential PM')
                    rewindMessages();
                    rewPrintf('TNT: PM: Calculate chi2\n');
                    [PM_ind_uq,~,PM_ind_uqx] = unique(PM_ind);
                    monoexp_decays = pfun_monoexpBG(taus_pm(PM_ind_uq),bgs_pm(PM_ind_uq));
                    % Calculate reduced chi2
                    PM_chi2 = sum((tcspc_cut-PM_np.*monoexp_decays(PM_ind_uqx,:)).^2./monoexp_decays(PM_ind_uqx,:),2)./PM_np./(numel(t)-3);
                    
                    clear monoexp_decays;
                    %% Save
                    options.outParamDescription = [options.outParamDescription;{'lt-tau';'lt-np';'lt-bg';'lt-chi2'}];
                    fitData = [PM_tau(:),PM_np(:),PM_bg(:),PM_chi2(:)];
                else
                    fitData = nan(size(tcspc_mol,1),4);
                    options.outParamDescription = [options.outParamDescription;{'lt-tau';'lt-amp';'lt-bg';'lt-chi2'}];
                    tol = options.tolerance; % 1e-8 default
                    steps = 600;
                    
                    attmpt = options.attempts;
                    attmpt_mod = 1+(-(attmpt-1)/2:(attmpt-1)/2)/attmpt; % Slightly vary tau_init when performing multiple attemts
                    attmpt_p = nan(3,attmpt);
                    attmpt_mle = nan(1,attmpt);
                    
                    PM_tau(isnan(PM_tau)) = nanmean(PM_tau); % Replace nans
                    amp_init = (1-PM_bg).*PM_np ./ exp(-options.cutoff./PM_tau);
                    bg_init = max(PM_bg.*PM_np,1);
                    for idx = 1:size(tcspc_cut,1)
                        % Output process every 0.5 seconds
                        mainTime = toc(mainTime_start);
                        if( (mainTime-lastElapsedTime) > 0.5)
                            rewindMessages();
                            rewPrintf('TNT: Fitting TCSPC %i/%i\n',idx,size(tcspc_cut,1));
                            lastElapsedTime = mainTime;
                        end
                        
                        if nphot(idx)<options.minAmp
                            fitData(idx,1:4) = NaN;
                        else
                            % Performe requested number of attempts
                            for attmpt_idx = 1:attmpt
                                % MLE Fit using distfit as initial guess
                                [attmpt_p(:,attmpt_idx),~,~,attmpt_mle(attmpt_idx)] = Simplex_Handle(@ExpFunMono, [bg_init(idx),PM_tau(idx)*attmpt_mod(attmpt_idx),amp_init(idx)], [0 options.minLT 0],[inf options.maxLT inf],tol,-steps,t,tcspc_cut(idx,:)',true);
                            end
                            [~, attmpt_best] = min(attmpt_mle);
                            [lt_err] = ExpFunMono(attmpt_p(:,attmpt_best),t,tcspc_cut(idx,:)',false);                                                                                                 % chi2 for output
                            fitData(idx,1:4) = [attmpt_p(2,attmpt_best), attmpt_p(3,attmpt_best), attmpt_p(1,attmpt_best), lt_err];
                        end
                    end
                end
                
            case {'average DistFit'}
                fitData = nan(size(tcspc_mol,1),4);
                options.outParamDescription = [options.outParamDescription;{'lt-tau';'lt-amp';'lt-bg';'lt-chi2'}];
                
                [tcspc_cut,t,tail_ind] = tcspc_apply_cutoff(tcspc_mol,options.cutoff,resolution);
                taus_guess = logspace(log10(options.minLT),log10(min(options.maxLT,t(end))),10);
                
                for idx = 1:size(tcspc_cut,1)
                    % Output process every 0.5 seconds
                    mainTime = toc(mainTime_start);
                    if( (mainTime-lastElapsedTime) > 0.5)
                        rewindMessages();
                        rewPrintf('TNT: Fitting TCSPC %i/%i\n',idx,size(tcspc_cut,1));
                        lastElapsedTime = mainTime;
                    end
                    
                    if nphot(idx)<options.minAmp
                        fitData(idx,1:4) = NaN;
                    else
                        %% Make a DistFit and take weighted average
                        [~,c]   = ExpFunLSQ(taus_guess,t,tcspc_cut(idx,:)',false);
                        tau_dist = sum(c(2:end))/sum(c(2:end)'./taus_guess);
                        % Refine guess
                        [lt_err,c_ref] = ExpFunLSQ(tau_dist,t,tcspc_cut(idx,:)',false);
                        bg_dist = c_ref(1); %in photons, ensures that we are not starting at the boundry
                        amp_dist = c_ref(2); % in photons (analytically normalised, incluedes photons cut in the front)

                        fitData(idx,1:4) = [tau_dist, amp_dist, bg_dist, lt_err];
                    end
                end
            case 'maximum correlation'
                rewindMessages();
                rewPrintf('TNT: Correlating TCSPCs \n');
                % Exclude tcspc with less than minAmp photons
                fitind = nphot>=options.minAmp;
                % Cut tcspc and genetate time axis
                [tcspc_cut,t,tail_ind] = tcspc_apply_cutoff(tcspc_mol(fitind,:),options.cutoff,resolution);
                % Generate log spaced candidate lifetimes and calculate the
                % decay curves. Add constant as lower thershold for
                % correlation. (constant = 0 correlation)
                taus_guess = logspace(log10(options.minLT),log10(min(options.maxLT,t(end))),500);
                taus_z = [ones(numel(t),1), exp(-t(:)./taus_guess)];% Normalization does not matter for correlation
                taus_guess_bg = [NaN taus_guess];
                % Correlate each tcspc with all decays and select the decay
                % time with the highest correlation
                allcorr = corr(tcspc_cut',taus_z);
                [allcorr_max,allcorr_maxind] = max(allcorr,[],2);
                tau_corr = taus_guess_bg(allcorr_maxind);
                
                %% Save
                fitData = nan(size(tcspc_mol,1),2);
                options.outParamDescription = [options.outParamDescription;{'lt-tau';'lt-corr'}];
                fitData(fitind,1:2) = [tau_corr(:), allcorr_max(:)];
                
            otherwise
                if ~strcmpi(options.fitType,'fast lifetime')
                    warning off backtrace
                    warning('Unrecognized fit type. Switching to fast lifetime.');
                    warning on backtrace
                end
                % Recalculate fast lifetime with cut TCSPC
                [tcspc_cut,t,tail_ind] = tcspc_apply_cutoff(tcspc_mol,options.cutoff,resolution);
                nphot_cut = sum(tcspc_cut,2);
                tau_fast_cut = real(sqrt(((squeeze(sum(t.^2.*double(tcspc_cut),2)))-((squeeze(sum(t.*double(tcspc_cut),2))).^2)./nphot_cut)./(nphot_cut-1))); % Normalise with (N-1) to make it consitent with (my)var
                fitData = tau_fast_cut(:);
                options.outParamDescription = [options.outParamDescription;{'lt-tau'}];
        end
        rewindMessages();
        % save results
        if options.sumTCSPC
            if isempty(fitData)
                postprocData = [trackingData, tau_fast(trackingData(:,1),1), nphot(trackingData(:,1),1)];
            else
                postprocData = [trackingData, tau_fast(trackingData(:,1),1), nphot(trackingData(:,1),1), fitData(trackingData(:,1),:)];
            end
        else
            postprocData = [trackingData, tau_fast(:), nphot(:), fitData];
        end
        
        if options.exportTCSPC
            options.TCSPC = tcspc_mol;
            options.TCSPC_t = t;
            options.TCSPC_tail = tail_ind;
        end
    else
        warning('TNT: The plugin ''%s'' only supports single photon data. \n',options.plugin_name);
    end
    %% Local functions for status printing
    function rewPrintf(msg, varargin)
        % Rewindable message printing: Print msg and cache it.
        % Usage is analogous to sprintf.
        msg = sprintf(msg, varargin{:});
        msgAccumulator = [msgAccumulator, msg];
        fprintf(msg);
    end

    function rewindMessages()
        % Remove cached messages from command line, reset cache
        reverseStr = repmat(sprintf('\b'), 1, length(msgAccumulator));
        fprintf(reverseStr);
        
        msgAccumulator = '';
    end
end

% Cuts the tcspc cutoff bins/ns after the position off the maximum in the sum
% of all curves with a intensity in the [0.25 0.75] quantile. This excludes
% curves potentially containing dead time artefacts or only background. Exspect
% tcspc in the form [mol gate chan]. If cutoff is in ns resolution has to be
% given. As second output a cut time axis centred at the peak is returned.
function [tcspc,t,ind] = tcspc_apply_cutoff(tcspc,cutoff,resolution)
    if nargin<2 && ~isempty(resolution)
        resolution = 1;
    end
    cutoff = round(cutoff/resolution);
    tcspc_int = sum(sum(tcspc,3),2);
    int_thres = quantile(tcspc_int,[0.25 0.75]);
    [~,max_pos] = max(sum(sum(tcspc(int_thres(1) <= tcspc_int & tcspc_int <= int_thres(2),:,:),3),1));
    ind = (max_pos+cutoff):size(tcspc,2);
    tcspc = tcspc(:,ind,:);
    t = ((1:size(tcspc,2))+cutoff-1)*resolution;
end
