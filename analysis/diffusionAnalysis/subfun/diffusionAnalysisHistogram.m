function [msd_result, velocity_result, plotdata_result,use_v] = diffusionAnalysisHistogram(tracks_processed,fitParam,experimParam)
% [msd_result, velocity_result, plotdata_result] =
% diffusionAnalysisHistogram(tracks_processed,fitParam,experimParam)
% Perform displacement histogram analysis of trajectories and calculate MSD
% curves.
% 
% INPUT:
%     tracks_processed: Cell row array of trajectories with one trajectory
%     row array per cell. Trajectory segments have to be separated by a NaN
%     line.
%     
%     fitParam: Struct of fitting options, see RunDiffusion Analysis.m for
%     details.
%     
%     experimParam: Struct of experimental parameters, see RunDiffusion
%     Analysis.m for details.
%     
%     
% OUTPUT:
%     msd_result: Cell array of msd curves. Every cell contains a double
%     array [t,sigma_i,error_sigma_i,weight_i,error_weight_i,outcome].
%     Here, t is the time in [s], sigma are the MSD values of each
%     diffusing subspecies in [µm^2] and weight are the weights.
% 
%     The row marks the number of subspecies and the column the dimension.
%     For example, msd_result(3,2) contains the MSD values for 3 diffusing
%     species for the y dimension (second dimension). Usually, the number
%     of possible subspecies tends to decrease for higher times as all fast
%     species vanish due to absent trajectories. The variable "outcome"
%     then marks if the outcome of the fit for these particular MSD values
%     at this particular time point is the best possible one (true) or if a
%     better one (higher number of species) exists (false).
%     One can then choose to only count "true" fit values when fitting an MSD curve.
%     
%     velocity_result: Double array of velocities
%     [[v_x,v_y];[error_v_x,error_v_y]] in [µm/s] taken as the unweighted
%     mean over all fits.
%     
%     plotdata_result: Struct array containing double arrays
%     [delta_r,histogram,fit] of histograms for plotting in the units
%     [px,1/px,1/px]. The struct has two entries, plotdata_x and plotdata_y
%     for both possible dimensions. For non-isotropicity, plotdata_y is
%     empty.
% 
%     use_v: Boolean switch for velocity fit, will be disabled
%     automatically for non appropriate inputs.


% Parse initial parameters
px = experimParam.pixel/experimParam.magnification;
dt = experimParam.timestep;
D_initial = experimParam.diffusionCoeffInitial/px^2*dt; %diffusion coefficient in pixel^2/frame
max_frame = fitParam.maximumSkip;
use_corr = fitParam.correlated;
fit_curve = fitParam.distanceMethod;
hist_bin = fitParam.histogramBinWidth;
use_v = fitParam.fitVelocity;
use_iso = fitParam.isotropic;
plot_fit = fitParam.plotFit;
bcond = fitParam.breakCondition;
alpha = fitParam.fitAlpha;

order = numel(D_initial);


if ~use_iso && strcmp(fit_curve,'jump');
    fprintf('Cannot consider non-isotropicity for jump distance histogram. Aborting. \n');
    return;
end

if use_v
    if strcmp(fit_curve,'jump')
        fprintf('Cannot fit velocity with jump distance histogram. Disabling velocity fit. \n');
        use_v = false;
    end
    
    if use_iso
        fprintf('Cannot fit velocity for one dimension. Disabling velocity fit. \n');
        use_v = false;
    end
end

% Start main loop
fit_result = repmat(struct('p_x',[],'p_y',[],'order_x',[],'order_y',[]),max_frame,1);
plotdata_result = repmat(struct('plotdata_x',[],'plotdata_y',[]),max_frame,1);
fit_aborted = false;
order_max = 1;

for iFrame=1:max_frame
    
    %format all trajectories with appropriate skip frame length
    [traj] = tracksSkip(tracks_processed,iFrame,use_corr);
    
    if size(traj,1)<bcond(1) %break if number of displacements is too low
        fit_aborted = true;
        break
    end
    
    %create displacement histograms
    if strcmp(fit_curve,'jump')
        data_hist_in = sqrt(traj(:,1).^2+traj(:,2).^2);
        [bins_out,hist_out,~] = Hist1D(data_hist_in,hist_bin,false,'freedman',true,false,false);
        hist_out = hist_out/sum(hist_out*(bins_out(2)-bins_out(1)));
    else
        if use_iso
            data_hist_in = [traj(:,1);traj(:,2)];
            [bins_out,hist_out,~] = Hist1D(data_hist_in,hist_bin,false,'freedman',true,true,false);
            hist_out = hist_out/sum(hist_out*(bins_out(2)-bins_out(1)));
        else
            [bins_out,hist_out,~] = Hist1D(traj(:,1),hist_bin,false,'freedman',true,true,false);
            hist_out = hist_out/sum(hist_out*(bins_out(2)-bins_out(1)));
            [bins_out_y,hist_out_y,~] = Hist1D(traj(:,2),hist_bin,false,'freedman',true,true,false);
            hist_out_y = hist_out_y/sum(hist_out_y*(bins_out_y(2)-bins_out_y(1)));
        end
    end
    
    if length(bins_out)<bcond(2) %break if number of histogram bins is too small
        fit_aborted = true;
        break
    elseif exist('bins_out_y','var') && length(bins_out_y)<bcond(2)
        fit_aborted = true;
        break
    end
    
    
    if use_v
        p0 = [ones(order-1,1)/order;2*D_initial(:)*iFrame;0]; %parameter vector: weight, sigma^2, velocity. We only need N-1 weights as sum w_i = 1
    else
        p0 = [ones(order-1,1)/order;2*D_initial(:)*iFrame];
    end
    
    % Fit histogram
    [p,order_new,plotdata] = fitDisplacementHistogram(bins_out,hist_out,p0,use_v,order,fit_curve,alpha,plot_fit,iFrame,px,dt); % Fit Gaussian to displacement histogram
    fit_result(iFrame).p_x = p; 
    fit_result(iFrame).order_x = order_new;
    plotdata_result(iFrame).plotdata_x = plotdata;
    
    if order_max<order_new
        order_max = order_new;
    end
    
    if ~use_iso
        if use_v
            p0 = [ones(order-1,1)/order,2*D_initial(:)*iFrame;0]; %parameter vector: weight, sigma^2. We only need N-1 weights as sum w_i = 1
        else
            p0 = [ones(order-1,1)/order,2*D_initial(:)*iFrame];
        end
        
        [p,order_new,plotdata] = fitDisplacementHistogram(bins_out_y,hist_out_y,p0,use_v,order,fit_curve,alpha,plot_fit,iFrame,px,dt); % Fit Gaussian to displacement histogram
        fit_result(iFrame).p_y = p;
        fit_result(iFrame).order_y = order_new;
        plotdata_result(iFrame).plotdata_y = plotdata;
        
        if order_max<order_new
            order_max = order_new;
        end
    end
end

if fit_aborted
    iFrame = iFrame-1;
end
fit_result = fit_result(1:iFrame);
plotdata_result = plotdata_result(1:iFrame);

[msd_result,velocity_result] = prepareResult(fit_result,use_v,use_iso,order_max,px,dt);
end %END MAINFUN



function [msd_result,velocity_result] = prepareResult(fit_result,use_v,use_iso,order_max,px,dt)
max_frame = size(fit_result,1);

%velocity
velocity_result = [];
if use_v
    v_x = zeros(2,max_frame);
    v_y = v_x;
    for iFrame=1:max_frame
        %choose velocity from final/best fit --> p_x{end}
        v_x(1+(iFrame-1)*2:iFrame*2) = [fit_result(iFrame).p_x{end}(1,end);fit_result(iFrame).p_x{end}(2,end)]/iFrame*px/dt; %normalize by amount of frames
        if ~use_iso
            v_y(1+(iFrame-1)*2:iFrame*2) = [fit_result(iFrame).p_y{end}(1,end);fit_result(iFrame).p_y{end}(2,end)]/iFrame*px/dt;
        end
    end
    v_x = [mean(v_x(1,:));sqrt(sum(v_x(2,:).^2)/max_frame^2+std(v_x(1,:))/max_frame)]; %avg(v_i) +/- sqrt(sigma_sys^2+sigma_stat^2)
    if ~use_iso
        v_y = [mean(v_y(1,:));sqrt(sum(v_y(2,:).^2)/max_frame^2+std(v_y(1,:))/max_frame)];
    else
        v_y = [];
    end
    
    velocity_result = [v_x,v_y]; %[[v_x,v_y];[error_v_x,error_v_y]]
end

%msd curve
% msd_result containts one cell per order. in each cell, there is a list [t,sigma_i,error_sigma_i,weight_i,error_weight_i,...,switch]. "switch" signifies if this specific combination of diffusion coefficient for this specific timepoint was the correct one
if ~use_iso
    msd_result = cell(order_max,2);
else
    msd_result = cell(order_max,1);
end

for iFrame=1:max_frame
    p_x = fit_result(iFrame).p_x;
    order = numel(p_x);
    for jOrder = 1:order
        if ~isempty(p_x{jOrder})
            msd_part = prepareResultMSD(p_x,iFrame,jOrder,order,px,dt);
            msd_result(jOrder,1) = {[msd_result{jOrder,1};msd_part]};
        else
            msd_result(jOrder,1) = {[]};
        end
    end
    
    if ~use_iso
        p_y = fit_result(iFrame).p_y;
        order = numel(p_y);
        for jOrder = 1:order
            if ~isempty(p_y{jOrder})
                msd_part = prepareResultMSD(p_y,iFrame,jOrder,order,px,dt);
                msd_result(jOrder,2) = {[msd_result{jOrder,2};msd_part]};
            else
                msd_result(jOrder,2) = {[]};
            end
        end
    end
end %iFrame

end %prepareResult



function [msd_part] = prepareResultMSD(p,iFrame,jOrder,order,px,dt)

if jOrder==1
    sigsq = p{jOrder}(:,jOrder:2*jOrder-1).'*px^2;
    weights = [1,0];
    msd_part = [iFrame*dt,sigsq,weights,jOrder==order];
else
    sigsq = p{jOrder}(:,jOrder:2*jOrder-1)*px^2; %[[sigma_1;error_sigma_1],...]
    weights = abs(p{jOrder}(:,1:jOrder-1)); weights = [weights,[1-sum(abs(weights(1,:)));sqrt(sum(weights(2,:).^2)/jOrder+std(weights(1,:))^2/jOrder)]]; %[[w_1;error_w_1],...]
    [~,idx] = sort(sigsq(1,:),'ascend');
    sigsq = sigsq(:,idx); %sort, starting at lowest diffusion coefficient
    weights = weights(:,idx);
    
    msd_part = [iFrame*dt,zeros(1,4*jOrder),jOrder==order]; %[t,sigma1,error_sigma1,weight1,error_weight1,...,is this correct order?]
    msd_part(2:4:end-1) = sigsq(1,:);
    msd_part(3:4:end-1) = sigsq(2,:);
    msd_part(4:4:end-1) = weights(1,:);
    msd_part(5:4:end-1) = weights(2,:);
end

end
