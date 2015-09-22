function [bleachData,blinkData,intensityData] = surfaceFluorophoreAnalysis(dataObject,timestep,enablePlot)
% [bleachData,blinkData,intensityData] = surfaceFluorophoreAnalysis(filename,cameraOptions,enablePlot)
% Perform analysis of fluorescent emitters on a surface, returning
% bleaching lifetime, blinking lifetime, average intensity in time, and
% intensity distribution
%
% INPUT:
%     dataObject: Can be either full path to TrackNTrace file to analyze OR trajectoryData array which comes out of TrackNTrace.
%
%     timestep: Camera acquisition time in seconds.
%
%     enablePlot: boolean, enables plotting and saving plots
%
% OUTPUT:
%     bleachData: Struct containing fitted bleach lifetime curves and fit
%     parameters.
%         -.distribution: Nx2 double array where N is the number of data
%         points. x and y axis (time and number) are stored in columns
%         -.residuals: Nxdim double array where dim is the number of
%         exponential functions fitted. Obtain fit curves by adding
%         respective residuals to distribution. -.param: Cell array of
%         fitted params, one per multi-exponential fit. [N1,tau1,B] for
%         mono-exponential, [N1,tau1,N2,tau2,B] for bi-exponential fit and
%         so on. 
%         -.framesDistribution: Nx3 double array with columns time,
%         number of total molecules, number of off molecules
%
%     blinkData: Struct containing fitted blinking lifetime curves and fit
%     parameters. Organized as bleachData.
%
%     intensityData: Struct containing intensity distribution (distribution
%     of intensities, averaged per trajectory) and trace (intensity trace
%     averaged over all molecules in respective frame w.r.t. first frame of
%     occurence, normalized to maximum intensity per trajectory). Also
%     contains distributions for total intensity per trajectory, background
%     and signal-to-background.
%

%% Parse some important options
if ischar(dataObject)
    filename = dataObject;
    [folder,file,~] = fileparts(filename);
    filename_base = [folder,filesep,file];
    saveToFile = true;
    load(filename,'trajectoryData');
else
    time = clock;
    timestamp = sprintf('%i-m%02i-d%02i-%ih%i',time(1),time(2),time(3),time(4),time(5));
    filename_base = ['surfaceAnalysis_',timestamp];
    trajectoryData = dataObject;
    saveToFile = false;
end
        
dt = timestep;

%load file

nrFrames = max(trajectoryData(:,2)); 
nrTraj = max(trajectoryData(:,1));

%Integrate intensities (factory 2pi sigma^2)
if size(trajectoryData,2)==7
    trajectoryData(:,5) = trajectoryData(:,5)*2*pi.*trajectoryData(:,7).^2;
end

%% Intitialize arrays
molecules_in_frame_relative = zeros(nrFrames,1); % number of molecules in relative image frame
molecules_in_frame = zeros(nrFrames,1); % number of molecules in absolute image frame

blink_distr = zeros(nrFrames,1); %collection of all lengths of blink off times (triplet state lifetime)
blink_in_frame = zeros(nrFrames,1); % number of emitters in off-state per frame

intensity_total_distr = zeros(nrTraj,1); %array of all summed up trajectory intensities (integrated)
intensity_trace_mean = zeros(nrFrames,1); %average intensity trace relative to first frame of appearance
intensity_trace_mean_counter = zeros(nrFrames,1); %number of traces contributing to respective element in avg trace

%% Collect data from trajectories
for iTraj = 1:nrTraj
    %Get current trajectory, increment active molecules counter
    trajectory = trajectoryData(trajectoryData(:,1)==iTraj,:);
    traj_start_frame = trajectory(1,2);
    traj_end_frame = trajectory(end,2);
    molecules_in_frame(traj_start_frame:traj_end_frame) = molecules_in_frame(traj_start_frame:traj_end_frame)+1;
    
    %Get active molecules relative to frame of first occurence
    trajectory(:,2) = trajectory(:,2)-traj_start_frame+1;
    molecules_in_frame_relative(1:trajectory(end,2)) = molecules_in_frame_relative(1:trajectory(end,2))+ones(trajectory(end,2),1);
    
    %Get blink histogram
    nr_blink_frames = trajectory(2:end,2)-trajectory(1:end-1,2)-1; %all gaps show up as a difference >1, subtract 1 and we get amount of inactive frames
    inactive_frames = (traj_start_frame:traj_end_frame).'; inactive_frames(trajectory([true;nr_blink_frames==0],2)) = []; %array of inactive frames. first frame is always active
    blink_in_frame(inactive_frames) = blink_in_frame(inactive_frames)+1; %increment blink count histogram
    nr_blink_frames = nr_blink_frames(any(nr_blink_frames,2)); %only get blink frames
    blink_distr(nr_blink_frames+1) = blink_distr(nr_blink_frames+1)+1; %increment histogram
    
    %Get intensity statistics
    intensity_total_distr(iTraj) = sum(trajectory(:,5));
    trajectory(:,5) = trajectory(:,5)/max(trajectory(:,5));
    intensity_trace_mean(trajectory(:,2)) = intensity_trace_mean(trajectory(:,2))+trajectory(:,5);
    intensity_trace_mean_counter(trajectory(:,2)) = intensity_trace_mean_counter(trajectory(:,2))+ones(size(trajectory,1),1);
end

%% Create data for analysis and plots

%Bleaching curve
t0 = find(diff([molecules_in_frame_relative(2:end),molecules_in_frame_relative(1:end-1)],1,2)>5); t0 = t0(1); %drop of point for bleaching molecules
bleach_taxis = ((0:(numel(molecules_in_frame_relative(t0:end))-1))*dt).'; %time axis, starts at 0
[bleach_param,bleach_res] = fitExponential(bleach_taxis,molecules_in_frame_relative(t0:end),true);
bleach_distr = [bleach_taxis,molecules_in_frame_relative(t0:end)];

%blink curve
blink_taxis = 0:nrFrames-1;
idx_blink = blink_distr>0;
blink_taxis = blink_taxis(idx_blink).'*dt;
[blink_param,blink_res] = fitExponential(blink_taxis,blink_distr(idx_blink),false);
blink_distr = [blink_taxis,blink_distr(idx_blink)];
blink_in_frame = [(0:nrFrames-1).'*dt,molecules_in_frame,blink_in_frame];

%Intensity curves  - we could use intensity filters for this!
[intensity_xbins,intensity_distr,~] = Hist1D(trajectoryData(:,5),0,false,'freedman',true,false,true);
intensity_distr = [intensity_xbins(1:end-1).',intensity_distr(1:end-1).'];

[intensity_total_xbins,intensity_total_distr,~] = Hist1D(intensity_total_distr,0,false,'freedman',true,false,true);
intensity_total_distr = [intensity_total_xbins(1:end-1).',intensity_total_distr(1:end-1).'];

[bg_xbins,bg_distr,~] = Hist1D(trajectoryData(:,6),0,false,'freedman',true,false,true);
bg_distr = [bg_xbins(1:end-1).',bg_distr(1:end-1).'];

sbr_distr = trajectoryData(:,5)./trajectoryData(:,6);
[sbr_xbins,sbr_distr,~] = Hist1D(sbr_distr(sbr_distr<100*median(sbr_distr)),0,false,'freedman',true,false,true);
sbr_distr = [sbr_xbins(1:end-1).',sbr_distr(1:end-1).'];

intensity_trace_mean = intensity_trace_mean./intensity_trace_mean_counter;

%% Save (and plot)
bleachData.distribution = bleach_distr;
bleachData.param = bleach_param;
bleachData.residuals = bleach_res;

blinkData.distribution = blink_distr;
blinkData.param = blink_param;
blinkData.residuals = blink_res;
blinkData.framesDistribution = blink_in_frame;

intensityData.distribution = intensity_distr;
intensityData.total_distribution = intensity_total_distr;
intensityData.background_distribution = bg_distr;
intensityData.sbr_distribution = sbr_distr;
intensityData.trace = intensity_trace_mean;
if saveToFile
    save(filename,'bleachData','blinkData','intensityData','-append');
end

if enablePlot
    clearvars -except bleachData blinkData intensityData filename_base
    
    filename_plot = [filename_base,'_bleach.png'];
    plotData(bleachData,'bleach');
    export_fig(filename_plot,'-png','-transparent','-r300');
    close all force
    
    filename_plot = [filename_base,'_blink.png'];
    plotData(blinkData,'blink');
    export_fig(filename_plot,'-png','-transparent','-r300');
    close all force
    
    filename_plot = [filename_base,'_blink_frames.png'];
    plotData(blinkData,'blink_frames');
    export_fig(filename_plot,'-png','-transparent','-r300');
    close all force
    
    filename_plot = [filename_base,'_int_distr.png'];
    plotData(intensityData,'intensity_distr');
    export_fig(filename_plot,'-png','-transparent','-r300');
    close all force
    
    filename_plot = [filename_base,'_int_total_distr.png'];
    plotData(intensityData,'intensity_total_distr');
    export_fig(filename_plot,'-png','-transparent','-r300');
    close all force
    
    filename_plot = [filename_base,'_bg_distr.png'];
    plotData(intensityData,'bg_distr');
    export_fig(filename_plot,'-png','-transparent','-r300');
    close all force
    
    filename_plot = [filename_base,'_sbr_distr.png'];
    plotData(intensityData,'sbr_distr');
    export_fig(filename_plot,'-png','-transparent','-r300');
    close all force
    
    filename_plot = [filename_base,'_int_trace.png'];
    plotData(intensityData,'intensity_trace');
    export_fig(filename_plot,'-png','-transparent','-r300');
    close all force
end


end %MAINFUN


function [param,res] = fitExponential(x_fit,y,fit_short)

% Initialize arrays, parse inputs
t_step = x_fit(2)-x_fit(1);
t_half = find(y<=y(1)/2); %half life for mono-decay
p0_exp1 = [y(1),t_half(1)*t_step/log(2),y(end)]; %guess for mono-decay
p0_exp2 = [y(1)/2,t_half(1)*t_step/log(2),y(1)/2,t_half(1)*t_step/log(2)*5,y(end)]; %guess for bi-decay

options = optimset('Jacobian','off',...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6); %options for fitter

% Fit mono- bi-, and tri-exponential decay
[p_exp1,~,res_exp1,~] = lsqnonlin(@fitExpRes,p0_exp1,[],[],options,...
    y,x_fit);
[p_exp2,~,res_exp2,~] = lsqnonlin(@fitExpRes,p0_exp2,[],[],options,...
    y,x_fit);
param = [{p_exp1};{p_exp2}];
if fit_short
    p0_exp3 = [p_exp2(1:4),10,0.025,p_exp2(end)]; %short glass lifetime
    [p_exp3,~,res_exp3,~] = lsqnonlin(@fitExpRes,p0_exp3,[],[],options,...
        y,x_fit);
    param = [param;{p_exp3}];
else
    res_exp3 = [];
end

res = [res_exp1,res_exp2,res_exp3];
end %fitExponential


function [F] = fitExpRes(p,y,x_fit)
F = zeros(numel(y),1);

for iExp = 1:(numel(p)-1)/2
    F = F+p(1+(iExp-1)*2)*exp(-x_fit/p(iExp*2));
end
F = F+p(end)-y;

end %fitExpRes


function [] = plotData(data_struct,plot_type)

switch plot_type
    case {'bleach', 'blink'}
        plot_short = numel(data_struct.param)==3;
        distr = data_struct.distribution;
        param = data_struct.param;
        res = data_struct.residuals;
        
        
        legend_strings = [{'Mono-exp.'},{'Bi-exp.'}];
        
        subplot(6,1,1:5);
        plot(distr(:,1),distr(:,2)+res(:,1),'b-',distr(:,1),distr(:,2)+res(:,2),'g.-'); hold on;
        if plot_short
            plot(distr(:,1),distr(:,2)+res(:,3),'r--');
            legend_strings = [legend_strings,{'Tri-exp'}];
        end
        plot(distr(:,1),distr(:,2),'ko');
        hold off;
        legend_strings = [legend_strings,{'Data'}];
        
        xlabel('t [s]');
        if plot_short
            ylabel('No. of emitters');
        else
            ylabel('No. of emitters in off-state');
        end
        legend(gca,legend_strings,'Location','North');
        for i=1:2+plot_short
            text_string = ['tau_{',int2str(i),'} = '];
            for j=1:i-1
                text_string = [text_string, num2str(param{i}(2*j),3), ', ']; %#ok<AGROW>
            end
            text_string = [text_string, num2str(param{i}(2*i),3), ' s'];  %#ok<AGROW>
            text(0.6,0.65-i*0.05,text_string,'Units','normalized');
        end
        if plot_short
            title('Emitter bleaching lifetime');
        else
            title('Blinking state lifetime');
        end
        
        h_res = subplot(6,1,6);
        plot(distr(:,1),res(:,1),'b-',distr(:,1),res(:,2),'g.-'); hold on;
        if plot_short
            plot(distr(:,1),res(:,3),'r--');
        end
        hold off;
        set(h_res,'XTick',[]);
        
    case 'intensity_distr'
        bar(data_struct.distribution(:,1),data_struct.distribution(:,2));
        xlabel('Photons');
        ylabel('Frequency');
        title('Molecule intensity distribution');
        
    case 'intensity_total_distr'
        bar(data_struct.total_distribution(:,1),data_struct.total_distribution(:,2));
        xlabel('Photons');
        ylabel('Frequency');
        title('Trajectory sum intensity distribution');
        
    case 'bg_distr'
        bar(data_struct.background_distribution(:,1),data_struct.background_distribution(:,2));
        xlabel('Photons');
        ylabel('Frequency');
        title('Molecule background distribution');
        
    case 'sbr_distr'
        bar(data_struct.sbr_distribution(:,1),data_struct.sbr_distribution(:,2));
        xlabel('SBR [-]');
        ylabel('Frequency');
        title('Molecule signal-to-background distribution');
        
    case 'intensity_trace'
        plot(1:numel(data_struct.trace), data_struct.trace,'b-');
        xlabel('Frame');
        ylabel('Normalized mean intensity');
        
    case 'blink_frames'
        plot(data_struct.framesDistribution(:,1),data_struct.framesDistribution(:,3)./data_struct.framesDistribution(:,2),'ko');
        xlabel('t [s]');
        ylabel('Frequency');
        title('Percentage of molecules in off-state');
        
        
end
end %plotData
