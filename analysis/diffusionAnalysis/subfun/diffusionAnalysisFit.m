function [D_final] = diffusionAnalysisFit(msd_result,experimParam,fitParam,printParam)
% [D_final] = diffusionAnalysisFit(msd_result,experimParam,fitParam,printParam) 
% Obtain diffusion coefficient from MSD curve obtained from MSD histograms.
% See diffusionAnalysisHistogram.m.
% 
% INPUT:
%     msd_result: [OxD] cell array where O is the order and D is the
%     dimension. Every cell contains a double array
%     [t,sigma_i,error_sigma_i,weight_i,error_weight_i,...,actual].
%     Here, t is the time in seconds, sigma is the mean squared
%     displacement in [µm^2] and weight is the weight of the diffusion
%     species. See diffusionAnalysisHistogram.m for details.
%     
%     experimParam: Struct of experimental parameters. See
%     RunDiffusionAnalysis.m for details.
%     
%     fitParam: Struct of fitting parameters. See RunDiffusionAnalysis.m
%     for details.
%     
%     print: Struct of print/plotting parameters. See
%     RunDiffusionAnalysis.m for details.
%     
%     
% OUTPUT:
%     D_final: Struct line array, one line per diffusing species order. It
%     contains a double array for D, the MSD y-intersection and the optimal
%     number of fitting parameters. The first line contains the values, the
%     second line the errors.


order = size(msd_result,1);
D_final = repmat(struct('D_x',[],'MSD_t0_x',[],'N_opt_x',[],'D_y',[],'MSD_t0_y',[],'N_opt_y',[]),order,1);

%normalization to 20°C
if ~isempty(experimParam.temperature)
    T = 273.15+experimParam.temperature;
    T_norm = 293.15;
    visc_table = load(experimParam.viscosityTable); visc_table = visc_table.viscosity;
    eta = interp1(visc_table(:,1),visc_table(:,2),T,'spline');
    eta_norm = interp1(visc_table(:,1),visc_table(:,2),T_norm,'spline');
    norm_factor = T_norm/T*eta/eta_norm;
end


for dim=1:size(msd_result,2)
    for iOrder=1:order
        msd_curve = msd_result{iOrder,dim};
        min_data_points = max(fitParam.minimumSkip,3); %need at least 3 points to fit
        
        if fitParam.fitTrueMSDValues
            msd_curve = msd_curve(msd_curve(:,end)==1,:); %discard orders which did not lead to best fit result
        end
        max_data_points = size(msd_curve,1);
        
        if max_data_points>=min_data_points
            D_error = zeros(max_data_points-min_data_points+1,2);
            
            %try all MSD curves
            for jNopt = min_data_points:max_data_points
                [D,msd_t0,~] = fitMSD(msd_curve(1:jNopt,:),fitParam,[]); %printParam is empty, we don't plot yet

                D_error(jNopt-min_data_points+1,:) =[mean(abs(D(:,2)./D(:,1))),mean(abs(msd_t0(:,1)))]; %error score: relative error of D and y-axis intersection
            end
            
            score = size(D_error,1):-1:1;
            [~,idx] = sort(D_error(:,1)); D_error(:,1) = score(idx);
            [~,idx] = sort(D_error(:,2)); D_error(:,2) = score(idx);
            
            [~,N_opt] = max(sum(D_error,2)); N_opt = N_opt-1+min_data_points; %optimal fit is the one with the lowes possible error for D and the loest y-intersection (as this is also a measure for the error)
        else
            N_opt = max_data_points;
        end
        
        %final fit
        if N_opt>=min_data_points
            [D,msd_t0,~] = fitMSD(msd_curve(1:N_opt,:),fitParam,printParam);
        else
            D = [];
            msd_t0 = [];
        end
        
        %renormalize diffusion coefficient to 20°C using provided viscosity table
        if ~isempty(experimParam.temperature)
            D = D*norm_factor;
            msd_t0 = msd_t0*norm_factor;
        end
        
        %save values
        if dim==1
            D_final(iOrder).D_x = D;
            D_final(iOrder).MSD_t0_x = msd_t0;
            D_final(iOrder).N_opt_x = N_opt;
        else
            D_final(iOrder).D_y = D;
            D_final(iOrder).MSD_t0_y = msd_t0;
            D_final(iOrder).N_opt_y = N_opt;
        end
    end %END iOrder
end %end dim

