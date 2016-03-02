function [D,msd_t0,v] = fitMSD(msd_curve,fitParam,printParam)
% [D,msd_t0,v] = fitMSD(msd_curve,fitParam,printParam)
% Fit provided MSD curve and return diffusion coefficient D, y-intersection msd_t0, and velocity v.
% 
% INPUT:
%     msd_curve: Double array
%     [t,sigma_i,error_sigma_i,weight_i,error_weight_i,...,actual]. See
%     diffusionAnalysis(MSD)Fit.m for details.
%         
%     fitParam: Struct array of fitting parameters, see
%     RunDiffusionAnalysis.m for details.
%     
%     printParam: Struct array of print/plotting parameters, see
%     RunDiffusionAnalysis.m for details.
%     
%     
% OUTPUT:
%     D: Double array for diffusion coefficient, one line per order. First
%     column contains value, second column contains error. Units are in
%     native units of msd_curve.
%     
%     msd_t0: Double array of y-intersection values.
%     
%     v: Double array of velocities.

order = (size(msd_curve,2)-2)/4;
x = msd_curve(:,1);
fit_weights = fitParam.weightedFit;
v = [];
if strcmp(fitParam.method,'msd')
    use_v = fitParam.fitVelocity;
else
    use_v = false;
end

D = zeros(order,2);
msd_t0 = zeros(order,2);
R = zeros(order,1);
h_plot = zeros(order,1);

first_plot = true;
markervec = {'o','*','.','^','v','>','<','x','+','s','d','p','h',};
colorvec = distinguishable_colors(order,'w');

for iOrder=1:order
    y = msd_curve(:,2+(iOrder-1)*4);
    if fit_weights
        y_weight = 1./msd_curve(:,3+(iOrder-1)*4);
    else
        y_weight = ones(size(y,1),1);
    end
    
    if use_v
        [p,~] = polyfit3(x,y,2,[],y_weight);
    else
        [p,~] = polyfit3(x,y,1,[],y_weight);
    end
    
    res = polyval(p,x)-y;
    if use_v
        J = [x.^2,x,ones(size(x,1),1)];
    else
        J = [x,ones(size(x,1),1)];
    end
    dof = size(x,1)-numel(p)-1;
    chi2 = sum(res.^2)/dof;
    pVarCovMat = chi2*inv(J'*J); %#ok<MINV>
    dp = sqrt(diag(pVarCovMat));
    
    if use_v && p(1)>0
        v = [sqrt(p(1)),dp(1)/(2*p(1))];
        D(iOrder,:) = [p(2)/2,dp(2)/2];
        msd_t0(iOrder,:) = [p(3),dp(3)];
    else 
        D(iOrder,:) = [p(1)/2,dp(1)/2];
        msd_t0(iOrder,:) = [p(2),dp(2)];
    end
    if strcmp(fitParam.method,'msd') && strcmp(fitParam.distanceMethod,'jump')
        D = D/2; %in this case, MSD slope would be 4*D (from x^2+y^2 instead of just x^2) instead of 2*D, correct this now.
    end  
    R(iOrder) = 1-chi2/var(y);
      
    if ~isempty(printParam)
        if first_plot
            h = figure;
            ymax = max(y);
        else
            ymax = max([ymax;y]);
        end
        
        if ~fit_weights
            plot(x*1e3,y,markervec{iOrder},'color',colorvec(iOrder,:),'MarkerFaceColor',colorvec(iOrder,:)); %t in ms
        else
            errorbar(x*1e3,y,1./y_weight,markervec{iOrder},'color',colorvec(iOrder,:),'MarkerFaceColor',colorvec(iOrder,:));
        end
        
        
        if first_plot
            hold on;
            first_plot = false;
        end
        
        h_plot(iOrder) = plot([0;x*1e3],[p(end);res+y],'-','color',colorvec(iOrder,:));
        if iOrder==order
            hold off;
        end
    end %END plot fit
end %end iOrder


if ~isempty(printParam)
    ylim([0,max(y)*1.5]);
    legend_strings = cell(0,0);
    
    for iOrder=1:order
        legend_strings = [legend_strings,{['Order ',int2str(iOrder)]}];
        if printParam.latex
            text(0.1,0.7-(iOrder-1)*0.11,['$R_',int2str(iOrder),' = ' num2str(R(iOrder),3),'$'],'Units','normalized');
        else
            text(0.15,0.64-(iOrder-1)*0.11,['R_',int2str(iOrder),' = ' num2str(R(iOrder),3)],'Units','normalized');
            text(0.1,0.7-(iOrder-1)*0.11,['D_',int2str(iOrder),' = ' num2str(D(iOrder,1),3),' \pm ',num2str(D(iOrder,2),2),' µm^2/s'],'Units','normalized');
        end
    end
    
    legend(h_plot,legend_strings,'Location','NorthEast');
    box on;
    
    if printParam.latex
        xlabel('$\Delta t$ [ms]','Interpreter','none');
        ylabel('$\sigma^2$ [$\micro\metre^2$]','Interpreter','none');
        set(h,'Position',[0 0 40*printParam.size(1:2)],'Units','centimeters');
        matlabfrag([printParam.filename,'-msd',num2str(fitplot(1)),num2str(fitplot(2))],'epspad',[0 0 0 0]);
    else
        xlabel('\Delta t [ms]','fontsize',12);
        ylabel('\sigma^2 [µm^2]','fontsize',12);
    end
end

