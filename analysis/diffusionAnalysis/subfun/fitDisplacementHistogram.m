function [p_total,order,plotdata] = fitDisplacementHistogram(x,y,p0,use_v,order,fit_curve,alpha,plot_fit,frame,px,dt)

% lsqnonlin fit parameters
options = optimset('Jacobian','off',...
    'MaxFunEvals', 1e6, ...
    'MaxIter', 1e6, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

% Loop variables
first_iteration = true;
keep_fitting = true;

% Fit function handle
if strcmp(fit_curve,'jump')
    fitfun = @fitMultiRGaussian;
else
    fitfun = @fitMultiGaussian;
end

if use_v %careful, v has dimension [px] at first_iteration
    lb = [zeros(order-1,1);zeros(order,1);-1000]; %bounds: [weights;sigmas;velocity]
    ub = [ones(order-1,1);inf(order,1);1000];
else
    lb = [zeros(order-1,1);zeros(order,1)]; %bounds: [weights; sigmas]
    ub = [ones(order-1,1);inf(order,1)];
end


% Main fit loop
p_total = cell(100,1); %initialize very large array
while keep_fitting
    if ~first_iteration
        % store old parameters
        dof = dof_new;
        chi2 = chi2_new;
        p = p_new;
        dp = dp_new;
        residuals = residuals_new;
        
        % ... and prepare new ones
        [p0,lb,ub,order] = prepareParameters(p,use_v,order);
    end
    
    [p_new,~,residuals_new,~,~,~,J] = lsqnonlin(fitfun,p0,lb,ub,options,x,y,use_v);
    
    dof_new = size(x,1)-size(p_new,1)-1; % degrees of freedom
    J = full(J); % Jacobian, transform from sparse to full
    chi2_new = sum(residuals_new.^2)/dof_new; % reduced chi² for unknown weights
    pVarCovMat = chi2_new*inv(J'*J); %#ok<MINV> % variance-covariance matrix of parameters
    dp_new = sqrt(diag(pVarCovMat)); % parameter uncertainties, vector
    
    p_total(order) = {[p_new(:).';dp_new(:).']};
    
    if first_iteration 
        first_iteration = false;
    else
        % compare against previous residuals with an F-test
        %1-sided F-test: H0=: F = 1 (null hypothesis), H1: F < 1
        paramF = chi2_new/chi2; % new/old, will always be <1
        pValue = fcdf(paramF,dof,dof_new);
        
        if pValue < alpha
            keep_fitting = true;
        else
            keep_fitting = false;
            order = order-1;
            p_total = p_total(1:order);
        end
    end
end %END while keep_fitting


% Return final parameters and plot
plotdata = [x,y,y-residuals];
if plot_fit
    weight = [1,0];
    if order>1
        weight = [[p(1:order-1);1-sum(abs(p(1:order-1)))],[dp(1:order-1);sqrt(sum(dp(1:order-1).^2))]];
    end
    sigma_sq = p(order:2*order-1);
    velocity = 0;
    if use_v
        velocity = p(end);
    end
        
    figure;
    subplot(3,1,1:2);
    bar(x,y);
    
    hold on;
    plot(x,plotdata(:,3),'r','LineWidth',3);
    for i=1:order
        if strcmp(fit_curve,'jump')
            plot(x,weight(i,1)/sigma_sq(i)*x.*exp(-(x).^2/(2*sigma_sq(i))),'g','LineWidth',2);
        else
            plot(x,weight(i,1)/sqrt(2*pi*sigma_sq(i))*exp(-(x-velocity).^2/(2*sigma_sq(i))),'g','LineWidth',2);
        end
    end
    
    xlabel('Displacement [px]','fontsize',12);
    ylabel('Frequency [-]','fontsize',12);
    title(['Result of \chi^2 fit for \Delta_t = ',int2str(frame),' frames.'],'fontsize',12);
    
    text(0.7,0.75,['\chi^2_{\nu-1} = ' num2str(chi2)],'Units','normalized');
    for i=1:order
        text(0.7,0.7-i*0.1,['D_',int2str(i),'= ',num2str(sigma_sq(i)/(2*frame)*px^2/dt,3) '\pm ', num2str(dp(order-1+i)/(2*frame)*px^2/dt,3)],'Units','normalized');
        text(0.75,0.65-i*0.1,['w_',int2str(i),'= ',num2str(weight(i,1),3) '\pm ', num2str(weight(i,2),3)],'Units','normalized');
    end
    
    subplot(3,1,3);
    plot(x,zeros(size(x)),'r-');
    hold on;
    plot(x,residuals,'b.');
    title({'Residuals'},'fontsize',12);
end %plot_fit end

end %FUNCTION END


function [p,lb,ub,order_new] = prepareParameters(p0,use_v,order)
%insert new weight and diffusion coefficient

velocity = [];
weight = p0(1:order-1);
sigma_sq = p0(order:2*order-1);
if use_v
    velocity = p0(end);
end

if order == 1
    p = [0.5;sigma_sq/2;sigma_sq*2];
else
    p = [[weight-1/(2+2*order);1/(order+1)];[sigma_sq;1.5*max(sigma_sq)]];
end

p = [p;velocity];
order_new = order+1;

if use_v
    lb = [zeros(order_new-1,1);zeros(order_new,1);-1000];
    ub = [ones(order_new-1,1);inf(order_new,1);1000];
else
    lb = [zeros(order_new-1,1);zeros(order_new,1)];
    ub = [ones(order_new-1,1);inf(order_new,1)];
end

end

