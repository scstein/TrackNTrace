function [traj_D,traj_v] = fitMSDSingle(traj,use_v,use_iso,use_displacement,fitParam)

dim = 1+~use_iso;
traj_D = zeros(1,2*dim);
traj_v = zeros(1,2*dim);

msd_length = min(size(traj,1)-1,fitParam.maximumSkip);
if msd_length<(3+use_v)
    return; %in this case, we cannot fit an msd curve, so abort immediately
end

if use_iso
    msd = zeros(msd_length,2);
else
    msd = zeros(msd_length,4);
end


for iFrame = 1:msd_length
    if fitParam.correlated
        displace = (traj(1+iFrame:iFrame:end,:)-traj(1:iFrame:end-iFrame,:)).^2;
    else
        displace = (traj(1+iFrame:1:end,:)-traj(1:1:end-iFrame,:)).^2;
    end
    
    if use_iso
        displace = displace(:);
    end
    
    if use_displacement
        for jDim=1:dim
            msd(iFrame,1+(jDim-1)*2:jDim*2) = [mean(displace(:,jDim)),std(displace(:,jDim))];
        end
    else
        msd(iFrame,:) = [mean(sum(displace.^2,2)),std(sum(displace.^2,2))];
    end
end

t = 1:msd_length;

for jDim=1:dim
    [D,v] = MSDFit(t(:),msd(:,1+(jDim-1)*2:jDim*2),use_v,fitParam.weightedFit,use_displacement);
    traj_D(1+(jDim-1)*2:jDim*2) = D;
    if use_v
        if ~isempty(v)
            traj_v(1+(jDim-1)*2:jDim*2) = v;
        else
            traj_v(1+(jDim-1)*2:jDim*2) = [Inf,Inf];
        end
    end
end
end %MAINFUN


function [D,v] = MSDFit(x,msd_curve,use_v,fit_weights,use_displacement)
if fit_weights
    y_weight = 1./msd_curve(:,2);
else
    y_weight = ones(size(msd_curve,1),1);
end

if use_v
    [p,~] = polyfit3(x,msd_curve(:,1),2,[],y_weight);
else
    [p,~] = polyfit3(x,msd_curve(:,1),1,[],y_weight);
end

res = polyval(p,x)-msd_curve(:,1);
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
    D = [p(2)/2,dp(2)/2];
else
    v = [];
    D = [p(1)/2,dp(1)/2];
end
if ~use_displacement
    D = D/2; %in this case, MSD slope would be 4*D (from x^2+y^2 instead of just x^2) instead of 2*D, correct this now.
end

end %MSDFIT


