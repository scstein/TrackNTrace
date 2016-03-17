function totm1 = CalculateDark(dark)
% Dark correction image, see Hirsch et. al. -  A Stochastic Model for Electron Multiplication Charge-Coupled Devices – From Theory to Practice 
	
totm1 = mean(dark,3);
rowm = repmat(mean(totm1,2),1,size(dark,2));
colm = repmat(mean(totm1,1),size(dark,1),1);
totm1 = repmat(mean(totm1(:)),[size(dark,1),size(dark,2)]);
totm1 = totm1-colm-rowm;

end