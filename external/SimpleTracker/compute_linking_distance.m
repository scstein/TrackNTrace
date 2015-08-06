function max_linking_distance = compute_linking_distance(F1, F2)
% Added by Simon Christoph Stein
% We estimate the linking distance by producing a histogram over all
% distances from points in frame F1 to those in frame F2. The maximum in
% this distribution is taken as an estimate of the linking distance.


% Compute a^2 and b^2
F1_NormVec = sum(F1.^2,2);
F2_NormVec = sum(F2.^2,2);

F1NormMat = repmat(F1_NormVec, 1, size(F2,1) );
F2NormMat = repmat(F2_NormVec',  size(F1,1), 1);

% Compute MixTerm
MixTermMat = F1*F2';

% Overall euclidean distance squared
Dist2Mat = F1NormMat + F2NormMat - 2*MixTermMat;
DistMat = sqrt(Dist2Mat);

% Bin the distances, find the maximum
[n, xout] = hist(DistMat(:),1000); % This may be bad, as resolution depends linearly on the sensor pixel size..
[m,i] = max(n);
max_linking_distance = xout(  min(i+1,size(xout,2) )  );

end