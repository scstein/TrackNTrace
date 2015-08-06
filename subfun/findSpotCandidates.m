function candidatePos = findSpotCandidates(avg, sigma, CORR_THRESH, plot_enable)

%% Candidate selection by normalized cross correlation
model = gaussian2d(round(10*sigma),sigma,1);
pattRows = size(model,1);
pattCols = size(model,2);

MIN_DIST = round(6*sigma);
[ match_data, match_img ] = matchSpot( avg, model, CORR_THRESH, MIN_DIST);

%% Post processing selection
% amplitudes = avg(logical(match_img));
% [amplitudes, idx] = sort(amplitudes);
% match_data = match_data(idx,:);
%
% fitPos = (round(size(amplitudes)/2)-round(size(amplitudes)*0.25) : round(size(amplitudes)/2)+round(size(amplitudes)*0.25))';
% ampli_fit = amplitudes(fitPos);
%
% % Fit linear function to center part of sorted amplitudes
% p = polyfit(fitPos,ampli_fit,1);
% yPoly = polyval(p,(1:numel(amplitudes))');
%
% % Kick out outliers
% margin = 0.025;
% residuals = abs(amplitudes-yPoly);
% inliers = residuals < margin*yPoly;
%
% % % Plotting
% % xInl = find(inliers);
% % ampInl = amplitudes(inliers);
% % figure;
% % plot( 1:numel(amplitudes),amplitudes,'r.', 1:numel(amplitudes), yPoly,'b', 1:numel(amplitudes), yPoly*(1+margin),'b--',1:numel(amplitudes), yPoly*(1-margin),'b--', xInl,ampInl,'g');
%
%
% match_data = match_data(inliers,:);
% amplitudes = amplitudes(inliers,:);

%% convert to pixel position
match_data(:,1) = match_data(:,1)+ceil(pattRows/2);
match_data(:,2) = match_data(:,2)+ceil(pattCols/2);

% ATTENTION! Here we switch the coordinates
% candidatePos(:,1) is columns (x-coordinate)
% candidatePos(:,2) is rows    (y-coordinate)
candidatePos = match_data(:,[2,1]);

%% plot
if plot_enable
    h=figure;
    mim(avg);
    hold on
    plot(match_data(:,2),match_data(:,1),'go');
    hold off
    title(gca,'Normalized cross correlation');
end

end


