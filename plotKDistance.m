function epsilon = plotKDistance(data, k)
% PLOTKDISTANCE Fits a polynomial curve to k-distance plot and identifies the
% point where second derivative falls below threshold to find epsilon for DBSCAN.

% Compute the Euclidean distance matrix between all pairs of data points
distMatrix = pdist2(data, data, 'euclidean');

% Find the distance to the k-th nearest neighbor for each data point
kDist = zeros(size(data, 1), 1);
for i = 1:size(data, 1)
    sortedDists = sort(distMatrix(i,:), 'ascend');
    kDist(i) = sortedDists(k + 1);  % Exclude the point itself (distance is 0)
end

% Sort the k-distances in descending order for fitting and plotting
[sortedKDist, ~] = sort(kDist, 'descend');
indices = 1:length(sortedKDist);

% Fit a 5th degree polynomial to the k-distance plot
[fitResult, ~] = fit(indices', sortedKDist, 'poly5');

% Calculate the first derivative of the fitResult
deriv1 = differentiate(fitResult, indices);

% Calculate the second derivative as the derivative of first derivatives
deriv2 = diff(deriv1);

% Plotting the original k-distance plot and polynomial fit
figure;
subplot(3, 1, 1);
plot(indices, sortedKDist, 'b.-');
hold on;
plot(indices, feval(fitResult, indices), 'r-', 'LineWidth', 2);
title('k-Distance Plot with Polynomial Fit');
xlabel('Point Index');
ylabel('k-distance');
legend('k-distance', 'Polynomial Fit');
grid on;

% Plotting the first derivative of the polynomial fit
subplot(3, 1, 2);
plot(indices, deriv1, 'k.-');
title('First Derivative (Slope)');
xlabel('Point Index');
ylabel('First Derivative');
grid on;

% Plotting the second derivative of the polynomial fit
subplot(3, 1, 3);
plot(indices(1:end-1), deriv2, 'r.-'); % Reduced index size by one due to diff
title('Second Derivative (Concavity Change)');
xlabel('Point Index');
ylabel('Second Derivative');
grid on;

% Find the first index where the second derivative is below the threshold
threshold = -2e-4;
belowThresholdIdx = find(deriv2 > threshold, 1, 'last')+1;
if isempty(belowThresholdIdx)
    epsilon = NaN; % In case no such point is found, return NaN
else
    candidateIdx = belowThresholdIdx + 1; % Reduced index size by one due to diff
    epsilon = mean([feval(fitResult, candidateIdx)...% Use the smoothed value from the fit
        sortedKDist(candidateIdx)]); % Use rawdata
    
    % Plot a marker on the smoothed curve where the selected epsilon is located
    subplot(3, 1, 1);
    hold on;
    plot(candidateIdx, epsilon, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
    
    % Also plot a marker on the second derivative plot
    subplot(3, 1, 3);
    hold on;
    plot(candidateIdx, deriv2(belowThresholdIdx), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
end

hold off;
end