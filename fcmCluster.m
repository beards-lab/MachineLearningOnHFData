function clusterIdx = fcmCluster(X, K)
    % Perform Fuzzy C-Means clustering
    [~, U, ~] = fcm(X, K);
    % Assign each data point to the cluster with the highest membership value
    [~, clusterIdx] = max(U,[],1); % Ensure clusterIdx is a column vector
    clusterIdx = clusterIdx';
    uniqueClusters = unique(clusterIdx);
    for i = 1:length(uniqueClusters)
        clusterIdx(clusterIdx == uniqueClusters(i)) = i;
    end
end