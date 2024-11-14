%% **********************************************************************************
%           H F p E F / H F r E F   Self Organizing Mapping A N A L Y S I S
% ***********************************************************************************
%
%   This script is based on Kayvan's suggestions. Using SOM unsupervised
%   neourn network maybe better to interpret the high dimentional data
%   which is totally independent.like model params.
%
%   Script created on:  04/22/2024 
%   Last modified on:   
%
%   Developed by        Feng Gu
%                       Physiological Systems Dyanmics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************


clear all

%% Load optimized parameter values
% Load optimized parameters
ParamsT = readtable("Features.xlsx",'VariableNamingRule','preserve');
% unsupervisedslot = [(1:8) (13:23) (26:28) (37:40) (45:58)];
% unsupervisedslot = (1:38);
% log transformation
% ParamsT{:,3:end} = log(ParamsT{:,3:end});

% unsupervisedslot = [(1:8) (13:23) (26:40) (45:50)];
% ParamsT = ParamsT(:,unsupervisedslot);
% check out std and delete the repeat ones
[~,locs] = unique(round(table2array(mean(ParamsT(:,3:end))),4));
locs = locs+2;
locs = sort(locs);
locs = [1;2;locs];
ParamsT = ParamsT(:,locs);
% replace the out of bound values by 95% confidence lower bound and upper
% bound 
data_normalized = zscore(ParamsT{:,3:end});
z_threshold = 3;
outliers1 = data_normalized > z_threshold;
outliers2 = data_normalized < -z_threshold;
[~,FeatureNumber] = size(data_normalized);
for i = 1:FeatureNumber
    ParamsT{:,i+2}(outliers1(:,i)) = max(ParamsT{:,i+2}(~outliers1(:,i)));
    ParamsT{:,i+2}(outliers2(:,i)) = max(ParamsT{:,i+2}(~outliers2(:,i)));
end
ANorm_Optp = zscore(ParamsT{:,3:end});

% checking distribution of the data
% for i = 1:30
%    figure(i);clf; histogram(ParamsT{:,i+2});
% end
%% checking covariance
A = cov(ANorm_Optp);
% Define the number of colors for each part (blue, white, red)
numColors = 100;
numBlue = round(numColors * 0.5);  % 50% of colormap for blue gradient
numWhite = numColors - 2 * numBlue;  % The rest is split between blue and red gradients -> 20% white
numRed = numBlue;                   % Equal amount of red and blue

% Create the blue to white gradient
blueToWhite = [linspace(0,1,numBlue)', linspace(0,1,numBlue)', linspace(0.85,1,numBlue)'];

% Create the white color in the middle
whiteMiddle = ones(numWhite, 3);  % Just white

% Create the white to red gradient
whiteToRed = [linspace(1,0.85,numRed)', linspace(1,0,numRed)', linspace(1,0,numRed)'];

% Combine them to form the new colormap
customColormap = [blueToWhite; whiteMiddle; whiteToRed];


% Use the new colormap with the heatmap
figure(11);clf;
hm = heatmap(A);
clim([-1 1]);
% Lable = {'C_S_A','C_S_V','C_P_A','C_P_V','G_S_A','G_P_A','G_t_S_A','G_t_P_A','kact_L_V','kact_R_V','kpas_L_V','kpas_R_V',...
%      'Vw_L_V','Vw_R_V','Vw_S_E_P','V_a_u','V_a_c','G_a','V4c','K_p_e_r_i','B_p_e_r_i','G_t_o','G_t_c','G_p_o','G_p_c','G_m_o','G_m_c','G_a_o','G_a_c','V_v_s'};
% hm.XDisplayLabels = Lable;
% hm.YDisplayLabels = Lable;
hm.Colormap = customColormap;
hm.CellLabelFormat = '%.2f';
clim(hm, [-1 1]); 

%% Use SOM to cluster patients

% if epcho unknow use comment code to calculate echo first
x = ANorm_Optp';  
optimalEpoch = 50000;

% Create a Self-Organizing Map
dimension1 = 10;
dimension2 = 10;

net = selforgmap([dimension1 dimension2]);
% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotsomtop','plotsomnc','plotsomnd', ...
    'plotsomplanes', 'plotsomhits', 'plotsompos'};
net.trainParam.epochs = optimalEpoch;

% Train the Network
[net,tr] = train(net,x);

%% DBscan to cluster the neourn based on distance which I built
% Get the similarity of SOM neurons, here I am using a weight+geometery
% distance as new distance martix
weights = net.IW{1,1};
[num_neurons, num_features] = size(weights);
distanceMatrix = zeros(num_neurons);

for i = 1:num_neurons
    for j = i+1:num_neurons
        % Calculate the distance between weights
        weight_distance = norm(weights(i, :) - weights(j, :));
        
        % Get the grid row and column position of the neuron
        [row_i, col_i] = ind2sub(net.layers{1}.dimensions, i);
        [row_j, col_j] = ind2sub(net.layers{1}.dimensions, j);
        
        % If it is an even-numbered row, increment the column coordinate by 0.5
        if mod(row_i, 2) == 0
            col_i = col_i + 0.5;
        end
        if mod(row_j, 2) == 0
            col_j = col_j + 0.5;
        end
        
        % Calculate the geometric distance in the hexagonal structure using the actual grid distances
        % The vertical distance takes into account the hexagon height
        hex_height = sqrt(3) / 2;
        grid_distance = sqrt((col_i - col_j)^2 + (hex_height * (row_i - row_j))^2);

        % Combined distance
        lambda = 1; % Controls the weight of the grid distance influence, adjust based on situation
        combined_distance = weight_distance + lambda * grid_distance;

        distanceMatrix(i, j) = combined_distance;
        distanceMatrix(j, i) = combined_distance;  % The distance matrix is symmetrical
    end
end
% distanceMatrix = weights;
epsilon = plotKDistance(distanceMatrix, 4);
% DBSCAN clustering
% epsilon = 9.5; % The neighborhood size (this value needs to be determined experimentally)
minpts = 4;   % The minimum number of points required to form a cluster (this value also needs to be experimentally adjusted)

% DBSCAN execution may require the Statistics Toolbox. If absent, you might need to download the appropriate implementation.
db_labels = dbscan(distanceMatrix, epsilon, minpts);

% The output of DBSCAN, db_labels, assigns a cluster index to each neuron. -1 denotes noise points.

% Reshape the labels into a 2D grid format corresponding to the SOM grid
gridLabels = reshape(db_labels, [dimension1, dimension2]);
gridLabels = gridLabels';
% Create a colormap, including a gray color and various other colors
numClusters = max(max(gridLabels)); % Number of clusters
visualizeSOMGrid(gridLabels, dimension1, dimension2, numClusters, 'DBSCAN Neuron Clusters',1);
%% Evaluate Optimal Number of Clusters for Kmeans and Hierarchical Clustering
% find a reasonable cluster number
rng("default") % Set the random number generator to default for reproducibility
maxClusters = 100;  % Maximum number of clusters to evaluate

% K-Means Clustering Evaluations
% Calculate Within-Cluster Sum of Square (elbow method) for Kmeans
wcssEvaluation = zeros(maxClusters, 1);
for k = 1:maxClusters
    [~, ~, sumd] = kmeans(distanceMatrix, k);
    wcssEvaluation(k) = sum(sumd);
end

% Calculate silhouette scores for Kmeans
silhouetteEvaluation = evalclusters(distanceMatrix, "kmeans", "silhouette", "KList", 1:maxClusters);
% Calculate Davies-Bouldin Index for Kmeans
dbiEvaluation = evalclusters(distanceMatrix, "kmeans", "DaviesBouldin", "KList", 1:maxClusters);
% Calculate Calinski-Harabasz Index for Kmeans
chiEvaluation = evalclusters(distanceMatrix, "kmeans", "CalinskiHarabasz", "KList", 1:maxClusters);

% Hierarchical Clustering Evaluations
% Calculate silhouette scores for Hierarchical Clustering
silhouetteHierEvaluation = evalclusters(distanceMatrix, "linkage", "silhouette", "KList", 1:maxClusters);
% Calculate Davies-Bouldin Index for Hierarchical Clustering
dbiHierEvaluation = evalclusters(distanceMatrix, "linkage", "DaviesBouldin", "KList", 1:maxClusters);
% Calculate Calinski-Harabasz Index for Hierarchical Clustering
chiHierEvaluation = evalclusters(distanceMatrix, "linkage", "CalinskiHarabasz", "KList", 1:maxClusters);

% Plot the evaluation metrics to assist with determining the optimal number of clusters
figure(2);

range = [0 20];
optimalNumClusters = 5;
% Plot the evaluation metrics to assist with determining the optimal number of clusters
figure(102);

% Plot Within-Cluster Sum of Square (elbow method) for Kmeans
subplot(2,4,1);
plot(1:maxClusters, wcssEvaluation, 'b-o');
title('Kmeans Elbow Method (WCSS)');
xlabel('Number of Clusters');
ylabel('Within-Cluster Sum of Square');
xlim(range);
hold on; 
scatter(optimalNumClusters, wcssEvaluation(optimalNumClusters), 'r', 'filled');
hold off
% Plot silhouette scores for Kmeans
subplot(2,4,2);
plot(1:maxClusters, silhouetteEvaluation.CriterionValues, 'b-o');
title('Kmeans Silhouette Coefficient');
xlabel('Number of Clusters');
ylabel('Silhouette Score');
xlim(range);
hold on; 
scatter(optimalNumClusters, silhouetteEvaluation.CriterionValues(optimalNumClusters), 'r', 'filled');
hold off
% Plot the Davies-Bouldin Index for Kmeans
subplot(2,4,3);
plot(1:maxClusters, dbiEvaluation.CriterionValues, 'b-o');
title('Kmeans Davies-Bouldin Index');
xlabel('Number of Clusters');
ylabel('DBI Score');
xlim(range);
hold on; 
scatter(optimalNumClusters, dbiEvaluation.CriterionValues(optimalNumClusters), 'r', 'filled');
hold off
% Plot the Calinski-Harabasz Index for Kmeans
subplot(2,4,4);
plot(1:maxClusters, chiEvaluation.CriterionValues, 'b-o');
title('Kmeans Calinski-Harabasz Index');
xlabel('Number of Clusters');
ylabel('CHI Score');
xlim(range);
hold on; 
scatter(optimalNumClusters, chiEvaluation.CriterionValues(optimalNumClusters), 'r', 'filled');
hold off
% Plot silhouette scores for Hierarchical Clustering
subplot(2,4,5);
plot(1:maxClusters, silhouetteHierEvaluation.CriterionValues, 'b-o');
title('Hierarchical Silhouette Coefficient');
xlabel('Number of Clusters');
ylabel('Silhouette Score');
xlim(range);
hold on; 
scatter(optimalNumClusters, silhouetteHierEvaluation.CriterionValues(optimalNumClusters), 'r', 'filled');
hold off
% Plot the Davies-Bouldin Index for Hierarchical Clustering
subplot(2,4,6);
plot(1:maxClusters, dbiHierEvaluation.CriterionValues, 'b-o');
title('Hierarchical Davies-Bouldin Index');
xlabel('Number of Clusters');
ylabel('DBI Score');
xlim(range);
hold on; 
scatter(optimalNumClusters, dbiHierEvaluation.CriterionValues(optimalNumClusters), 'r', 'filled');
hold off
% Plot the Calinski-Harabasz Index for Hierarchical Clustering
subplot(2,4,7);
plot(1:maxClusters, chiHierEvaluation.CriterionValues, 'b-o');
title('Hierarchical Calinski-Harabasz Index');
xlabel('Number of Clusters');
ylabel('CHI Score');
xlim(range);
hold on; 
scatter(optimalNumClusters, chiHierEvaluation.CriterionValues(optimalNumClusters), 'r', 'filled');
hold off

% Adjust the layout to prevent labels from overlapping
set(gcf, 'Position', [100, 100, 1024, 768]);  % Resize figure to make it wider
%% Perform K-means Clustering and Hierarchical Clustering on SOM Neurons and Visualize
% Perform K-means Clustering
% Choose the number of clusters (numClusters) for K-means (best determined with previous evaluation methods)
numClusters = 5;
rng('default');
% Perform K-means Clustering on the distance matrix
[kmeansClusters, ~] = kmeans(distanceMatrix, numClusters,'Replicates',100);

% Reshape the K-means cluster indices to match the SOM grid
kmeansClusterIndexMap = reshape(kmeansClusters, [dimension1, dimension2]);
kmeansClusterIndexMap = kmeansClusterIndexMap';
% Create a new mapping from original indices to new indices
kmeansIndexMapping = [1 2 3 4 5];  % For example, original index 2 maps to new index 3, and so on

% Replace original indices with new indices for K-means clustering
remappedKmeansClusterIndexMap = arrayfun(@(x) kmeansIndexMapping(x), kmeansClusterIndexMap);
% Visualize K-means clusters on the SOM grid
visualizeSOMGrid(remappedKmeansClusterIndexMap, dimension1, dimension2, numClusters, 'K-means SOM Neuron Clusters',3);

% Perform Hierarchical Clustering
% Cluster using 'ward' linkage method on the distance matrix
wardLinkage = linkage(distanceMatrix, 'ward');

% Generate clusters from the hierarchical clustering dendrogram
hierarchicalClusters = cluster(wardLinkage, 'Maxclust', numClusters);

% Reshape the hierarchical cluster indices to match the SOM grid
hierarchicalClusterIndexMap = reshape(hierarchicalClusters, [dimension1, dimension2]);
hierarchicalClusterIndexMap = hierarchicalClusterIndexMap';
% Create a new mapping from original indices to new indices
hierarchicalIndexMapping = [3 2 5 4 1];  % For example, original index 2 maps to new index 3, and so on

% Replace original indices with new indices for hierarchical clustering
remappedHierarchicalClusterIndexMap = arrayfun(@(x) hierarchicalIndexMapping(x), hierarchicalClusterIndexMap);

% Visualize hierarchical clusters on the SOM grid
visualizeSOMGrid(remappedHierarchicalClusterIndexMap, dimension1, dimension2, numClusters, 'Hierarchical SOM Neuron Clusters',4);

%% Visualize SOM based on label
% Predefined color values
Color_HFType = {
    [0.84, 0.08, 0.18],  % Red for HFrEF
    [0.93, 0.69, 0.13],  % Yellow for HFmrEF
    [0, 1, 0.725],        % Green
    [0, 0.725, 1],        % Cyan
    [0, 0.45, 0.74]       % Blue for HFpEF
};

% Determine the size of the color map
colorMapSize = 100;  % Total number of colors

% Determine the size of each gradient interval
halfMapSize = colorMapSize / 2;       % Red to Yellow occupies half the map

% Create gradients
gradientR2Y = zeros(halfMapSize, 3);  % Gradient from red to yellow
gradientY2G = zeros(17, 3);           % Gradient from yellow to green
gradientG2C = zeros(17, 3);           % Gradient from green to cyan
gradientC2B = zeros(16, 3);           % Gradient from cyan to blue

% Build red to yellow gradient
for i = 1:3
    gradientR2Y(:, i) = linspace(Color_HFType{1}(i), Color_HFType{2}(i), halfMapSize);
end

% Build yellow to green gradient
for i = 1:3
    gradientY2G(:, i) = linspace(Color_HFType{2}(i), Color_HFType{3}(i), 17);
end

% Build green to cyan gradient
for i = 1:3
    gradientG2C(:, i) = linspace(Color_HFType{3}(i), Color_HFType{4}(i), 17);
end

% Build cyan to blue gradient
for i = 1:3
    gradientC2B(:, i) = linspace(Color_HFType{4}(i), Color_HFType{5}(i), 16);
end

% Concatenate all gradients to form the final color map
colorMap = [gradientR2Y; gradientY2G; gradientG2C; gradientC2B];

% The colorMap now contains a gradient from red to yellow, then yellow to green,
% green to cyan, and cyan to blue.
% Assuming you have a SOM network 'net', and a dataset 'data'

% Apply the SOM network to the dataset, mapping data vectors to neurons in the SOM grid
mappings = vec2ind(net(x));  % Get the index of the most activated neuron for each data point

% Initialize a matrix for neuron colors and counters matching the size of the SOM grid
neuronCounts = zeros(dimension1, dimension2);  % Count matrix with dimensions of the SOM grid
neuronColors = zeros(dimension1, dimension2, 3);  % Color matrix with dimensions of the SOM grid
% Initialize matrix for neuron HF scores
neuronHFScore = zeros(dimension1, dimension2);

for i = 1:length(x)
    % Get the row and column indices of the neuron for each data point
    [col, row] = ind2sub(net.layers{1}.dimensions, mappings(i));
    % Increase counts
    neuronCounts(row, col) = neuronCounts(row, col) + 1;
    % Add score based on HF type (HFrEF = 0, HFmrEF = 0.5, HFpEF = 1)
    if strcmp('HFrEF', ParamsT.HFtype(i))
        neuronHFScore(row, col) = neuronHFScore(row, col) + 0;
    elseif strcmp('HFmrEF', ParamsT.HFtype(i))
        neuronHFScore(row, col) = neuronHFScore(row, col) + 0.5;
    elseif strcmp('HFpEF', ParamsT.HFtype(i))
        neuronHFScore(row, col) = neuronHFScore(row, col) + 1;
    end
end

% Calculate normalized scores and map to colormap indices
neuronColorIndices = round((neuronHFScore ./ max(neuronCounts, 1)) * (colorMapSize - 1)) + 1;

% Apply mapping to get the color for each neuron
for r = 1:dimension1
    for c = 1:dimension2
        if neuronCounts(r, c) > 0
            % Find the corresponding color using the computed index
            neuronColors(r, c, :) = colorMap(neuronColorIndices(r, c), :);
        end
    end
end

% Visualize the SOM grid and colors
figure(99); clf; hold on;
% Hexagon side length
hexRadius = 0.5;

% Calculate the hexagon's horizontal and vertical spacing based on hexagon side length and their relationship
% Since the hexagons are rotated, the width should be reduced, as the width is now point-to-point (not edge-to-edge)
hexWidth = hexRadius * sqrt(3);  % Distance from one edge of the hexagon to the opposite edge
horizSpacing = hexWidth;         % Horizontal spacing between centers of hexagons
vertSpacing = hexRadius * 1.5;   % Vertical spacing between hexagons, equal to 1.5 times the radius

for r = 1:dimension1
    for c = 1:dimension2
        % Calculate the current hexagon’s center position
        if mod(r, 2) == 0
            % Even rows are offset horizontally to center the hexagons in each row
            centerPos = [(c + 0.5) * horizSpacing, r * vertSpacing];
        else
            % Odd rows are not offset
            centerPos = [c * horizSpacing, r * vertSpacing];
        end
        
        % Calculate the coordinates of the six vertices of the hexagon
        % Start angle is -π/6 to have a point upwards
        angles = (-pi/6) + (0:5) * (2*pi/6);
        hexagonX = centerPos(1) + hexRadius * cos(angles);
        hexagonY = centerPos(2) + hexRadius * sin(angles);
        
        % Draw and fill hexagon with color
        color = squeeze(neuronColors(r, c, :))'; % Extract color
        patch(hexagonX, hexagonY, color, 'EdgeColor', 'white', 'LineWidth', 0.01);
        
        % Add text for counts at the hexagon center
        text(centerPos(1), centerPos(2), num2str(neuronCounts(r, c)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20);
    end
end
% Add a color bar after drawing is complete
colormap(colorMap);  % Specify the colormap
c = colorbar;
c.Limits = [0 1];  % Set the range of the color bar to match the data range

% If needed, set the ticks and tick labels for the color bar
% Here we assume the color bar is divided into 'HFrEF', 'HFmrEF', 'HFpEF'
c.Ticks = [0, 0.5, 1];  % Set the locations of ticks on the color bar
c.TickLabels = {'HFrEF', 'HFmrEF', 'HFpEF'};  % Set tick labels
c.FontSize = 20;
axis equal;
axis off;
%% Find patient locations in the SOM mapping
[colPositions, rowPositions] = patientPositionsInSOM(net, x);
KmeansClusterAssignments = NaN(length(colPositions), 1);
HierarchicalClusterAssignments = NaN(length(colPositions), 1);

% Populate the cluster assignments from each clustering algorithm for each patient
for i = 1:length(colPositions)
    KmeansClusterAssignments(i) = remappedKmeansClusterIndexMap(rowPositions(i), colPositions(i));
    HierarchicalClusterAssignments(i) = remappedHierarchicalClusterIndexMap(rowPositions(i), colPositions(i));
end

% Combine the cluster assignments into one summary matrix
SumClusters = [KmeansClusterAssignments HierarchicalClusterAssignments];

% Initialize the result vectors
rowsEqual = false(size(SumClusters, 1), 1);            % Rows where all elements are equal
rowsDifferent = false(size(SumClusters, 1), 1);        % Rows where all elements are different

% Check the relationships between elements in each row
for i = 1:size(SumClusters, 1)
    if SumClusters(i, 1) == SumClusters(i, 2)
        rowsEqual(i) = true;
    else
        rowsDifferent(i) = true;
    end
end

% Get indices of rows for each condition
rowsWithAllEqual = find(rowsEqual);            % Indices of rows with all elements equal
rowsWithAllDifferent = find(rowsDifferent);    % Indices of rows with all elements different

% Calculate the number of occurrences for each condition
numRowsAllEqual = sum(rowsEqual);              % Number of rows where all elements are equal
numRowsAllDifferent = sum(rowsDifferent);      % Number of rows where all elements are different

% Display the results
fprintf('Number of rows with all elements equal: %d\n', numRowsAllEqual);
fprintf('Indices of rows: %s\n', mat2str(rowsWithAllEqual'));

fprintf('Number of rows with all elements different: %d\n', numRowsAllDifferent);
fprintf('Indices of rows: %s\n', mat2str(rowsWithAllDifferent'));
%% Initialize risk level column array
% Initialize the risk level array
riskLevel = zeros(size(SumClusters, 1), 1);

% Iterate over each row in SumClusters
for i = 1:size(SumClusters, 1)
    % Get the current row values
    values = SumClusters(i, :);
    
    % Check if both values are low risk
    if all(values == 1 | values == 3)
        riskLevel(i) = 1; % Low risk
    % Check if both values are high risk
    elseif all(values == 4 | values == 5)
        riskLevel(i) = 3; % High risk
    % Check if both values are medium risk
    elseif all(values == 2)
        riskLevel(i) = 2; % Medium risk
    % Check if there is one medium risk and one high/low risk
    else
        % Medium and high risk
        if any(values == 2) && any(values == 4 | values == 5)
            riskLevel(i) = 3; % High risk
        % Medium and low risk
        elseif any(values == 2) && any(values == 1 | values == 3)
            riskLevel(i) = 1; % Low risk
        % High risk and low risk
        elseif any(values == 1 | values == 3) && any(values == 4 | values == 5)
            riskLevel(i) = 2; % Medium risk
        end
    end
end
% Add the risk level array to SumClusters matrix to create a new table
SumClustersWithRisk = [SumClusters, riskLevel];

%% Find the portion of each cluster
% Replace -1 in DBSCAN clustering results with 5
DBSCANClusterAssignments(DBSCANClusterAssignments == -1) = 6;

% Obtain all possible diagnoses types
uniqueDiagnoses = unique(ParamsT.HFtype);

% Initialize tables with rows representing different cluster numbers and columns representing different diagnoses
numClusters = max([DBSCANClusterAssignments; KmeansClusterAssignments; HierarchicalClusterAssignments]);
diagnosisTableDBSCAN = array2table(zeros(numClusters, numel(uniqueDiagnoses)), 'VariableNames', matlab.lang.makeValidName(uniqueDiagnoses));
diagnosisTableKmeans = array2table(zeros(numClusters, numel(uniqueDiagnoses)), 'VariableNames', matlab.lang.makeValidName(uniqueDiagnoses));
diagnosisTableHierarchical = array2table(zeros(numClusters, numel(uniqueDiagnoses)), 'VariableNames', matlab.lang.makeValidName(uniqueDiagnoses));

% Calculate the composition of diagnoses for each clustering method
for i = 1:length(DBSCANClusterAssignments)
    diagnosis = ParamsT.HFtype(i);  % Current patient's diagnosis
    diagnosisIndex = find(strcmp(uniqueDiagnoses, diagnosis));  % Find the column index for the diagnosis

    % Update the count of diagnoses for the respective cluster number
    dbscanClusterIdx = DBSCANClusterAssignments(i);
    kmeansClusterIdx = KmeansClusterAssignments(i);
    hierarchicalClusterIdx = HierarchicalClusterAssignments(i);

    diagnosisTableDBSCAN{dbscanClusterIdx, diagnosisIndex} = diagnosisTableDBSCAN{dbscanClusterIdx, diagnosisIndex} + 1;
    diagnosisTableKmeans{kmeansClusterIdx, diagnosisIndex} = diagnosisTableKmeans{kmeansClusterIdx, diagnosisIndex} + 1;
    diagnosisTableHierarchical{hierarchicalClusterIdx, diagnosisIndex} = diagnosisTableHierarchical{hierarchicalClusterIdx, diagnosisIndex} + 1;
end

% Display tables for review
disp('DBSCAN Cluster Diagnoses:');
disp(diagnosisTableDBSCAN);
disp('Kmeans Cluster Diagnoses:');
disp(diagnosisTableKmeans);
disp('Hierarchical Cluster Diagnoses:');
disp(diagnosisTableHierarchical);