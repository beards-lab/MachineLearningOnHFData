%% Script Summary:
% This script applies knnsearch to cluster patients who are not consistently clustered in script ClusterOnDT.m.
% The input is the "Features.xlsx" file, containing 25 functional parameters and 4 simulated regurgitation fractions, and clustering results from ClusterOnDT.m.
% Outputs include clustering labels

% Created by Feng Gu
% Last modified: 10/29/2024

clear
Color_Clust = {[0.4940 0.1840 0.5560], [0 0.75 0.75], [0.8500 0.3250 0.0980], [0.5 0.5 0.5], [0.4660 0.6740 0.1880]};
Color_Clust{end+1} = [1.0 0.84 0.0]; 
Color_Clust{end+1} = [0.0 1.0 0.0];  

%% Load optimized parameter values and Preprocessing

ParamsT = readtable("Features.xlsx",'VariableNamingRule','preserve',Sheet='Sheet1');
ParamsT.MVr_S = (ParamsT.MVr_S-1).*20;
ParamsT.AVr_S = (ParamsT.AVr_S-1).*20;
ParamsT.PVr_S = (ParamsT.PVr_S-1).*15;
ParamsT.TVr_S = (ParamsT.TVr_S-1).*17.5;

[~,locs] = unique(round(table2array(mean(ParamsT(:,3:end))),4));
locs = locs+2;
locs = sort(locs);
locs = [1;2;locs];
ParamsT = ParamsT(:,locs);

data_normalized = zscore(ParamsT{:,3:end});
z_threshold = 3;
outliers1 = data_normalized > z_threshold;
outliers2 = data_normalized < -z_threshold;
[~,FeatureNumber] = size(data_normalized);
for i = 1:FeatureNumber
    ParamsT{:,i+2}(outliers1(:,i)) = max(ParamsT{:,i+2}(~outliers1(:,i)));
    ParamsT{:,i+2}(outliers2(:,i)) = min(ParamsT{:,i+2}(~outliers2(:,i)));
end

UCCParamsT  = ParamsT(ParamsT.Risk==4,:); % Risk = 4 consistently denotes patients who are inconsistently clustered.
UCClusterResults = UCCParamsT.Risk;
ParamsT  = ParamsT(~(ParamsT.Risk==4),:); 
ClusterResults = ParamsT.Risk;
UCCParamsT  = UCCParamsT(:,1:width(UCCParamsT)-1); % patients who are inconsistently clustered
ParamsT = ParamsT(:,1:width(ParamsT)-1); % patients who are consistently clustered

supervisedslot = [(4:26) 28 29]; % These differences in features, identified through ANOVA multiple comparisons, are significant among the consistently clustered patients.
supervised = 1;
if supervised == 1
    slot = supervisedslot;
else
    slot = (3:width(ParamsT));
end

% Prepare for mapping UCC patients
Optp_Names = ParamsT.Properties.VariableNames(slot) ;
A_Optp = ParamsT{:,slot};
UCC_Optp = UCCParamsT{:,slot};
mu = mean(A_Optp);
sigma = std(A_Optp);
ANorm_Optp = zscore(A_Optp);
%% Run the PCA on optimized parameters only
% Running the SVD with rows as each patient and columns as each optim param
[U_Optp,S_Optp,V_Optp] = svd(ANorm_Optp);       % Singular value decomposition
% In this case since A = U*S*V', U will represent the rotation in the patient
%  dimension and V' will represent the rotation in the optim param dimension.
%  Therefore our PCA score will be U*S and our PCA loadings are columns of V
AV_Optp = U_Optp * S_Optp;
% Now check to see the total variance is the same before and after the PCA
sigma_Optp = diag(S_Optp);                      % Singular values
TotVar_Optp = norm(ANorm_Optp,'fro')^2 ;         % Total variance of optim params
rho_Optp = norm(sigma_Optp,'fro')^2  ;                 % Variance in S matrix
% Get projection of the UCC ones
X_new_standardized = (UCC_Optp - mu) ./ sigma;
X_new_pca = X_new_standardized * V_Optp;
%% Make loadings plot
% There is something strange here; I'm not sure why. The MATLAB functions biplot 
% and quiver are reversed compared to the PCA plot, so I have to flip them.

% Adjust the direction of loadings (flip x and y axes)
loadings = V_Optp(:, 1:2);  % Flip both the first and second principal components
loadings(:, 1) = -loadings(:, 1);  % Flip the first principal component
loadings(:, 2) = -loadings(:, 2);  % Flip the second principal component (if needed)
% Create a biplot               % Getting screen size
% Consider the initial position of coordinates (flip direction as needed for the plot)
figure(88); clf; hold on;
ScrSize = get(0, 'ScreenSize');                  % Get screen size
set(gcf, 'Position', ...                         % Set window position
    [ScrSize(3) / 30 ScrSize(4) / 30 ...        
    ScrSize(3) / 1.5 ScrSize(4) / 1.5]);
% Draw arrows
quiver(zeros(size(loadings, 1), 1), zeros(size(loadings, 1), 1), ...
       loadings(:, 1), loadings(:, 2), 0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.05);
% Create biplot
biplot(loadings, 'VarLabels', Optp_Names, 'MarkerSize', 15, 'Color', 'b', 'Marker', 'none');

% Optionally adjust axis labels and tick directions
ax = gca;
ax.XDir = 'reverse';  % Use 'reverse' to flip the X-axis direction
ax.YDir = 'reverse';  % Use 'reverse' to flip the Y-axis direction
xline(0, '-', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '-', 'LineWidth', 1.5, 'Color', 'k'); hold on;

xlim([-0.48 0.52])
ylim([-0.45 0.55])
xticks = -10:5:20;
yticks = -10:5:20;
set(gca, 'Xtick', xticks, 'XTickLabel', []);
set(gca, 'Ytick', yticks, 'YTickLabel', []);
set(gca, 'FontSize', 30)
% Add axis labels
xlabel(['PC 1 (' num2str(round(100 * sigma_Optp(1)^2 / TotVar_Optp, 2)) '% variance)']);
ylabel(['PC 2 (' num2str(round(100 * sigma_Optp(2)^2 / TotVar_Optp, 2)) '% variance)']);
pbaspect([1.33 1 1]);
box on;

%% Plotting out the score (A*V) for two PCs of the optimized params

AVOptp_2PCFig = figure(1);
clf
set(gcf,'Position', ...              % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);
PCxidx = 1;
PCyidx = 2;
VarExp2PC_Optp = ...                                % Variation explained
    norm([sigma_Optp(PCxidx) sigma_Optp(PCyidx)])^2 / TotVar_Optp;          %  by first 1 and 2 PCs
% Extract PC for subgroup
Groups = {'C1','C2','C3','UCC','C4'};
% Plot UCC one first, gray
for i = 1:height(X_new_pca)
    scatter(X_new_pca(i,PCxidx), X_new_pca(i,PCyidx), 'Marker','o','SizeData',200,'MarkerFaceColor',Color_Clust{4},'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end
% Plot consistent one,other color
combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [];
for i = length(Groups) :-1 :1
    xdata = AV_Optp(ClusterResults == i,PCxidx);
    ydata = AV_Optp(ClusterResults  == i,PCyidx);
    len(i) = length(xdata);
    combinedPCX = [combinedPCX; xdata];
    combinedPCY = [combinedPCY; ydata];
    combinedColors = [combinedColors; Color_Clust{i}.* ones(length(xdata), 1)];
end

idx = randperm(len(1) + len(2) + len(3)+ len(4));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));

VarExpStr_Optp = ['Variance explained ', ...
    num2str(round(VarExp2PC_Optp,2)*100),'%'];

xline(0, '-', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '-', 'LineWidth', 1.5, 'Color', 'k'); hold on;

title('PCA Analysis','FontSize',10);
xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
xticks = -10:5:20;
yticks = -10:5:20;
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30)
box on;
pbaspect([1.33 1 1]);
text(max(AV_Optp(:,PCxidx))-0.45*xrange,max(AV_Optp(:,PCyidx))+0.1*yrange,VarExpStr_Optp,'FontSize',26, 'FontWeight','bold','Color',[0 0 0]);

%% Assign UCC one to the cluster center by knnsearch
numClusters = 3;
maxIter = 20;

% Extract only the PC1 and PC2 scores of AV_Optp for clustering
scores_AVOpt = AV_Optp(:, [PCxidx, PCyidx]);

bestClusterIdx_AVOpt = [];
bestCentroids = [];
lowestSumD = inf;

% Run K-means clustering
for i = 1:maxIter
    [clusterIdx_AVOpt, centroids, sumD] = kmeans(scores_AVOpt, numClusters, 'Replicates', 20);
    if sum(sumD) < lowestSumD
        lowestSumD = sum(sumD);
        bestClusterIdx_AVOpt = clusterIdx_AVOpt;
        bestCentroids = centroids;
    end
end

% Assign X_new_pca data to the existing clusters
X_new_pca_scores = X_new_pca(:, [PCxidx, PCyidx]);
assignedClusters_Xnew = knnsearch(bestCentroids, X_new_pca_scores);

% Plot the results
figure(3);
clf
set(gcf, 'Position', ...
    [ScrSize(3) / 30 ScrSize(4) / 30 ...
    ScrSize(3) / 1.5 ScrSize(4) / 1.5]);
Color_Clust = {[0.8500 0.3250 0.0980], [0 0.75 0.75], [0.4940 0.1840 0.5560], [0.5 0.5 0.5], [0.4660 0.6740 0.1880]};

% Plot the K-means clustering results
hold on;
for i = 1:numClusters
    scatter(scores_AVOpt(bestClusterIdx_AVOpt == i, 1), scores_AVOpt(bestClusterIdx_AVOpt == i, 2), ...
        'Marker', 'o', 'SizeData', 200, ...
        'MarkerFaceColor', Color_Clust{i}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1); % If you want to project new patients, set the MarkerFaceAlpha to a low value; otherwise, set it to a high value.
end

% Plot the assignment results for the new data
for i = 1:numClusters
    scatter(X_new_pca_scores(assignedClusters_Xnew == i, 1), X_new_pca_scores(assignedClusters_Xnew == i, 2), ...
        'Marker', 'o', 'SizeData', 200, ...
        'MarkerFaceColor', Color_Clust{i}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
end

% Plot the cluster centers
scatter(bestCentroids(:, 1), bestCentroids(:, 2), 100, 'kx', 'LineWidth', 2);

xrange = max(AV_Optp(:, PCxidx)) - min(AV_Optp(:, PCxidx));
yrange = max(AV_Optp(:, PCyidx)) - min(AV_Optp(:, PCyidx));

% Draw axis lines
xline(0, '-', 'LineWidth', 1.5, 'Color', 'k'); 
yline(0, '-', 'LineWidth', 1.5, 'Color', 'k');

% Set plot properties
xlim([min(AV_Optp(:, PCxidx)) - 0.3 * xrange, max(AV_Optp(:, PCxidx)) + 0.3 * xrange])
ylim([min(AV_Optp(:, PCyidx)) - 0.15 * yrange, max(AV_Optp(:, PCyidx)) + 0.15 * yrange])
xticks = -10:5:20;
yticks = -10:5:20;
set(gca, 'XtickLabel', []);
set(gca, 'YtickLabel', []);
set(gca, 'FontSize', 30);
box on;
pbaspect([1.33 1 1]);

%% See transition between clusters in 7 patients with 2 model windows
% Load new data from the second sheet
SecParamsT = readtable("Features.xlsx", 'VariableNamingRule', 'preserve', Sheet='Sheet2');
SecParamsT.MVr_S = (SecParamsT.MVr_S - 1) * 20;
SecParamsT.AVr_S = (SecParamsT.AVr_S - 1) * 20;
SecParamsT.PVr_S = (SecParamsT.PVr_S - 1) * 15;
SecParamsT.TVr_S = (SecParamsT.TVr_S - 1) * 17.5;

% Find unique feature positions for normalization
[~, locs] = unique(round(table2array(mean(SecParamsT(:, 3:end-3))), 4));
locs = locs + 2;
locs = sort(locs);
locs = [1; 2; locs];
SecParamsT = SecParamsT(:, locs);

% Normalize the new data
data_normalized = zscore(SecParamsT{:, 3:end});
z_threshold = 3;
outliers1 = data_normalized > z_threshold;
outliers2 = data_normalized < -z_threshold;
[~, FeatureNumber] = size(data_normalized);

% Replace outliers with max/min values
for i = 1:FeatureNumber
    SecParamsT{:, i + 2}(outliers1(:, i)) = max(SecParamsT{:, i + 2}(~outliers1(:, i)));
    SecParamsT{:, i + 2}(outliers2(:, i)) = min(SecParamsT{:, i + 2}(~outliers2(:, i)));
end

% Standardize new data (using mean and standard deviation from the first sheet)
SecData_standardized = (SecParamsT{:, slot} - mu) ./ sigma;

% Project new data into 2D PCA space
SecData_pca = SecData_standardized * V_Optp;
SecData_pca_scores = SecData_pca(:, [PCxidx, PCyidx]);
SecClusters = knnsearch(bestCentroids, SecData_pca_scores);

% Retrieve specific patient IDs from the second table
selectedPatients = SecParamsT.Patient_NO;  
ParamsT = [ParamsT; UCCParamsT]; 
clusterIdx = [bestClusterIdx_AVOpt; assignedClusters_Xnew];

% Match patients by ID between the first and second tables
[~, idxInParamsT] = ismember(selectedPatients, ParamsT.Patient_NO);

% Extract data for the matched patients
selectedPatientData = ParamsT(idxInParamsT, :);

% Standardize feature data for selected patients
selectedA_Optp = selectedPatientData{:, slot};
selectedANorm_Optp = (selectedA_Optp - mu) ./ sigma;

% Project selected patient data into PCA space
selectedAV_Optp = selectedANorm_Optp * V_Optp;

% Plot projections for data in Table 1 and Table 2, adding arrows

% Initialize figure
figure(3); % New figure
% clf;
% hold on;
ScrSize = get(0, 'ScreenSize'); 
set(gcf, 'Position', [ScrSize(3)/30 ScrSize(4)/30 ScrSize(3)/1.5 ScrSize(4)/1.5]);

% Plot data from Table 1 with colors based on clusters
for i = 1:7
    scatter(selectedAV_Optp(i, PCxidx), selectedAV_Optp(i, PCyidx), 'Marker', 'o', 'SizeData', 200, ...
        'MarkerFaceColor', Color_Clust{clusterIdx(idxInParamsT(i))}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
end

% Plot data from Table 2 with colors based on nearest cluster results
for i = 1:height(SecData_pca)
    scatter(SecData_pca(i, PCxidx), SecData_pca(i, PCyidx), ...
        'Marker', 'o', 'SizeData', 200, ...
        'MarkerFaceColor', Color_Clust{SecClusters(i)}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
end

% Find common patients between Table 1 and Table 2
[commonPatients, idxInSecParamsT, idxInParamsT] = intersect(SecParamsT.Patient_NO, ParamsT(idxInParamsT, :).Patient_NO);

% Retrieve projection coordinates for these patients in both tables
points1 = selectedAV_Optp(idxInParamsT, [PCxidx, PCyidx]);  % Projection points from Table 1
points2 = SecData_pca(idxInSecParamsT, [PCxidx, PCyidx]);   % Projection points from Table 2

% Draw arrows from points in Table 1 to Table 2
for i = 1:size(points1, 1)
    % quiver(x1, y1, deltaX, deltaY) to draw an arrow
    quiver(points1(i, 1), points1(i, 2), points2(i, 1) - points1(i, 1), points2(i, 2) - points1(i, 2), ...
        'Color', [0, 0, 0], 'LineWidth', 2, 'MaxHeadSize', 1.5); % Black arrows, adjusted head size
end

% Set axis limits and other figure properties
xrange = max(AV_Optp(:, PCxidx)) - min(AV_Optp(:, PCxidx));
yrange = max(AV_Optp(:, PCyidx)) - min(AV_Optp(:, PCyidx));

xline(0, '-', 'LineWidth', 1.5, 'Color', 'k'); 
yline(0, '-', 'LineWidth', 1.5, 'Color', 'k');

xlim([min(AV_Optp(:, PCxidx)) - 0.3 * xrange, max(AV_Optp(:, PCxidx)) + 0.3 * xrange]);
ylim([min(AV_Optp(:, PCyidx)) - 0.15 * yrange, max(AV_Optp(:, PCyidx)) + 0.15 * yrange]);

set(gca, 'Xtick', xticks, 'XTickLabel', []);
set(gca, 'Ytick', yticks, 'YTickLabel', []);
set(gca, 'FontSize', 30);
pbaspect([1.33 1 1]);
box on;