%% **********************************************************************************
%           H F p E F / H F r E F   P C A   k M E A N S   A N A L Y S I S
% ***********************************************************************************
%
%   This script takes uses SVD to perform PCA and k-Means clustering on both the
%   raw data and then the optimized parameters for heart failure patient data with
%   preserved and reduced ejection fraction (HFpEF and HFrEF).
%
%   Script created on:  16 October  2019
%   Last modified on:   10 December 2019
%
%   Developed by        Edith Jones and Brian Carlson
%                       Physiological Systems Dyanmics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************

%close all
%clc
clear all
%% Option flags for running script

%     Color_HFType = {[1 0 0], [0 1 0]};   % Colors for HFrEF and HFpEF
%     Color_HFType = {[1 0 0], [0 0.5 0]};   % Colors for HFrEF and HFpEF
Color_HFType = {[0.84,0.08,0.18], ...  % orange for HFrEF
    [0,0.45,0.74]};   % blue for HFpEF



Color_Clust = {[0.4660 0.6740 0.1880],[0 0.75 0.75],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],[0.5 0.5 0.5]};


%% Load optimized parameter values
% Load optimized parameters
ParamsT = readtable("Features.xlsx",'VariableNamingRule','preserve');
% unsupervisedslot = [(1:8) (13:23) (26:28) (37:40) (45:58)];
unsupervisedslot = (1:38); % use G
% unsupervisedslot  = [(1:4) (11:21) (24:26) (35:38) (43:48) (50:56)];

% unsupervisedslot = [(1:8) (13:23) (26:40) (45:50)];
ParamsT = ParamsT(:,unsupervisedslot);

% check out std and delete the repeat ones
[~,locs] = unique(round(table2array(mean(ParamsT(:,3:end))),4));
locs = locs+2;
locs = sort(locs);
locs = [1;2;locs];
ParamsT = ParamsT(:,locs);
% % 提取第24至27列元素的一部分
% subMatrix = ParamsT(:, [24:30]);
% 
% % 仅对这些元素取对数
% subMatrix = log2(subMatrix+1);
% 
% % 将修改后的元素放回原矩阵
% ParamsT(:, [24:30]) = subMatrix;
% replace the out of bound values by 95% confidence lower bound and upper
% bound 
data_normalized = zscore(ParamsT{:,3:end});
z_threshold = 3;
outliers1 = data_normalized > z_threshold;
outliers2 = data_normalized < -z_threshold;
[~,FeatureNumber] = size(data_normalized);
for i = 1:FeatureNumber
    ParamsT{:,i+2}(outliers1(:,i)) = max(ParamsT{:,i+2}(~outliers1(:,i)));
    ParamsT{:,i+2}(outliers2(:,i)) = min(ParamsT{:,i+2}(~outliers2(:,i)));
end

% feature selection based on mRMR 
% this is supervised it is a good approach but will not use it except for
% we want to compromise ourself
data_normalized = zscore(ParamsT{:,3:end});
[idx,scores] = fscmrmr(data_normalized,ParamsT.HFtype);
a = find(scores == 0);
figure(198);clf
bar(scores(idx))
xlabel('Predictor rank')
ylabel('Predictor importance score')
Names = ParamsT.Properties.VariableNames;
supervisedslot = idx(1:end-length(a))+2;

%%
% slot = [(3:4) 5 7 9 11 13 15 17 19 (21:27)  (28:31) 32 34 (36:38) 39 41 43 45 47 49 51 53 (55:58)]; % choose all params
supervised = 0;
if supervised == 1
    slot = supervisedslot;
else
    slot = (3:width(ParamsT));
end
 % choose which params you want to put into PCA, the ideal one should be mRMR.
% but couldn't be figured it out.

Optp_Names = ParamsT.Properties.VariableNames(slot) ;
A_Optp = ParamsT{:,slot};
HFType_Optp = ParamsT.HFtype;          % Heart failure type
PatNum_Optp = ParamsT.Patient_NO;          % Patient number

NumPats_Optp = size(A_Optp,1);                  % Number of patients

NumPatsHFrEF_Optp = 0;
NumPatsHFmrEF_Optp = 0;
NumPatsHFpEF_Optp = 0;
for i = 1:NumPats_Optp
    if ismember('HFrEF',HFType_Optp(i))
        NumPatsHFrEF_Optp = ...
            NumPatsHFrEF_Optp + 1;
    elseif ismember('HFmrEF',HFType_Optp(i))
        NumPatsHFmrEF_Optp = ...
            NumPatsHFmrEF_Optp + 1;
    elseif ismember('HFpEF',HFType_Optp(i))
        NumPatsHFpEF_Optp = ...
            NumPatsHFpEF_Optp + 1;
    end
end

for i = 1:NumPats_Optp

end
Num_Optp = size(A_Optp,2);                      % Number of optimized params


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

%% Plotting out the score (A*V) for two PCs of the optimized params
ScrSize = get(0,'ScreenSize');                  % Getting screen size
AVOptp_2PCFig = figure(1);
clf
set(gcf,'Position', ...              % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);
PCxidx = 1;
PCyidx = 3;
VarExp2PC_Optp = ...                                % Variation explained
    norm([sigma_Optp(PCxidx) sigma_Optp(PCyidx)])^2 / TotVar_Optp;          %  by first 1 and 3 PCs
% extract PC for subgroup
Groups = {'HFrEF','HFpEF'};



combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [];
for i = length(Groups) :-1 :1
    pcx = AV_Optp(strcmp(Groups{i},HFType_Optp),PCxidx);
    pcy = AV_Optp(strcmp(Groups{i},HFType_Optp),PCyidx);
    len(i) = length(pcx);
    [ellipse_x, ellipse_y] = confidence_ellipse(pcx, pcy, 0.9,1);%
    plot(ellipse_x, ellipse_y, 'LineWidth', 2, 'Color', Color_HFType{i} );hold on;
    combinedPCX = [combinedPCX; pcx];
    combinedPCY = [combinedPCY; pcy];
    combinedColors = [combinedColors; Color_HFType{i}.* ones(length(pcx), 1)];
end
idx = randperm(len(1) + len(2));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',270,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end
    % Sc(i) = scatter(pcx,pcy, ...            % Plotting pairwise
    %     'Marker','o','SizeData',200, ...           %  patient optimized
    %     'MarkerFaceColor', ...                     %  parameters
    %     Color_HFType{i}, ...
    %     'MarkerFaceAlpha', 0.75, ...
    %     'MarkerEdgeColor', ...
    %     'none'); hold on;

xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));


VarExpStr_Optp = ['Variance explained ', ...
    num2str(round(VarExp2PC_Optp,2)*100),'%'];





xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('PCA Analysis','FontSize',10);
xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
xticks = -10:5:20;
yticks = -10:5:20;
% legend(Sc,{'HFrEF', 'HFmrEF','HFpEF'},'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30)
% text(max(AV_Optp(:,PCxidx))-0.38*xrange,max(AV_Optp(:,PCyidx))+0.15*yrange,VarExpStr_Optp,'FontSize',26, 'FontWeight','bold','Color',[0 0 0]);

%% tsne
rng default % for reproducibility
[normData, meanValue, stdValue] = zscore(A_Optp);
Y = tsne(ANorm_Optp,'Algorithm','exact','Distance','euclidean');
ScrSize = get(0,'ScreenSize');                  % Getting screen size
TSNEFig = figure(2);
clf
set(gcf,'Position', ...              % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);

Groups = {'HFrEF','HFpEF'};

for i = length(Groups) :-1 :1
    xdata = Y(strcmp(Groups{i},HFType_Optp),1);
    ydata = Y(strcmp(Groups{i},HFType_Optp),2);
    Sc(i) = scatter(xdata,ydata, ...            % Plotting pairwise
        'Marker','o','SizeData',500, ...           %  patient optimized
        'MarkerFaceColor', ...                     %  parameters
        Color_HFType{i}, ...
        'MarkerFaceAlpha', 0.66, ...
        'MarkerEdgeColor', ...
        'none'); hold on;
end

xrange = max(xdata)-min(xdata);
yrange = max(ydata)-min(ydata);


% xlabel('tSNE 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('tSNE 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('t-SNE','FontSize',10);
xlim([min(xdata)-0.15*xrange  max(xdata)+0.15*xrange])
ylim([min(ydata)-0.15*yrange  max(ydata)+0.15*yrange])
xticks = -100:25:100;
yticks = -100:25:100;
% hLegend = legend(Sc,{'HFrEF', 'HFmrEF','HFpEF'},'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',15)
pbaspect([1.3,1,1]);
box on;
%% umap
pyenv('Version', 'C:\Users\fenggu\AppData\Local\Programs\Python\Python311\python.exe');
reduced_data = py.umap.UMAP(pyargs('n_neighbors',int32(15),'n_components',int32(2),'metric','euclidean')).fit_transform(py.numpy.array(normData));
Y = single(reduced_data);
ScrSize = get(0,'ScreenSize');                  % Getting screen size
UMAPFig = figure(3);
clf
set(gcf,'Position', ...              % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);

Groups = {'HFrEF','HFpEF'};

for i = length(Groups) :-1 :1
    xdata = Y(strcmp(Groups{i},HFType_Optp),1);
    ydata = Y(strcmp(Groups{i},HFType_Optp),2);
    Sc(i) = scatter(xdata,ydata, ...            % Plotting pairwise
        'Marker','o','SizeData',500, ...           %  patient optimized
        'MarkerFaceColor', ...                     %  parameters
        Color_HFType{i}, ...
        'MarkerFaceAlpha', 0.66, ...
        'MarkerEdgeColor', ...
        'none'); hold on;
end

xrange = max(xdata)-min(xdata);
yrange = max(ydata)-min(ydata);


% xlabel('UMAP1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('UMAP2', ...
%     'FontSize',30,'FontWeight','bold')

% title('UMAP','FontSize',10);
xlim([min(xdata)-0.25*xrange  max(xdata)+0.25*xrange])
ylim([min(ydata)-0.25*yrange  max(ydata)+0.25*yrange])
xticks = -100:2:100;
yticks = -100:2:100;
% hLegend = legend(Sc,{'HFrEF', 'HFmrEF','HFpEF'},'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',15)
pbaspect([1.3,1,1]);
box on;
%% DBscan to cluster the patients and visabltized
epsilon = plotKDistance(ANorm_Optp, 4);
% epsilon = 3.5;
% DBSCAN clustering
% epsilon =4; % The neighborhood size (this value needs to be determined experimentally)
minpts = 4;   % The minimum number of points required to form a cluster (this value also needs to be experimentally adjusted)

% DBSCAN execution may require the Statistics Toolbox. If absent, you might need to download the appropriate implementation.
db_labels = dbscan(ANorm_Optp, epsilon, minpts);
% Now use the 2PC plot to visualize the kmeans groups
% Creating the figure and plotting their positions
figure(11);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);


db_labels(db_labels ==-1) = 5;
Groups = {'Unclusted','DBcsan Cluster 1','DBcsan Cluster 2','DBcsan Cluster 3','DBcsan Cluster 4'};
combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [ ];
DBscanIndexMappingColor = [2 1 4 3 5];  % For example, original index 2 maps to new index 3, and so on
% Replace original indices with new indices for K-means clustering
remappedDBscanClusterIndexMapColor = Color_Clust(DBscanIndexMappingColor);
for i = length(Groups) :-1 :1
    xdata = AV_Optp(db_labels == i,PCxidx);
    ydata = AV_Optp(db_labels == i,PCyidx);
    len(i) = length(xdata);
    combinedPCX = [combinedPCX; xdata];
    combinedPCY = [combinedPCY; ydata];
    combinedColors = [combinedColors; remappedDBscanClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+ len(4)+ len(5));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('K-Means Clustering','FontSize',10);
xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
xticks = -10:5:10;
yticks = -10:5:10;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;
%% Evaluate Optimal Number of Clusters for Kmeans and Hierarchical Clustering
% find a reasonable cluster number
rng("default") % Set the random number generator to default for reproducibility
maxClusters = 100;  % Maximum number of clusters to evaluate

% K-Means Clustering Evaluations
% Calculate Within-Cluster Sum of Square (elbow method) for Kmeans
wcssEvaluation = zeros(maxClusters, 1);
for k = 1:maxClusters
    [~, ~, sumd] = kmeans(ANorm_Optp, k);
    wcssEvaluation(k) = sum(sumd);
end

% Calculate silhouette scores for Kmeans
silhouetteEvaluation = evalclusters(ANorm_Optp, "kmeans", "silhouette", "KList", 1:maxClusters);
% Calculate Davies-Bouldin Index for Kmeans
dbiEvaluation = evalclusters(ANorm_Optp, "kmeans", "DaviesBouldin", "KList", 1:maxClusters);
% Calculate Calinski-Harabasz Index for Kmeans
chiEvaluation = evalclusters(ANorm_Optp, "kmeans", "CalinskiHarabasz", "KList", 1:maxClusters);

% Hierarchical Clustering Evaluations
% Calculate silhouette scores for Hierarchical Clustering
silhouetteHierEvaluation = evalclusters(ANorm_Optp, "linkage", "silhouette", "KList", 1:maxClusters);
% Calculate Davies-Bouldin Index for Hierarchical Clustering
dbiHierEvaluation = evalclusters(ANorm_Optp, "linkage", "DaviesBouldin", "KList", 1:maxClusters);
% Calculate Calinski-Harabasz Index for Hierarchical Clustering
chiHierEvaluation = evalclusters(ANorm_Optp, "linkage", "CalinskiHarabasz", "KList", 1:maxClusters);
range = [0 20];
optimalNumClusters = 4;
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

%% Grouping into clusters using kmeans on optimized parameters

rng('default')
[KMnsClust_Optp,KMnsCntrd_Optp, ...                 % KMeans clust numbers
    KMnsSumDst_Optp] = ...                          %  centroids, sum of dist
    kmeans(ANorm_Optp,4, ...            %  Norm data, # of clusts
    'Distance','sqeuclidean', ...                     %  Use L1 distance
    'Display','final', ...                          %  Display results
    'Replicates',20);                              %  Run 200 times
for i  = 1:20
    if abs(mean(KMnsClust_Optp(strcmp('HFpEF',HFType_Optp)))-2) < 0.5     
    else
        [KMnsClust_Optp,KMnsCntrd_Optp, ...
            KMnsSumDst_Optp] = ...
            kmeans(ANorm_Optp,4, ...
            'Distance','sqeuclidean', ...
            'Display','final', ...
            'Replicates',20);
    end
end

% Assign cluster back to table
ParamsT.KeamsCluster = KMnsClust_Optp;
structname = ParamsT.Properties.VariableNames;
for i  = 3:length(structname)
    slot = ParamsT.(structname{i});
    DataC1 = slot(ParamsT.KeamsCluster==1);
    DataC2 = slot(ParamsT.KeamsCluster==2);
    DataC3 = slot(ParamsT.KeamsCluster==3);
    DataC4 = slot(ParamsT.KeamsCluster==4);
    DataC5 = slot(ParamsT.KeamsCluster==5);
    H = max([length(DataC1) length(DataC2) length(DataC3) length(DataC4) length(DataC5)]);
    newT  = table(nan(H,1),nan(H,1),nan(H,1),nan(H,1),nan(H,1),...
        'VariableNames', {'C1', 'C2', 'C3','C4','C5'});
    newT.C1(1:length(DataC1)) = DataC1;
    newT.C2(1:length(DataC2)) = DataC2;
    newT.C3(1:length(DataC3)) = DataC3;
    newT.C4(1:length(DataC4)) = DataC4;
    newT.C5(1:length(DataC5)) = DataC5;
    ViolinDataKmeans.(structname{i}) = newT;
end
%
% Now use the 2PC plot to visualize the kmeans groups
% Creating the figure and plotting their positions
AVOptp_2PCKMnsFig = figure(12);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);



Groups = {'Kmeans Cluster 1','Kmeans Cluster 2','Kmeans Cluster 3','Kmeans Cluster 4'};
combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [ ];
KmeansIndexMappingColor = [2 4 3 1];  % For example, original index 2 maps to new index 3, and so on
% Replace original indices with new indices for K-means clustering
remappedKmeansClusterIndexMapColor = Color_Clust(KmeansIndexMappingColor);
KmsIndexRenumber = [2 1 4 3];
NewKMnsClust_Optp = arrayfun(@(x) KmsIndexRenumber(x), KMnsClust_Optp);

for i = length(Groups) :-1 :1
    xdata = AV_Optp(KMnsClust_Optp == i,PCxidx);
    ydata = AV_Optp(KMnsClust_Optp == i,PCyidx);
    len(i) = length(xdata);
    combinedPCX = [combinedPCX; xdata];
    combinedPCY = [combinedPCY; ydata];
    combinedColors = [combinedColors; remappedKmeansClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+ len(4));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('K-Means Clustering','FontSize',10);
xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
xticks = -10:5:10;
yticks = -10:5:10;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;

%% Grouping into clusters using hierachical clustering on optimized parameters

% Creating hierachical tree where we use the Ward metric
%  when joining clusters which considers the increase in
%  the within cluster sum of squares. This is just the
%  squared distance between each element in the cluster
%  and it's centroid
HCLink_Optp = linkage(ANorm_Optp,'ward');           % Making the tree
HCClust_Optp = cluster(HCLink_Optp, ...             % Cluster number to display
    'Maxclust',4);

% Assign cluster back to table
ParamsT.HCCCluster = HCClust_Optp;
structname = ParamsT.Properties.VariableNames;
for i  = 3:length(structname)
    slot = ParamsT.(structname{i});
    DataC1 = slot(ParamsT.HCCCluster==1);
    DataC2 = slot(ParamsT.HCCCluster==2);
    DataC3 = slot(ParamsT.HCCCluster==3);
    DataC4 = slot(ParamsT.HCCCluster==4);
    DataC5 = slot(ParamsT.HCCCluster==5);
    H = max([length(DataC1) length(DataC2) length(DataC3) length(DataC4) length(DataC5)]);
    newT  = table(nan(H,1),nan(H,1),nan(H,1),nan(H,1),nan(H,1),...
        'VariableNames', {'C1', 'C2', 'C3','C4','C5'});
    newT.C1(1:length(DataC1)) = DataC1;
    newT.C2(1:length(DataC2)) = DataC2;
    newT.C3(1:length(DataC3)) = DataC3;
    newT.C4(1:length(DataC4)) = DataC4;
    newT.C5(1:length(DataC5)) = DataC5;
    ViolinDataHCC.(structname{i}) = newT;
end
%% 
%  third-from-last linkages.
HCCutOff_Optp = median([HCLink_Optp(end-3,3) ...
    HCLink_Optp(end-2,3)]);
% % assign the color I want
% H = dendrogram(HCLink_Optp, 0, 'ColorThreshold', HCCutOff_Optp, 'Labels', num2str(PatNum_Optp));
% for i = 1:length(H)
%     x = get(H(i), 'XData');
%     THier = sort(HCClust_Optp);
%     group = THier(round(x(1))); 
%     set(H(i), 'Color', Color_HCCClust{group}); 
% end



%%heatmap

cgo=clustergram (ANorm_Optp,'Linkage','ward','Cluster',1, 'Dendrogram',HCCutOff_Optp,'Colormap',redbluecmap);
set(cgo,'RowLabels',PatNum_Optp,'ColumnLabels',{'C_S_A','C_S_V','C_P_A','C_P_V','G_S_A','G_P_A','G_t_S_A','G_t_P_A','kact_L_V','kact_R_V','kpas_L_V','kpas_R_V',...
     'Vw_L_V','Vw_R_V','Vw_S_E_P','V_a_u','V_a_c','G_a','V4c','K_p_e_r_i','B_p_e_r_i','G_t_o','G_t_c','G_p_o','G_p_c','G_m_o','G_m_c','G_a_o','G_a_c','V_v_s'},'ColumnLabelsRotate', 45)



%%
% Now use the 2PC plot to visualize the HCC groups
% Creating the figure and plotting their positions
AVOptp_2PCHCFig = figure(13);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);


Groups = {'Hierchical Cluster 1','Hierchical Cluster 2','Hierchical Cluster 3','Hierchical Cluster 4'};

combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [ ];
HierIndexMappingColor = [3 1 2 4];  % For example, original index 2 maps to new index 3, and so on
remappedHierClusterIndexMapColor = Color_Clust(HierIndexMappingColor);
HierIndexRenumber = [4 3 2 1];
NewHierClust_Optp = arrayfun(@(x) HierIndexRenumber(x), HCClust_Optp);
for i = length(Groups) :-1 :1
    xdata = AV_Optp(HCClust_Optp == i,PCxidx);
    ydata = AV_Optp(HCClust_Optp == i,PCyidx);
    len(i) = length(xdata);
    combinedPCX = [combinedPCX; xdata];
    combinedPCY = [combinedPCY; ydata];
    combinedColors = [combinedColors; remappedHierClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+len(4));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;

xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('Hierchical Clustering','FontSize',10);
xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
xticks = -10:5:10;
yticks = -10:5:10;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;

%% Initialize risk level column array
SumClusters = [NewKMnsClust_Optp NewHierClust_Optp];
riskLevel = zeros(size(SumClusters, 1), 1);  % Default is medium risk, represented by 2

% 重新计算 riskLevel
for i = 1:size(SumClusters, 1)
    if all(SumClusters(i, :) == 3 | SumClusters(i, :) == 4) % 如果一行里都是3或者4
        riskLevel(i) = 1;  % 低风险，表示为 1
    elseif all(SumClusters(i, :) == 2) % 如果一行里都是2
        riskLevel(i) = 2;  % 中等风险，表示为 2
    elseif all(SumClusters(i, :) == 1) % 如果一行里都是1
        riskLevel(i) = 3;  % 高风险，表示为 3
    elseif any(SumClusters(i, :) == 1) && any(SumClusters(i, :) == 2) % 如果一行有1和有2
        riskLevel(i) = 3;  % 高风险，表示为 3
    elseif (any(SumClusters(i, :) == 3) || any(SumClusters(i, :) == 4)) && any(SumClusters(i, :) == 2) % 如果一行有3或4，并且有2
        riskLevel(i) = 1;  % 低风险，表示为 1
    elseif any(SumClusters(i, :) == 1) && (any(SumClusters(i, :) == 3) || any(SumClusters(i, :) == 4)) % 如果一行有1，并且有3或4
        riskLevel(i) = 2;  % 中等风险，表示为 2
    end
end

% Add the risk level array to SumClusters matrix to create a new table
SumClustersWithRisk = [SumClusters, riskLevel];

%% calculate percentange of patients in each group
HierC1 = HFType_Optp( HCClust_Optp==1 );
HierC2 = HFType_Optp( HCClust_Optp==2 );
HierC3 = HFType_Optp( HCClust_Optp==3 );
HierC4 = HFType_Optp( HCClust_Optp==4 );
HierC5 = HFType_Optp( HCClust_Optp==5 );
C1portion = [length(HierC1(strcmp('HFrEF',HierC1))) length(HierC1(strcmp('HFmrEF',HierC1))) length(HierC1(strcmp('HFpEF',HierC1)))];
C2portion = [length(HierC2(strcmp('HFrEF',HierC2))) length(HierC2(strcmp('HFmrEF',HierC2))) length(HierC2(strcmp('HFpEF',HierC2)))];
C3portion = [length(HierC3(strcmp('HFrEF',HierC3))) length(HierC3(strcmp('HFmrEF',HierC3))) length(HierC3(strcmp('HFpEF',HierC3)))];
C4portion = [length(HierC4(strcmp('HFrEF',HierC4))) length(HierC4(strcmp('HFmrEF',HierC4))) length(HierC4(strcmp('HFpEF',HierC4)))];
C5portion = [length(HierC5(strcmp('HFrEF',HierC5))) length(HierC5(strcmp('HFmrEF',HierC5))) length(HierC5(strcmp('HFpEF',HierC5)))];
THier = array2table([C1portion;C2portion;C3portion;C4portion;C5portion]);
THier.Properties.VariableNames = {'HFrEF','HFmrEF','HFpEF'};
THier.Properties.RowNames = {'C1','C2','C3','C4','C5'};
KmsC1 = HFType_Optp( KMnsClust_Optp==1 );
KmsC2 = HFType_Optp( KMnsClust_Optp==2 );
KmsC3 = HFType_Optp( KMnsClust_Optp==3 );
C1portion = [length(KmsC1(strcmp('HFrEF',KmsC1))) length(KmsC1(strcmp('HFmrEF',KmsC1))) length(KmsC1(strcmp('HFpEF',KmsC1)))];
C2portion = [length(KmsC2(strcmp('HFrEF',KmsC2))) length(KmsC2(strcmp('HFmrEF',KmsC2))) length(KmsC2(strcmp('HFpEF',KmsC2)))];
C3portion = [length(KmsC3(strcmp('HFrEF',KmsC3))) length(KmsC3(strcmp('HFmrEF',KmsC3))) length(KmsC3(strcmp('HFpEF',KmsC3)))];
TKms = array2table([C1portion;C2portion;C3portion]);
TKms.Properties.VariableNames = {'HFrEF','HFmrEF','HFpEF'};
TKms.Properties.RowNames = {'C1','C2','C3'};


function [ellipse_x, ellipse_y] = confidence_ellipse(x, y, confidence_level, adjust_factor)

s = chi2inv(confidence_level, 2);

mu = [mean(x) mean(y)];

sigma = cov(x, y);

[eigenvec, eigenval] = eig(sigma);

[max_eigenval, max_eigenvec_index] = max(diag(eigenval));
if max_eigenvec_index == 2
    eigenvec = fliplr(eigenvec);
end

a = sqrt(max_eigenval * s); % 
b = sqrt(min(diag(eigenval)) * s); % 


angle = atan2(eigenvec(2,1), eigenvec(1,1));


theta_grid = linspace(0, 2*pi);
ellipse_x_r = a * cos( theta_grid );
ellipse_y_r = b * sin( theta_grid );

ellipse_y_r = ellipse_y_r * adjust_factor;


R = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];
r_ellipse = [ellipse_x_r; ellipse_y_r]' * R;


if mu(2) < -0.5
    ellipse_y = r_ellipse(:,2) + mu(2) ;
    ellipse_x = r_ellipse(:,1) + mu(1);
elseif mu(2) > 0.5
    ellipse_y = r_ellipse(:,2) + mu(2);
    ellipse_x = r_ellipse(:,1) + mu(1)-0.5;
else
    ellipse_y = r_ellipse(:,2) + mu(2);
    ellipse_x = r_ellipse(:,1) + mu(1);
end

end





