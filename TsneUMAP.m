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



Color_Clust = {[0.4940 0.1840 0.5560], [0 0.75 0.75], [0.8500 0.3250 0.0980], [0.5 0.5 0.5], [0.4660 0.6740 0.1880]};

% 新的对比强烈的颜色
Color_Clust{end+1} = [1.0 0.84 0.0]; % 金黄色
Color_Clust{end+1} = [0.0 1.0 0.0]; % 绿色

%% Load optimized parameter values
% Load optimized parameters
ParamsT = readtable("Features.xlsx",'VariableNamingRule','preserve');
ParamsT.MVr_S = (ParamsT.MVr_S-1).*20;
ParamsT.AVr_S = (ParamsT.AVr_S-1).*20;
ParamsT.PVr_S = (ParamsT.PVr_S-1).*15;
ParamsT.TVr_S = (ParamsT.TVr_S-1).*17.5;
ParamsT = ParamsT(:,1:width(ParamsT)-1);
% AlldataT = readtable("Predictors.xlsx", 'VariableNamingRule', 'preserve');
% unsupervisedslot = [(1:8) (13:23) (26:28) (37:40) (45:58)];
% unsupervisedslot = (1:38); % use G
% unsupervisedslot  = [(1:30) 32 34 36];
% ParamsT = ParamsT(:, {'Patient_NO','HFtype','R_m_o','k_pas_LV','K_P','Vw_RV','Vh0','RAV0u','RAV0c','Vw_SEP',...
%     'Vw_LV','B_P','R_RA','V_SV_s','k_act_LV','C_PA','C_PV','C_SV','MVr_S','TVr_S','LVpowerIndex_S',...
%     'LVMD_S','StrainRV_S','StrainLV_S','StrainSEP_S','StressRV_S','RVpower_S','RVpowerIndex_S'});

% unsupervisedslot = [(1:8) (13:23) (26:40) (45:50)];
% ParamsT = ParamsT(:,unsupervisedslot);
% shuntlist = [34 41 54 61 83 116 183 231 268 278 312];
% AlldataT(shuntlist,:) = [];
% ParamsT(shuntlist,:) = [];
% d_columns = find(endsWith(AlldataT.Properties.VariableNames, '_D'));
% if ~isempty(d_columns)
%     first_d_col = d_columns(1);
%     last_d_col = d_columns(end);
%     ParamsT = AlldataT(:,  [1,4,8,(first_d_col:last_d_col)]);
% end
% only feed RHC MRI AND ECHO AND MRI
% check out std and delete the repeat ones
[~,locs] = unique(round(table2array(mean(ParamsT(:,3:end))),4));
locs = locs+2;
locs = sort(locs);
locs = [1;2;locs];
ParamsT = ParamsT(:,locs);
% ParamsT.R_m_c = log(ParamsT.R_m_c)/log(1e3);
% ParamsT.R_a_c = log(ParamsT.R_a_c)/log(1e3);
% ParamsT.R_t_c = log(ParamsT.R_t_c)/log(1e3);
% ParamsT.R_p_c = log(ParamsT.R_p_c)/log(1e3);

data_normalized = zscore(ParamsT{:,3:end});
z_threshold = 3;
outliers1 = data_normalized > z_threshold;
outliers2 = data_normalized < -z_threshold;
[~,FeatureNumber] = size(data_normalized);
for i = 1:FeatureNumber
    ParamsT{:,i+2}(outliers1(:,i)) = max(ParamsT{:,i+2}(~outliers1(:,i)));
    ParamsT{:,i+2}(outliers2(:,i)) = min(ParamsT{:,i+2}(~outliers2(:,i)));
end
% % Separate numeric and non-numeric columns
% numericCols = varfun(@isnumeric, ParamsT, 'OutputFormat', 'uniform');
% numericParamT = ParamsT(:, numericCols);
% nonNumericParamT = ParamsT(:, ~numericCols);
% 
% % Remove columns with more than 25% missing values from numeric columns
% missingThreshold = 0.25;
% colsToRemove = varfun(@(x) mean(isnan(x)) > missingThreshold, numericParamT, 'OutputFormat', 'uniform');
% numericParamT(:, colsToRemove) = [];
% 
% % Use built-in kNN imputation on numeric columns
% numericParamTArray = table2array(numericParamT);
% numericParamTArray = knnimpute(numericParamTArray);
% numericParamT = array2table(numericParamTArray, 'VariableNames', numericParamT.Properties.VariableNames);
% 
% % Combine numeric and non-numeric columns back into one table
% ParamsT = [nonNumericParamT, numericParamT];
% if ismember('PatID', ParamsT.Properties.VariableNames)
%     ParamsT.Properties.VariableNames{'PatID'} = 'Patient_NO';
% else
%     error('Table does not contain a column named ''patID''.');
% end
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
PCyidx = 2;
VarExp2PC_Optp = ...                                % Variation explained
    norm([sigma_Optp(PCxidx) sigma_Optp(PCyidx)])^2 / TotVar_Optp;          %  by first 1 and 3 PCs
% extract PC for subgroup
Groups = {'HFrEF','HFpEF'};



combinedX = [ ];
combinedY = [ ];
combinedColors = [];
for i = length(Groups) :-1 :1
    pcx = AV_Optp(strcmp(Groups{i},HFType_Optp),PCxidx);
    pcy = AV_Optp(strcmp(Groups{i},HFType_Optp),PCyidx);
    len(i) = length(pcx);
    [ellipse_x, ellipse_y] = confidence_ellipse(pcx, pcy, 0.9,1);%
    plot(ellipse_x, ellipse_y, 'LineWidth', 2, 'Color', Color_HFType{i} );hold on;
    combinedX = [combinedX; pcx];
    combinedY = [combinedY; pcy];
    combinedColors = [combinedColors; Color_HFType{i}.* ones(length(pcx), 1)];
end
idx = randperm(len(1) + len(2));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',270,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
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
Ytsne = tsne(ANorm_Optp,'Algorithm','exact','Distance','euclidean');
ScrSize = get(0,'ScreenSize');                  % Getting screen size
TSNEFig = figure(2);
clf
set(gcf,'Position', ...              % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);

Groups = {'HFrEF','HFpEF'};

for i = length(Groups) :-1 :1
    xdata = Ytsne(strcmp(Groups{i},HFType_Optp),1);
    ydata = Ytsne(strcmp(Groups{i},HFType_Optp),2);
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
seed = int32(8);
pyenv('Version', 'C:\Users\fenggu\AppData\Local\Programs\Python\Python311\python.exe');
reduced_data = py.umap.UMAP(pyargs('n_neighbors',int32(15),'n_components',int32(2),'metric','euclidean','random_state', seed, 'n_jobs', int32(1))).fit_transform(py.numpy.array(normData));
Yumap = single(reduced_data);
ScrSize = get(0,'ScreenSize');                  % Getting screen size
UMAPFig = figure(3);
clf
set(gcf,'Position', ...              % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);

Groups = {'HFrEF','HFpEF'};

for i = length(Groups) :-1 :1
    xdata = Yumap(strcmp(Groups{i},HFType_Optp),1);
    ydata = Yumap(strcmp(Groups{i},HFType_Optp),2);
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


xlabel('UMAP1', ...                 % Labelling axes with
    'FontSize',30,'FontWeight','bold')              %  PC labels
ylabel('UMAP2', ...
    'FontSize',30,'FontWeight','bold')

title('UMAP','FontSize',10);
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
% %% DBscan to cluster the patients and visabltized
% epsilon = plotKDistance(ANorm_Optp, 4);
% % epsilon = 3.5;
% % DBSCAN clustering
% % epsilon =4; % The neighborhood size (this value needs to be determined experimentally)
% minpts = 4;   % The minimum number of points required to form a cluster (this value also needs to be experimentally adjusted)
% 
% % DBSCAN execution may require the Statistics Toolbox. If absent, you might need to download the appropriate implementation.
% db_labels = dbscan(ANorm_Optp, epsilon, minpts);
% % Now use the 2PC plot to visualize the kmeans groups
% % Creating the figure and plotting their positions
% figure(11);
% clf
% set(gcf,'Position', ...          % Positioning the figure
%     [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
%     ScrSize(3)/1.5 ScrSize(4)/1.5]);
% 
% 
% db_labels(db_labels ==-1) = 5;
% Groups = {'Unclusted','DBcsan Cluster 1','DBcsan Cluster 2','DBcsan Cluster 3','DBcsan Cluster 4'};
% combinedPCX = [ ];
% combinedPCY = [ ];
% combinedColors = [ ];
% DBscanIndexMappingColor = [2 1 4 3 5];  % For example, original index 2 maps to new index 3, and so on
% % Replace original indices with new indices for K-means clustering
% remappedDBscanClusterIndexMapColor = Color_Clust(DBscanIndexMappingColor);
% for i = length(Groups) :-1 :1
%     xdata = AV_Optp(db_labels == i,PCxidx);
%     ydata = AV_Optp(db_labels == i,PCyidx);
%     len(i) = length(xdata);
%     combinedPCX = [combinedPCX; xdata];
%     combinedPCY = [combinedPCY; ydata];
%     combinedColors = [combinedColors; remappedDBscanClusterIndexMapColor{i}.* ones(length(xdata), 1)];
% end
% idx = randperm(len(1) + len(2) + len(3)+ len(4)+ len(5));
% for i = 1:length(idx)
%     scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
% end
% 
% xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
% yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));
% 
% xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% % xlabel('Principal Component 1', ...                 % Labelling axes with
% %     'FontSize',30,'FontWeight','bold')              %  PC labels
% % ylabel('Principal Component 2', ...
% %     'FontSize',30,'FontWeight','bold')
% 
% % title('K-Means Clustering','FontSize',10);
% xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
% ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
% xticks = -10:5:10;
% yticks = -10:5:10;
% % legend(Sc,Groups,'Location','southwest');
% % legend(Box="off");
% set(gca,'Xtick',xticks,'XTickLabel',[]);
% set(gca,'Ytick',yticks,'YTickLabel',[]);
% set(gca,'FontSize',30);box on;
% hold off;
%% Evaluate Optimal Number of Clusters for Kmeans and Hierarchical Clustering
% % Set maximum number of clusters
% maxClusters = 50;
% 
% % Evaluate K-means clustering
% silhouetteEvaluationKmeans = evalclusters(ANorm_Optp, "kmeans", "silhouette", "KList", 1:maxClusters);
% dbiEvaluationKmeans = evalclusters(ANorm_Optp, "kmeans", "DaviesBouldin", "KList", 1:maxClusters);
% chEvaluationKmeans = evalclusters(ANorm_Optp, "kmeans", "CalinskiHarabasz", "KList", 1:maxClusters);
% gapEvaluationKmeans = evalclusters(ANorm_Optp, "kmeans", "gap", "KList", 1:maxClusters);
% 
% % Evaluate Hierarchical clustering
% silhouetteEvaluationHier = evalclusters(ANorm_Optp, "linkage", "silhouette", "KList", 1:maxClusters);
% dbiEvaluationHier = evalclusters(ANorm_Optp, "linkage", "DaviesBouldin", "KList", 1:maxClusters);
% chEvaluationHier = evalclusters(ANorm_Optp, "linkage", "CalinskiHarabasz", "KList", 1:maxClusters);
% gapEvaluationHier = evalclusters(ANorm_Optp, "linkage", "gap", "KList", 1:maxClusters);
% 
% % Evaluate Gaussian Mixture Model clustering
% silhouetteEvaluationGMM = evalclusters(ANorm_Optp, "gmdistribution", "silhouette", "KList", 1:maxClusters);
% dbiEvaluationGMM = evalclusters(ANorm_Optp, "gmdistribution", "DaviesBouldin", "KList", 1:maxClusters);
% chEvaluationGMM = evalclusters(ANorm_Optp, "gmdistribution", "CalinskiHarabasz", "KList", 1:maxClusters);
% gapEvaluationGMM = evalclusters(ANorm_Optp, "gmdistribution", "gap", "KList", 1:maxClusters);
% 
% %% Plot the evaluation metrics to assist with determining the optimal number of clusters
% figure(101);clf;
% range = [0 25];
% fontSize = 14;  
% 
% % K-means Clustering
% % subplot(3, 4, 1);
% % plot(1:maxClusters, silhouetteEvaluationKmeans.CriterionValues, 'b-o');
% % title('K-means Silhouette', 'FontSize', fontSize);
% % xlabel('Number of Clusters', 'FontSize', fontSize);
% % ylabel('Silhouette Score', 'FontSize', fontSize);
% % xlim(range);
% 
% subplot(3, 2, 1);
% plot(1:maxClusters, dbiEvaluationKmeans.CriterionValues, 'b-o');hold on;
% title('K-means Davies-Bouldin', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('DBI Score', 'FontSize', fontSize);
% xlim(range);
% scatter(3,dbiEvaluationKmeans.CriterionValues(3),10,"red","filled")
% % subplot(3, 4, 3);
% % plot(1:maxClusters, chEvaluationKmeans.CriterionValues, 'b-o');
% % title('K-means Calinski-Harabasz', 'FontSize', fontSize);
% % xlabel('Number of Clusters', 'FontSize', fontSize);
% % ylabel('Calinski-Harabasz', 'FontSize', fontSize);
% % xlim(range);
% 
% subplot(3, 2, 2);
% plot(1:maxClusters, gapEvaluationKmeans.CriterionValues, 'b-o');hold on;
% title('K-means Gap Statistic', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('Gap Value', 'FontSize', fontSize);
% xlim(range);
% scatter(3,gapEvaluationKmeans.CriterionValues(3),10,"red","filled")
% % Hierarchical Clustering
% % subplot(3, 4, 5);
% % plot(1:maxClusters, silhouetteEvaluationHier.CriterionValues, 'b-o');
% % title('Hierarchical Silhouette', 'FontSize', fontSize);
% % xlabel('Number of Clusters', 'FontSize', fontSize);
% % ylabel('Silhouette Score', 'FontSize', fontSize);
% % xlim(range);
% 
% subplot(3, 2, 3);
% plot(1:maxClusters, dbiEvaluationHier.CriterionValues, 'b-o');hold on;
% title('Hierarchical Davies-Bouldin', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('DBI Score', 'FontSize', fontSize);
% xlim(range);
% scatter(4,dbiEvaluationHier.CriterionValues(4),10,"red","filled")
% % subplot(3, 4, 7);
% % plot(1:maxClusters, chEvaluationHier.CriterionValues, 'b-o');
% % title('Hierarchical Calinski-Harabasz', 'FontSize', fontSize);
% % xlabel('Number of Clusters', 'FontSize', fontSize);
% % ylabel('Calinski-Harabasz', 'FontSize', fontSize);
% % xlim(range);
% 
% subplot(3, 2, 4);
% plot(1:maxClusters, gapEvaluationHier.CriterionValues, 'b-o');hold on;
% title('Hierarchical Gap Statistic', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('Gap Value', 'FontSize', fontSize);
% xlim(range);
% scatter(3,gapEvaluationHier.CriterionValues(3),10,"red","filled")
% % GMM Clustering
% % subplot(3, 4, 9);
% % plot(1:maxClusters, silhouetteEvaluationGMM.CriterionValues, 'b-o');
% % title('GMM Silhouette', 'FontSize', fontSize);
% % xlabel('Number of Clusters', 'FontSize', fontSize);
% % ylabel('Silhouette Score', 'FontSize', fontSize);
% % xlim(range);
% 
% subplot(3, 2, 5);
% plot(1:maxClusters, dbiEvaluationGMM.CriterionValues, 'b-o');hold on;
% title('GMM Davies-Bouldin', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('DBI Score', 'FontSize', fontSize);
% xlim(range);
% scatter(3,dbiEvaluationGMM.CriterionValues(3),10,"red","filled")
% % subplot(3, 4, 11);
% % plot(1:maxClusters, chEvaluationGMM.CriterionValues, 'b-o');
% % title('GMM Calinski-Harabasz', 'FontSize', fontSize);
% % xlabel('Number of Clusters', 'FontSize', fontSize);
% % ylabel('Calinski-Harabasz', 'FontSize', fontSize);
% % xlim(range);
% 
% subplot(3, 2, 6);
% plot(1:maxClusters, gapEvaluationGMM.CriterionValues, 'b-o');hold on;
% title('GMM Gap Statistic', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('Gap Value', 'FontSize', fontSize);
% xlim(range);
% scatter(4,gapEvaluationGMM.CriterionValues(4),10,"red","filled")
% % Adjust the layout to prevent labels from overlapping
% set(gcf, 'Position', [100, 100, 1200, 800]);  % Resize figure to make it wider


%% Grouping into clusters using kmeans on optimized parameters

rng('default')
[KMnsClust_Optp,KMnsCntrd_Optp, ...                 % KMeans clust numbers
    KMnsSumDst_Optp] = ...                          %  centroids, sum of dist
    kmeans(ANorm_Optp,3, ...            %  Norm data, # of clusts
    'Distance','sqeuclidean', ...                     %  Use L1 distance
    'Display','final', ...                          %  Display results
    'Replicates',20);                              %  Run 200 times
for i  = 1:20
    if abs(mean(KMnsClust_Optp(strcmp('HFpEF',HFType_Optp)))-2) < 0.5     
    else
        [KMnsClust_Optp,KMnsCntrd_Optp, ...
            KMnsSumDst_Optp] = ...
            kmeans(ANorm_Optp,3, ...
            'Distance','sqeuclidean', ...
            'Display','final', ...
            'Replicates',20);
    end
end
TsneKMnsFig = figure(12);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);



Groups = {'Kmeans Cluster 1','Kmeans Cluster 2','Kmeans Cluster 3','Kmeans Cluster 4','Kmeans Cluster 5'};
combinedX = [ ];
combinedY = [ ];
combinedColors = [ ];
KmeansIndexMappingColor = [1 2 3 4 5];  % For example, original index 2 maps to new index 3, and so on
% Replace original indices with new indices for K-means clustering
remappedKmeansClusterIndexMapColor = Color_Clust(KmeansIndexMappingColor);
KmsIndexRenumber = [1 2 3 4 5];
NewKMnsClust_Optp = arrayfun(@(x) KmsIndexRenumber(x), KMnsClust_Optp);

for i = length(Groups) :-1 :1
    xdata = Ytsne(KMnsClust_Optp == i,1);
    ydata = Ytsne(KMnsClust_Optp == i,2);
    len(i) = length(xdata);
    combinedX = [combinedX; xdata];
    combinedY = [combinedY; ydata];
    combinedColors = [combinedColors; remappedKmeansClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+ len(4) + len(5));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(Ytsne(:,1))-min(Ytsne(:,1));
yrange = max(Ytsne(:,2))-min(Ytsne(:,2));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('K-Means Clustering','FontSize',10);
xlim([min(Ytsne(:,1))-0.3*xrange  max(Ytsne(:,1))+0.3*xrange])
ylim([min(Ytsne(:,2))-0.15*yrange  max(Ytsne(:,2))+0.15*yrange])
xticks = -100:25:100;
yticks = -100:25:100;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;

UmapKMnsFig = figure(22);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);



combinedX = [ ];
combinedY = [ ];
combinedColors = [ ];

for i = length(Groups) :-1 :1
    xdata = Yumap(KMnsClust_Optp == i,1);
    ydata = Yumap(KMnsClust_Optp == i,2);
    len(i) = length(xdata);
    combinedX = [combinedX; xdata];
    combinedY = [combinedY; ydata];
    combinedColors = [combinedColors; remappedKmeansClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+ len(4) + len(5));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(Yumap(:,1))-min(Yumap(:,1));
yrange = max(Yumap(:,2))-min(Yumap(:,2));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('K-Means Clustering','FontSize',10);
xlim([min(Yumap(:,1))-0.3*xrange  max(Yumap(:,1))+0.3*xrange])
ylim([min(Yumap(:,2))-0.15*yrange  max(Yumap(:,2))+0.15*yrange])
xticks = -100:2:100;
yticks = -100:2:100;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;
%% ISODATA
% [Z, ~, A, KMnsClust_Optp]=wfIsodata_ND(ANorm_Optp, 10, 2, 2000, 1, 50, 1, 1, 100);
% KMnsClust_Optp = KMnsClust_Optp';


%% Grouping into clusters using hierachical clustering on optimized parameters

% Creating hierachical tree where we use the Ward metric
%  when joining clusters which considers the increase in
%  the within cluster sum of squares. This is just the
%  squared distance between each element in the cluster
%  and it's centroid
HCLink_Optp = linkage(ANorm_Optp,'ward');           % Making the tree
HCClust_Optp = cluster(HCLink_Optp, ...             % Cluster number to display
    'Maxclust',3);
%
%  third-from-last linkages.
HCCutOff_Optp = median([HCLink_Optp(end-1,3) ...
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
% set(cgo,'RowLabels',PatNum_Optp,'ColumnLabels',{'C_S_A','C_S_V','C_P_A','C_P_V','kact_L_V','kact_R_V','kpas_L_V','kpas_R_V',...
%      'Vw_L_V','Vw_R_V','Vw_S_E_P','V_a_u','V_a_c','V_4_C','K_P_C','B_P_C','V_V_S','R_S_A','R_P_A','R_t_S_A','R_t_P_A','R_A_t','R_t_o','R_t_c','R_p_o','R_p_c','R_m_o','R_m_c','R_a_o','R_a_c',},'ColumnLabelsRotate', 45)

set(cgo,'RowLabels',PatNum_Optp,'ColumnLabels',Names(3:end),'ColumnLabelsRotate', 45)

TsneHCFig = figure(13);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);


Groups = {'Hierchical Cluster 1','Hierchical Cluster 2','Hierchical Cluster 3','Hierchical Cluster 4'};

combinedX = [ ];
combinedY = [ ];
combinedColors = [ ];
% HierIndexMappingColor = [2 4 1 3];  % For example, original index 2 maps to new index 3, and so on
% remappedHierClusterIndexMapColor = Color_Clust(HierIndexMappingColor);
% HierIndexRenumber = [2 4 1 3];
HierIndexMappingColor = [2 1 3 4];  % For example, original index 2 maps to new index 3, and so on
remappedHierClusterIndexMapColor = Color_Clust(HierIndexMappingColor);
HierIndexRenumber = [2 1 3 4];
NewHierClust_Optp = arrayfun(@(x) HierIndexRenumber(x), HCClust_Optp);
for i = length(Groups) :-1 :1
    xdata = Ytsne(HCClust_Optp == i,1);
    ydata = Ytsne(HCClust_Optp == i,2);
    len(i) = length(xdata);
    combinedX = [combinedX; xdata];
    combinedY = [combinedY; ydata];
    combinedColors = [combinedColors; remappedHierClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+len(4));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;

xrange = max(Ytsne(:,1))-min(Ytsne(:,1));
yrange = max(Ytsne(:,2))-min(Ytsne(:,2));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('Hierchical Clustering','FontSize',10);
xlim([min(Ytsne(:,1))-0.3*xrange  max(Ytsne(:,1))+0.3*xrange])
ylim([min(Ytsne(:,2))-0.15*yrange  max(Ytsne(:,2))+0.15*yrange])
xticks = -100:25:100;
yticks = -100:25:100;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;


UmapHCFig = figure(33);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);


Groups = {'Hierchical Cluster 1','Hierchical Cluster 2','Hierchical Cluster 3','Hierchical Cluster 4'};

combinedX = [ ];
combinedY = [ ];
combinedColors = [ ];

for i = length(Groups) :-1 :1
    xdata = Yumap(HCClust_Optp == i,1);
    ydata = Yumap(HCClust_Optp == i,2);
    len(i) = length(xdata);
    combinedX = [combinedX; xdata];
    combinedY = [combinedY; ydata];
    combinedColors = [combinedColors; remappedHierClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+len(4));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;

xrange = max(Yumap(:,1))-min(Yumap(:,1));
yrange = max(Yumap(:,2))-min(Yumap(:,2));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('Hierchical Clustering','FontSize',10);
xlim([min(Yumap(:,1))-0.3*xrange  max(Yumap(:,1))+0.3*xrange])
ylim([min(Yumap(:,2))-0.15*yrange  max(Yumap(:,2))+0.15*yrange])
xticks = -100:2:100;
yticks = -100:2:100;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;
%% GMM
% rng('default')
gmModel = fitgmdist(ANorm_Optp, 3,'RegularizationValue',1e-5,'Start', 'plus','Replicates',100);

GMMClust_Optp = cluster(gmModel, ANorm_Optp);

TsneGMMFig = figure(14);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...             o   %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);


Groups = {'GMM Cluster 1','GMM Cluster 2','GMM Cluster 3','GMM Cluster 4'};

combinedX = [ ];
combinedY = [ ];
combinedColors = [ ];
GMMIndexMappingColor = [2 4 1 3];  % For example, original index 2 maps to new index 3, and so on
remappedGMMClusterIndexMapColor = Color_Clust(GMMIndexMappingColor);
GMMIndexRenumber = [2 4 1 3];
NewGMMClust_Optp = arrayfun(@(x) GMMIndexRenumber(x), GMMClust_Optp);
for i = length(Groups) :-1 :1
    xdata = Ytsne(GMMClust_Optp == i,1);
    ydata = Ytsne(GMMClust_Optp == i,2);
    len(i) = length(xdata);
    combinedX = [combinedX; xdata];
    combinedY = [combinedY; ydata];
    combinedColors = [combinedColors; remappedGMMClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+len(4));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;

xrange = max(Ytsne(:,1))-min(Ytsne(:,1));
yrange = max(Ytsne(:,2))-min(Ytsne(:,2));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('Hierchical Clustering','FontSize',10);
xlim([min(Ytsne(:,1))-0.3*xrange  max(Ytsne(:,1))+0.3*xrange])
ylim([min(Ytsne(:,2))-0.15*yrange  max(Ytsne(:,2))+0.15*yrange])
xticks = -100:25:100;
yticks = -100:25:100;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;


UmapGMMFig = figure(44);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);


Groups = {'GMM Cluster 1','GMM Cluster 2','GMM Cluster 3','GMM Cluster 4'};

combinedX = [ ];
combinedY = [ ];
combinedColors = [ ];
NewGMMClust_Optp = arrayfun(@(x) GMMIndexRenumber(x), GMMClust_Optp);
for i = length(Groups) :-1 :1
    xdata = Yumap(GMMClust_Optp == i,1);
    ydata = Yumap(GMMClust_Optp == i,2);
    len(i) = length(xdata);
    combinedX = [combinedX; xdata];
    combinedY = [combinedY; ydata];
    combinedColors = [combinedColors; remappedGMMClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+len(4));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;

xrange = max(Yumap(:,1))-min(Yumap(:,1));
yrange = max(Yumap(:,2))-min(Yumap(:,2));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('Hierchical Clustering','FontSize',10);
xlim([min(Yumap(:,1))-0.3*xrange  max(Yumap(:,1))+0.3*xrange])
ylim([min(Yumap(:,2))-0.15*yrange  max(Yumap(:,2))+0.15*yrange])
xticks = -100:2:100;
yticks = -100:2:100;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;
%% Initialize final clusters
SumClusters = [NewKMnsClust_Optp NewHierClust_Optp];
FinalClusters = NaN(343,1);
% Recalculate riskLevel with new criteria
for i = 1:size(SumClusters, 1)
    if NewKMnsClust_Optp(i) == NewHierClust_Optp(i)
       FinalClusters(i) = NewKMnsClust_Optp(i);
    else
        FinalClusters(i) = 4;
    end
end

TsneFCFig = figure(8);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);



Groups = {'Cluster 1','Cluster 2','Cluster 3','Unconsistent'};
combinedX = [ ];
combinedY = [ ];
combinedColors = [ ];
IndexMappingColor = [1 2 3 4];  % For example, original index 2 maps to new index 3, and so on
% Replace original indices with new indices for K-means clustering
remappedClusterIndexMapColor = Color_Clust(IndexMappingColor);
IndexRenumber = [1 2 3 4];
NewClust_Optp = arrayfun(@(x) IndexRenumber(x), FinalClusters);

for i = length(Groups) :-1 :1
    xdata = Ytsne(FinalClusters == i,1);
    ydata = Ytsne(FinalClusters == i,2);
    len(i) = length(xdata);
    combinedX = [combinedX; xdata];
    combinedY = [combinedY; ydata];
    combinedColors = [combinedColors; remappedClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+ len(4) + len(5));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(Ytsne(:,1))-min(Ytsne(:,1));
yrange = max(Ytsne(:,2))-min(Ytsne(:,2));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('K-Means Clustering','FontSize',10);
xlim([min(Ytsne(:,1))-0.3*xrange  max(Ytsne(:,1))+0.3*xrange])
ylim([min(Ytsne(:,2))-0.15*yrange  max(Ytsne(:,2))+0.15*yrange])
xticks = -100:25:100;
yticks = -100:25:100;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;

UmapFCFig = figure(88);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);



combinedX = [ ];
combinedY = [ ];
combinedColors = [ ];

for i = length(Groups) :-1 :1
    xdata = Yumap(FinalClusters == i,1);
    ydata = Yumap(FinalClusters == i,2);
    len(i) = length(xdata);
    combinedX = [combinedX; xdata];
    combinedY = [combinedY; ydata];
    combinedColors = [combinedColors; remappedClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+ len(4) + len(5));
for i = 1:length(idx)
    scatter(combinedX(idx(i)), combinedY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(Yumap(:,1))-min(Yumap(:,1));
yrange = max(Yumap(:,2))-min(Yumap(:,2));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('K-Means Clustering','FontSize',10);
xlim([min(Yumap(:,1))-0.3*xrange  max(Yumap(:,1))+0.3*xrange])
ylim([min(Yumap(:,2))-0.15*yrange  max(Yumap(:,2))+0.15*yrange])
xticks = -100:2:100;
yticks = -100:2:100;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;
%% make violin for prism 

structname = ParamsT.Properties.VariableNames;
for i  = 3:length(structname)
    slot = ParamsT.(structname{i});
    DataC1 = slot(FinalClusters==1);
    DataC2 = slot(FinalClusters==2);
    DataC3 = slot(FinalClusters==3);
    H = max([length(DataC1) length(DataC2) length(DataC3)]);
    newT  = table(nan(H,1),nan(H,1),nan(H,1),...
        'VariableNames', {'C1', 'C2', 'C3'});
    newT.C1(1:length(DataC1)) = DataC1;
    newT.C2(1:length(DataC2)) = DataC2;
    newT.C3(1:length(DataC3)) = DataC3;
    ViolinData.(structname{i}) = newT;
end
%% calculate percentange of patients in each group
HierC1 = HFType_Optp( NewHierClust_Optp==1 );
HierC2 = HFType_Optp( NewHierClust_Optp==2 );
HierC3 = HFType_Optp( NewHierClust_Optp==3 );

C1portion = [length(HierC1(strcmp('HFrEF',HierC1))) length(HierC1(strcmp('HFpEF',HierC1)))];
C2portion = [length(HierC2(strcmp('HFrEF',HierC2))) length(HierC2(strcmp('HFpEF',HierC2)))];
C3portion = [length(HierC3(strcmp('HFrEF',HierC3))) length(HierC3(strcmp('HFpEF',HierC3)))];
THier = array2table([C1portion;C2portion;C3portion]);
THier.Properties.VariableNames = {'HFrEF','HFpEF'};
THier.Properties.RowNames = {'C1','C2','C3'};
% 对7个聚类进行处理
KmsC1 = HFType_Optp(NewKMnsClust_Optp == 1);
KmsC2 = HFType_Optp(NewKMnsClust_Optp == 2);
KmsC3 = HFType_Optp(NewKMnsClust_Optp == 3);
KmsC4 = HFType_Optp(NewKMnsClust_Optp == 4);
KmsC5 = HFType_Optp(NewKMnsClust_Optp == 5);
KmsC6 = HFType_Optp(NewKMnsClust_Optp == 6);
KmsC7 = HFType_Optp(NewKMnsClust_Optp == 7);

% 计算每个聚类中HFrEF和HFpEF的数量
C1portion = [length(KmsC1(strcmp('HFrEF', KmsC1))), length(KmsC1(strcmp('HFpEF', KmsC1)))];
C2portion = [length(KmsC2(strcmp('HFrEF', KmsC2))), length(KmsC2(strcmp('HFpEF', KmsC2)))];
C3portion = [length(KmsC3(strcmp('HFrEF', KmsC3))), length(KmsC3(strcmp('HFpEF', KmsC3)))];
C4portion = [length(KmsC4(strcmp('HFrEF', KmsC4))), length(KmsC4(strcmp('HFpEF', KmsC4)))];
C5portion = [length(KmsC5(strcmp('HFrEF', KmsC5))), length(KmsC5(strcmp('HFpEF', KmsC5)))];
C6portion = [length(KmsC6(strcmp('HFrEF', KmsC6))), length(KmsC6(strcmp('HFpEF', KmsC6)))];
C7portion = [length(KmsC7(strcmp('HFrEF', KmsC7))), length(KmsC7(strcmp('HFpEF', KmsC7)))];

% 将结果放入表格中
TKms = array2table([C1portion; C2portion; C3portion; C4portion; C5portion; C6portion; C7portion]);
TKms.Properties.VariableNames = {'HFrEF', 'HFpEF'};
TKms.Properties.RowNames = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};


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





