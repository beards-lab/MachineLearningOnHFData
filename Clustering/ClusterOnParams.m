%% Script Summary:
% This script applies a well-established clustering method to obtain clustering results and generate plots.
% The input is the "Features.xlsx" file, containing 25 functional parameters and 4 simulated regurgitation fractions.
% Outputs include clustering labels, which may require adjustment of the optimal cluster number and possibly changing the order in "Color_Clust" for better visualization.
% The results are also saved in the last two columns of "Features.xlsx".

% Created by Feng Gu
% Last modified: 10/29/2024
clear
%% Option flags for running script
% Colors for HFrEF and HFpEF
Color_HFType = {[0.84,0.08,0.18], ...  % orange for HFrEF
    [0,0.45,0.74]};   % blue for HFpEF
Color_Clust = {[0.4940 0.1840 0.5560], [0 0.75 0.75], [0.8500 0.3250 0.0980], [0.5 0.5 0.5], [0.4660 0.6740 0.1880]};

Color_Clust{end+1} = [1.0 0.84 0.0]; 
Color_Clust{end+1} = [0.0 1.0 0.0];  

%% Load optimized parameter values
% Load optimized parameters
ParamsT = readtable("Features.xlsx",'VariableNamingRule','preserve'); % Digital twins, replace valve backward resistant to regurigation fraction.
ParamsT.MVr_S = (ParamsT.MVr_S-1).*20;% Linear interpretation to assign each value the same magnitude
ParamsT.AVr_S = (ParamsT.AVr_S-1).*20;
ParamsT.PVr_S = (ParamsT.PVr_S-1).*15;
ParamsT.TVr_S = (ParamsT.TVr_S-1).*17.5;
ParamsT = ParamsT(:,1:width(ParamsT)-2);% The last two columns represent the clustering results, where a value of 4 indicates inconsistent clustering.
% Features Preprocessing
% Delete the replaced column if it exists.
[~,locs] = unique(round(table2array(mean(ParamsT(:,3:end))),4));
locs = locs+2;
locs = sort(locs);
locs = [1;2;locs];
ParamsT = ParamsT(:,locs);
% Kill out of bound values
data_normalized = zscore(ParamsT{:,3:end});
z_threshold = 3;
outliers1 = data_normalized > z_threshold;
outliers2 = data_normalized < -z_threshold;
[~,FeatureNumber] = size(data_normalized);
for i = 1:FeatureNumber
    ParamsT{:,i+2}(outliers1(:,i)) = max(ParamsT{:,i+2}(~outliers1(:,i)));
    ParamsT{:,i+2}(outliers2(:,i)) = min(ParamsT{:,i+2}(~outliers2(:,i)));
end
Names = ParamsT.Properties.VariableNames;

%% Preparing PCA
supervised = 0;
if supervised == 1
    slot = supervisedslot;
else
    slot = (3:width(ParamsT));
end

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
    norm([sigma_Optp(PCxidx) sigma_Optp(PCyidx)])^2 / TotVar_Optp;          %  by first 1 and 2 PCs
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


xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));


VarExpStr_Optp = ['Variance explained ', ...
    num2str(round(VarExp2PC_Optp,2)*100),'%'];

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
xlabel('Principal Component 1', ...                 % Labelling axes with
    'FontSize',30,'FontWeight','bold')              %  PC labels
ylabel('Principal Component 2', ...
    'FontSize',30,'FontWeight','bold')

title('PCA Analysis','FontSize',20);
xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
xticks = -10:5:20;
yticks = -10:5:20;
legend({'HFpEF', 'HFrEF'},'Location','southwest');
legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',20)
text(max(AV_Optp(:,PCxidx))-0.38*xrange,max(AV_Optp(:,PCyidx))+0.05*yrange,VarExpStr_Optp,'FontSize',26, 'FontWeight','bold','Color',[0 0 0]);

%% Evaluate Optimal Number of Clusters for Kmeans and Hierarchical Clustering
% % find a reasonable cluster number
% rng("default") % Set the random number generator to default for reproducibility
% maxClusters = 25;  % Maximum number of clusters to evaluate
% 
% % Calculate silhouette scores for Kmeans
% silhouetteEvaluation = evalclusters(ANorm_Optp, "kmeans", "silhouette", "KList", 1:maxClusters);
% % Calculate Davies-Bouldin Index for Kmeans
% dbiEvaluation = evalclusters(ANorm_Optp, "kmeans", "DaviesBouldin", "KList", 1:maxClusters);
% % Calculate GAP Index for Kmeans
% GAPEvaluation = evalclusters(ANorm_Optp, "kmeans", "gap", "KList", 1:maxClusters);
% 
% % Hierarchical Clustering Evaluations
% % Calculate silhouette scores for Hierarchical Clustering
% silhouetteHierEvaluation = evalclusters(ANorm_Optp, "linkage", "silhouette", "KList", 1:maxClusters);
% % Calculate Davies-Bouldin Index for Hierarchical Clustering
% dbiHierEvaluation = evalclusters(ANorm_Optp, "linkage", "DaviesBouldin", "KList", 1:maxClusters);
% % Calculate GAP Index for Hierarchical Clustering
% GAPHierEvaluation = evalclusters(ANorm_Optp, "linkage", "gap", "KList", 1:maxClusters);
% 
% % GMM Clustering Evaluations
% % Calculate silhouette scores for Hierarchical Clustering
% silhouetteGMMEvaluation = evalclusters(ANorm_Optp, "gmdistribution", "silhouette", "KList", 1:maxClusters);
% % Calculate Davies-Bouldin Index for Hierarchical Clustering
% dbiGMMEvaluation = evalclusters(ANorm_Optp, "gmdistribution", "DaviesBouldin", "KList", 1:maxClusters);
% % Calculate GAP Index for Hierarchical Clustering
% GAPGMMEvaluation = evalclusters(ANorm_Optp, "gmdistribution", "gap", "KList", 1:maxClusters);

%% Plot the evaluation metrics to assist with determining the optimal number of clusters
% figure(102);
% 
% range = [0 25];
% optimalNumClusters1 = 3;
% optimalNumClusters2 = 4;
% fontSize = 14;  
% 
% % Plot the Davies-Bouldin Index for Kmeans
% subplot(3,2,1);
% plot(1:maxClusters, dbiEvaluation.CriterionValues, 'b-o');
% title('Kmeans Davies-Bouldin Index', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('DBI Score', 'FontSize', fontSize);
% xlim(range);
% hold on; 
% scatter(optimalNumClusters1, dbiEvaluation.CriterionValues(optimalNumClusters1), 'r', 'filled');
% hold off
% 
% % Plot the GAP Index for Kmeans
% subplot(3,2,2);
% plot(1:maxClusters, GAPEvaluation.CriterionValues, 'b-o');
% title('Kmeans GAP Index', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('GAP Score', 'FontSize', fontSize);
% xlim(range);
% hold on; 
% scatter(optimalNumClusters1, GAPEvaluation.CriterionValues(optimalNumClusters1), 'r', 'filled');
% hold off
% 
% % Plot the Davies-Bouldin Index for Hierarchical Clustering
% subplot(3,2,3);
% plot(1:maxClusters, dbiHierEvaluation.CriterionValues, 'b-o');
% title('Hierarchical Davies-Bouldin Index', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('DBI Score', 'FontSize', fontSize);
% xlim(range);
% hold on; 
% scatter(optimalNumClusters2, dbiHierEvaluation.CriterionValues(optimalNumClusters2), 'r', 'filled');
% hold off
% 
% % Plot the GAP Index for Hierarchical Clustering
% subplot(3,2,4);
% plot(1:maxClusters, GAPHierEvaluation.CriterionValues, 'b-o');
% title('Hierarchical GAPIndex', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('GAP Score', 'FontSize', fontSize);
% xlim(range);
% hold on; 
% scatter(optimalNumClusters1, GAPHierEvaluation.CriterionValues(optimalNumClusters1), 'r', 'filled');
% hold off
% 
% % GMM is soft clustering and reflect proability, every time running results will be different. 
% % Plot the Davies-Bouldin Index for GMM Clustering
% subplot(3,2,5);
% plot(1:maxClusters, dbiGMMEvaluation.CriterionValues, 'b-o');
% title('GMM Davies-Bouldin Index', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('DBI Score', 'FontSize', fontSize);
% xlim(range);
% hold on; 
% scatter(optimalNumClusters1, dbiGMMEvaluation.CriterionValues(optimalNumClusters1), 'r', 'filled');
% hold off
% 
% % Plot the GAP Index for Hierarchical Clustering
% subplot(3,2,6);
% plot(1:maxClusters, GAPGMMEvaluation.CriterionValues, 'b-o');
% title('GMM GAPIndex', 'FontSize', fontSize);
% xlabel('Number of Clusters', 'FontSize', fontSize);
% ylabel('GAP Score', 'FontSize', fontSize);
% xlim(range);
% hold on; 
% scatter(optimalNumClusters2, GAPGMMEvaluation.CriterionValues(optimalNumClusters2), 'r', 'filled');
% hold off
% 
% 
% % Adjust the layout to prevent labels from overlapping
% set(gcf, 'Position', [100, 100, 1024, 768]);  % Resize figure to make it wider

%% Grouping into clusters using kmeans on optimized parameters

rng('default')
[KMnsClust_Optp,KMnsCntrd_Optp, ...                 % KMeans clust numbers
    KMnsSumDst_Optp] = ...                          %  centroids, sum of dist
    kmeans(ANorm_Optp,3, ...            %  Norm data, # of clusts, could be eithor 3 or 4
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

AVOptp_2PCKMnsFig = figure(12);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);

Groups = {'Kmeans Cluster 1','Kmeans Cluster 2','Kmeans Cluster 3','Kmeans Cluster 4','Kmeans Cluster 5','Kmeans Cluster 6','Kmeans Cluster 7'};
combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [ ];
KmeansIndexMappingColor = [1 2 3 4 5 6 7];  % For example, original index 2 maps to new index 3, and so on
% Replace original indices with new indices for K-means clustering
remappedKmeansClusterIndexMapColor = Color_Clust(KmeansIndexMappingColor);
KmsIndexRenumber = [1 2 3 4 5 6 7];
NewKMnsClust_Optp = arrayfun(@(x) KmsIndexRenumber(x), KMnsClust_Optp);

for i = length(Groups) :-1 :1
    xdata = AV_Optp(KMnsClust_Optp == i,PCxidx);
    ydata = AV_Optp(KMnsClust_Optp == i,PCyidx);
    len(i) = length(xdata);
    combinedPCX = [combinedPCX; xdata];
    combinedPCY = [combinedPCY; ydata];
    combinedColors = [combinedColors; remappedKmeansClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+ len(4) + len(5) + len(6)+ len(7));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
xlabel('Principal Component 1', ...                 % Labelling axes with
    'FontSize',30,'FontWeight','bold')              %  PC labels
ylabel('Principal Component 2', ...
    'FontSize',30,'FontWeight','bold')

title('K-Means Clustering','FontSize',10);
xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
xticks = -10:5:10;
yticks = -10:5:10;
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
    'Maxclust',3);

% Second or third-from-last linkages.
HCCutOff_Optp = median([HCLink_Optp(end-2,3) ...
    HCLink_Optp(end-1,3)]);
% % assign the color I want
% H = dendrogram(HCLink_Optp, 0, 'ColorThreshold', HCCutOff_Optp, 'Labels', num2str(PatNum_Optp));
% for i = 1:length(H)
%     x = get(H(i), 'XData');
%     THier = sort(HCClust_Optp);
%     group = THier(round(x(1))); 
%     set(H(i), 'Color', Color_HCCClust{group}); 
% end

cgo=clustergram (ANorm_Optp,'Linkage','ward','Cluster',1, 'Dendrogram',HCCutOff_Optp,'Colormap',redbluecmap);
% set(cgo,'RowLabels',PatNum_Optp,'ColumnLabels',{'C_S_A','C_S_V','C_P_A','C_P_V','kact_L_V','kact_R_V','kpas_L_V','kpas_R_V',...
%      'Vw_L_V','Vw_R_V','Vw_S_E_P','V_a_u','V_a_c','V_4_C','K_P_C','B_P_C','V_V_S','R_S_A','R_P_A','R_t_S_A','R_t_P_A','R_A_t','R_t_o','R_t_c','R_p_o','R_p_c','R_m_o','R_m_c','R_a_o','R_a_c',},'ColumnLabelsRotate', 45)


% I couldn't figure out how to improve the appearance of the dendrogram; it may require manual adjustments using the MATLAB GUI.
set(cgo,'RowLabels',PatNum_Optp,'ColumnLabels',Names(3:end),'ColumnLabelsRotate', 45)


%%
% Now use the 2PC plot to visualize the HCC groups
% Creating the figure and plotting their positions
AVOptp_2PCHCFig = figure(13);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);


Groups = {'Hierchical Cluster 1','Hierchical Cluster 2','Hierchical Cluster 3'};

combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [ ];
HierIndexMappingColor = [2 1 3 4];  % For example, original index 2 maps to new index 3, and so on
remappedHierClusterIndexMapColor = Color_Clust(HierIndexMappingColor);
HierIndexRenumber = [2 1 3 4];
NewHierClust_Optp = arrayfun(@(x) HierIndexRenumber(x), HCClust_Optp);
for i = length(Groups) :-1 :1
    xdata = AV_Optp(HCClust_Optp == i,PCxidx);
    ydata = AV_Optp(HCClust_Optp == i,PCyidx);
    len(i) = length(xdata);
    combinedPCX = [combinedPCX; xdata];
    combinedPCY = [combinedPCY; ydata];
    combinedColors = [combinedColors; remappedHierClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;

xrange = max(AV_Optp(:,PCxidx))-min(AV_Optp(:,PCxidx));
yrange = max(AV_Optp(:,PCyidx))-min(AV_Optp(:,PCyidx));

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
xlabel('Principal Component 1', ...                 % Labelling axes with
    'FontSize',30,'FontWeight','bold')              %  PC labels
ylabel('Principal Component 2', ...
    'FontSize',30,'FontWeight','bold')

title('Hierchical Clustering','FontSize',30);
xlim([min(AV_Optp(:,PCxidx))-0.3*xrange  max(AV_Optp(:,PCxidx))+0.3*xrange])
ylim([min(AV_Optp(:,PCyidx))-0.15*yrange  max(AV_Optp(:,PCyidx))+0.15*yrange])
xticks = -10:5:10;
yticks = -10:5:10;

set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;

%% Final output
% There is a trick here: once the optimal cluster number is fixed, the
% cluster number will remain consistent for hierarchical clustering.
% However, for k-means, the cluster numbers can change each time. This
% means that while the clustering results for k-means are exactly the same,
% the cluster labels (e.g., cluster 1) could correspond to different groups
% in different runs (e.g., cluster 2 in another run). After running
% multiple times, there is a chance that both clustering methods will
% assign the same group the same number. Once this happens, you can copy
% the cluster results from the environment (variable KMnsClust_Optp and
% HCClust_Optp) and paste them into the "Features.xlsx" file for the next
% script. The "risk" column corresponds to cluster number 3, and the
% "risk2" column corresponds to cluster number 4. Cluster number 4
% represents individuals who are not consistently clustered by k-means and
% hierarchical clustering.

%% Function used
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





