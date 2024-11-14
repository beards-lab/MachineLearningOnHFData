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



Color_Clust = {[0.4940 0.1840 0.5560],[0 0.75 0.75],[0.8500 0.3250 0.0980],[0.5 0.5 0.5],[0.4660 0.6740 0.1880]};


%% Load optimized parameter values
% Read the data and remove specified rows
ParamsT = readtable("Predictors.xlsx", 'VariableNamingRule', 'preserve');
shuntlist = [34 41 54 61 83 116 183 231 268 278 312];
ParamsT(shuntlist, :) = [];

% Keep columns up to the last one ending with _D
d_columns = find(endsWith(ParamsT.Properties.VariableNames, '_D'));
if ~isempty(d_columns)
    last_d_col = d_columns(end);
    ParamsT = ParamsT(:, 1:last_d_col);
end

% Convert categorical variables to binary

ParamsT.Race = strcmp(ParamsT.Race, 'Caucasian');
ParamsT.Race = double(ParamsT.Race);

ParamsT.Smoking = strcmp(ParamsT.Smoking, 'Current');
ParamsT.Smoking = double(ParamsT.Smoking);

% Convert Gender 2 to 0
ParamsT.Gender(ParamsT.Gender == 2) = 0;

ParamsT.Alcohol = strcmp(ParamsT.Alcohol, 'Yes');
ParamsT.Alcohol = double(ParamsT.Alcohol);

ParamsT.Drug = strcmp(ParamsT.Drug, 'Yes');
ParamsT.Drug = double(ParamsT.Drug);

% Convert cell arrays of strings to double
columnsToConvert = ParamsT.Properties.VariableNames;
for i = 5:length(columnsToConvert)
    colName = columnsToConvert{i};
    if iscell(ParamsT.(colName))
        ParamsT.(colName) = str2double(ParamsT.(colName));
    end
end

% check out std and delete the repeat ones
[~,locs] = unique(round(table2array(mean(ParamsT(:,5:end))),4));
locs = locs+4;
locs = sort(locs);
locs = [1;4;locs];
ParamsT = ParamsT(:,locs);

% Separate numeric and non-numeric columns
numericCols = varfun(@isnumeric, ParamsT, 'OutputFormat', 'uniform');
numericParamT = ParamsT(:, numericCols);
nonNumericParamT = ParamsT(:, ~numericCols);

% Remove columns with more than 25% missing values from numeric columns
missingThreshold = 0.25;
colsToRemove = varfun(@(x) mean(isnan(x)) > missingThreshold, numericParamT, 'OutputFormat', 'uniform');
numericParamT(:, colsToRemove) = [];

% Use built-in kNN imputation on numeric columns
numericParamTArray = table2array(numericParamT);
numericParamTArray = knnimpute(numericParamTArray);
numericParamT = array2table(numericParamTArray, 'VariableNames', numericParamT.Properties.VariableNames);

% Combine numeric and non-numeric columns back into one table
ParamsT = [nonNumericParamT, numericParamT];
if ismember('PatID', ParamsT.Properties.VariableNames)
    ParamsT.Properties.VariableNames{'PatID'} = 'Patient_NO';
else
    error('Table does not contain a column named ''patID''.');
end

%%
start_col = 3; 
[num_rows, num_cols] = size(ParamsT);
processedData = ParamsT;

for col = start_col:num_cols
    colData = ParamsT{:, col};
    if ~all(ismember(colData, [0, 1]))  
        processedData{:, col} = zscore(colData);
    else
        processedData{:, col} = colData; 
    end
end
supervised = 0;
if supervised == 1
    slot = supervisedslot;
else
    slot = (3:width(ParamsT));
end
 % choose which params you want to put into PCA, the ideal one should be mRMR.
% but couldn't be figured it out.

Optp_Names = ParamsT.Properties.VariableNames(slot) ;
HFType_Optp = ParamsT.HFtype;          % Heart failure type
PatNum_Optp = ParamsT.Patient_NO;          % Patient number

NumPats_Optp = size(processedData,1);                  % Number of patients

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
Num_Optp = size(processedData,2)-2;                      % Number of optimized params

ANorm_Optp = processedData{:,3:end};
%%
% 选择最佳簇数
rng('default');
options = statset('MaxIter', 200);
aic = zeros(1, 10);
for k = 1:10
    gmModel = fitgmdist(ANorm_Optp, k,'Options', options, 'RegularizationValue', 0.01);
    aic(k) = gmModel.AIC;
end

% 绘制 BIC 曲线
figure;
plot(1:10, aic, '-o');
xlabel('Number of Components');
ylabel('AIC');
title('AIC for Different Number of Components');
% 指定要拟合的高斯成分的数量（簇数）
numComponents = 3;

% 拟合高斯混合模型
gmModel = fitgmdist(ANorm_Optp, numComponents,'Options', options, 'RegularizationValue', 0.01);

% 获取每个数据点的聚类标签
clusterIdx = cluster(gmModel, ANorm_Optp);

% 可视化结果
figure;
gscatter(ANorm_Optp(:,20), ANorm_Optp(:,21), clusterIdx);
hold on;
% 可视化高斯成分的均值
plot(gmModel.mu(:,1), gmModel.mu(:,2), 'kx', 'LineWidth', 2, 'MarkerSize', 10);
legend('Cluster 1','Cluster 2','GMM Means');
title('Gaussian Mixture Model Clustering');
hold off;
%% Get clusters and FAMD from R 
resultTable = readtable('Final_Clustered_Data.xlsx',VariableNamingRule='preserve');
%% Visualized FAMD results
% Plotting out the score (A*V) for two PCs of the optimized params
ScrSize = get(0,'ScreenSize');                  % Getting screen size
AVOptp_2PCFig = figure(1);
clf
set(gcf,'Position', ...              % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);
PCxidx = 1;
PCyidx = 2;
% extract PC for subgroup
Groups = {'HFrEF','HFpEF'};



combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [];
for i = length(Groups) :-1 :1
    pcx = resultTable{strcmp(Groups{i},HFType_Optp),PCxidx};
    pcy = resultTable{strcmp(Groups{i},HFType_Optp),PCyidx};
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

xrange = max(resultTable{:,PCxidx})-min(resultTable{:,PCxidx});
yrange = max(resultTable{:,PCyidx})-min(resultTable{:,PCyidx});


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('PCA Analysis','FontSize',10);
xlim([min(resultTable{:,PCxidx})-0.3*xrange  max(resultTable{:,PCxidx})+0.3*xrange])
ylim([min(resultTable{:,PCyidx})-0.15*yrange  max(resultTable{:,PCyidx})+0.15*yrange])
xticks = -10:5:20;
yticks = -10:5:20;
% legend(Sc,{'HFrEF', 'HFmrEF','HFpEF'},'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30)
% text(max(AV_Optp(:,PCxidx))-0.38*xrange,max(AV_Optp(:,PCyidx))+0.15*yrange,VarExpStr_Optp,'FontSize',26, 'FontWeight','bold','Color',[0 0 0]);


%% Now use the 2PC plot to visualize the kmeans groups
% Creating the figure and plotting their positions
AVOptp_2PCKMnsFig = figure(12);
clf
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);
KProtoClust_Optp = clusterIdx;


Groups = {'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5'};
combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [ ];
KmeansIndexMappingColor = [3 2 1 4 5];  % For example, original index 2 maps to new index 3, and so on
% Replace original indices with new indices for K-means clustering
remappedKmeansClusterIndexMapColor = Color_Clust(KmeansIndexMappingColor);
KmsIndexRenumber = [3 1 2 4 5];
NewKMnsClust_Optp = arrayfun(@(x) KmsIndexRenumber(x), KProtoClust_Optp);

for i = length(Groups) :-1 :1
    xdata = resultTable{KProtoClust_Optp == i,PCxidx};
    ydata = resultTable{KProtoClust_Optp == i,PCyidx};
    len(i) = length(xdata);
    combinedPCX = [combinedPCX; xdata];
    combinedPCY = [combinedPCY; ydata];
    combinedColors = [combinedColors; remappedKmeansClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2) + len(3)+ len(4) +len(5));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end

xrange = max(resultTable{:,PCxidx})-min(resultTable{:,PCxidx});
yrange = max(resultTable{:,PCyidx})-min(resultTable{:,PCyidx});

xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('K-Means Clustering','FontSize',10);
xlim([min(resultTable{:,PCxidx})-0.3*xrange  max(resultTable{:,PCxidx})+0.3*xrange])
ylim([min(resultTable{:,PCyidx})-0.15*yrange  max(resultTable{:,PCyidx})+0.15*yrange])
xticks = -10:5:10;
yticks = -10:5:10;
% legend(Sc,Groups,'Location','southwest');
% legend(Box="off");
set(gca,'Xtick',xticks,'XTickLabel',[]);
set(gca,'Ytick',yticks,'YTickLabel',[]);
set(gca,'FontSize',30);box on;
hold off;



%% Now use the 2PC plot to visualize the HCC groups
% Creating the figure and plotting their positions
AVOptp_2PCHCFig = figure(13);
clf
HCClust_Optp = resultTable.HC_Cluster;
set(gcf,'Position', ...          % Positioning the figure
    [ScrSize(3)/30 ScrSize(4)/30 ...                %  on the screen
    ScrSize(3)/1.5 ScrSize(4)/1.5]);


Groups = {'Hierchical Cluster 1','Hierchical Cluster 2','Hierchical Cluster 3'};

combinedPCX = [ ];
combinedPCY = [ ];
combinedColors = [ ];
HierIndexMappingColor = [3 2 1 4];  % For example, original index 2 maps to new index 3, and so on
remappedHierClusterIndexMapColor = Color_Clust(HierIndexMappingColor);
HierIndexRenumber = [2 1 3 4];
NewHierClust_Optp = arrayfun(@(x) HierIndexRenumber(x), HCClust_Optp);
for i = length(Groups) :-1 :1
    xdata = resultTable{HCClust_Optp == i,PCxidx};
    ydata = resultTable{HCClust_Optp == i,PCyidx};
    len(i) = length(xdata);
    combinedPCX = [combinedPCX; xdata];
    combinedPCY = [combinedPCY; ydata];
    combinedColors = [combinedColors; remappedHierClusterIndexMapColor{i}.* ones(length(xdata), 1)];
end
idx = randperm(len(1) + len(2)+len(3));
for i = 1:length(idx)
    scatter(combinedPCX(idx(i)), combinedPCY(idx(i)), 'Marker','o','SizeData',200,'MarkerFaceColor',combinedColors(idx(i),:),'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.66);hold on;
end


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;

xrange = max(resultTable{:,PCxidx})-min(resultTable{:,PCxidx});
yrange = max(resultTable{:,PCyidx})-min(resultTable{:,PCyidx});


xline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
yline(0, '--', 'LineWidth', 1.5, 'Color', 'k'); hold on;
% xlabel('Principal Component 1', ...                 % Labelling axes with
%     'FontSize',30,'FontWeight','bold')              %  PC labels
% ylabel('Principal Component 2', ...
%     'FontSize',30,'FontWeight','bold')

% title('Hierchical Clustering','FontSize',10);
xlim([min(resultTable{:,PCxidx})-0.3*xrange  max(resultTable{:,PCxidx})+0.3*xrange])
ylim([min(resultTable{:,PCyidx})-0.15*yrange  max(resultTable{:,PCyidx})+0.15*yrange])
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

% K-means Clustering
KmsC1 = HFType_Optp(NewKMnsClust_Optp == 1);
KmsC2 = HFType_Optp(NewKMnsClust_Optp == 2);
KmsC3 = HFType_Optp(NewKMnsClust_Optp == 3);
KmsC4 = HFType_Optp(NewKMnsClust_Optp == 4);
KmsC5 = HFType_Optp(NewKMnsClust_Optp == 5);

C1portion = [length(KmsC1(strcmp('HFrEF',KmsC1))) length(KmsC1(strcmp('HFpEF',KmsC1)))];
C2portion = [length(KmsC2(strcmp('HFrEF',KmsC2))) length(KmsC2(strcmp('HFpEF',KmsC2)))];
C3portion = [length(KmsC3(strcmp('HFrEF',KmsC3))) length(KmsC3(strcmp('HFpEF',KmsC3)))];
C4portion = [length(KmsC4(strcmp('HFrEF',KmsC4))) length(KmsC4(strcmp('HFpEF',KmsC4)))];
C5portion = [length(KmsC5(strcmp('HFrEF',KmsC5))) length(KmsC5(strcmp('HFpEF',KmsC5)))];

TKms = array2table([C1portion; C2portion; C3portion; C4portion;C5portion]);
TKms.Properties.VariableNames = {'HFrEF', 'HFpEF'};
TKms.Properties.RowNames = {'C1', 'C2', 'C3', 'C4','C5'};


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





