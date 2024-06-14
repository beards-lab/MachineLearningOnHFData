clear
data = readtable('coxregression.xlsx',VariableNamingRule='preserve');
% % 排除不应该标准化的列索引
% cols_to_exclude = [16, 68, 70, 71, 72,73,74];
%
% % 获取所有列的索引
% all_cols = 1:width(data);
%
% % 除去不需要标准化的列，获得需要标准化的列的索引
% cols_to_standardize = setdiff(all_cols, cols_to_exclude);
%
% % 标准化数据：对标准化列索引中的每一列执行zscore
% for i = cols_to_standardize
%     % 仅对数值数据进行标准化
%     if isnumeric(data{:, i})
%         data{:, i} = zscore(data{:, i});
%     end
% end
T = data.Time;
C = data.Event;
C = 1 - C;
X0 = data{:,{'LVEF(MRI)','LVEF（echo）','Aging','SBP','Sex','BMI','Current smoker','Diabetes','Diagnosis of COPD'}};
HR0 = NaN(width(X0),1);
CIlower0 =  NaN(width(X0),1);
CIupper0 =  NaN(width(X0),1);
% Get HR for X0
for i = 1:width(X0)
    [b,logl,H,stats] = coxphfit(X0(:,i),T,'Censoring',C,'Baseline',0);
    HR0(i) = exp(b);
    CIlower0(i) = exp(b - 1.96 * stats.se);
    CIupper0(i) = exp(b + 1.96 * stats.se);
end
HR0T = table(HR0,CIlower0,CIupper0);
HR0T.Properties.RowNames = {'LVEF(MRI)','LVEF（echo）','Aging','SBP','Sex','BMI','Current smoker','Diabetes','Diagnosis of COPD'};
%%
% Initialize the result table
varNames = data.Properties.VariableNames(1:63);  % Assuming these are the names of the first 63 columns
resultsTable = table('Size', [63 3], 'VariableTypes', {'double', 'double', 'double'}, 'VariableNames', {'HR', 'CI_Lower', 'CI_Upper'});

% Loop through the first 63 columns
for i = 1:63
    % Extract the current column and add it to the existing covariates
    currentX = data{:, i};
    X = [ currentX];

    % Run the Cox regression
    [b, logl, H, stats] = coxphfit(X, T, 'Censoring', C, 'Baseline', 0);

    % Extract HR and CI for the new variable
    newVariableIdx = size(X, 2);  % index of the newly added variable in the regression model
    HR = exp(b(newVariableIdx));  % hazard ratio
    CI = exp(stats.beta(newVariableIdx) + [-1.96, 1.96] .* stats.se(newVariableIdx));  % 95% CI

    % Store the results into the table
    resultsTable.Properties.RowNames(i) = varNames(i);
    resultsTable.HR(i) = HR;
    resultsTable.CI_Lower(i) = CI(1);
    resultsTable.CI_Upper(i) = CI(2);
end

%% Create a forest plot
figure(1); clf;set(gcf, 'Position', [800, 50, 600, 1000])
for i = 1:height(HR0T)
    line([HR0T.CIlower0(i) HR0T.CIupper0(i)], [i i], 'Color', 'k'); % 95% CI line
    hold on;
    % Check if CI is entirely above 1 (significant positive relationship)
    if HR0T.CIlower0(i) > 1
        plotColor = 'r'; % Red for HR significantly > 1
        % Check if CI is entirely below 1 (significant negative relationship)
    elseif HR0T.CIupper0(i) < 1
        plotColor = 'b'; % Blue for HR significantly < 1
    else
        plotColor = 'k'; % Black for non-significant HR
    end

    % Plot the point with the determined color
    plot(HR0T.HR0(i), i, 'o', 'MarkerFaceColor', plotColor, 'MarkerEdgeColor', plotColor);
end

% Set y-axis labels and title
set(gca, 'YTick', 1:height(HR0T), 'YTickLabel', HR0T.Properties.RowNames, 'FontSize', 10);
xlabel('Hazard Ratio', 'FontSize', 14);

% Set logarithmic scale
set(gca, 'XScale', 'log');
% Set the X-axis ticks and tick labels
set(gca, 'XTick', [0.01 0.1 1 10 100]);
set(gca, 'XTickLabel', {'0.01', '0.1', '1', '10', '100'});
xlim([0.2 5]);

% Add a vertical line to represent HR=1
line([1, 1], [0, height(HR0T) + 2], 'Color', 'b', 'LineStyle', '--');
ylim([-1 11]);
set(gca, 'TickDir', 'out', 'YDir', 'reverse');
hold off;
figure(2);clf;set(gcf, 'Position', [200, 50, 600, 1000])
for i = 1:height(resultsTable)
    hold on;
    % Check if CI is entirely above 1 (significant positive relationship)
    if resultsTable.CI_Lower(i) > 1
        plotColor = 'r'; % Red for HR significantly > 1
        % Check if CI is entirely below 1 (significant negative relationship)
    elseif resultsTable.CI_Upper(i) < 1
        plotColor = 'b'; % Blue for HR significantly < 1
    else
        plotColor = 'k'; % Black for non-significant HR
    end

    % Plot the point with the determined color
    if i <= 21
        line([resultsTable.CI_Lower(i) resultsTable.CI_Upper(i)], [i-4 i-4], 'Color', 'k'); % 95% CI line
        plot(resultsTable.HR(i), i-4, 'o', 'MarkerFaceColor', plotColor, 'MarkerEdgeColor', plotColor);
    elseif i<= 35
        line([resultsTable.CI_Lower(i) resultsTable.CI_Upper(i)], [i-2 i-2], 'Color', 'k'); % 95% CI line
        plot(resultsTable.HR(i), i-2, 'o', 'MarkerFaceColor', plotColor, 'MarkerEdgeColor', plotColor);
    else
        line([resultsTable.CI_Lower(i) resultsTable.CI_Upper(i)], [i i], 'Color', 'k'); % 95% CI line
        plot(resultsTable.HR(i), i, 'o', 'MarkerFaceColor', plotColor, 'MarkerEdgeColor', plotColor);
    end

end

% Set y-axis labels and title
set(gca, 'YTick', [], 'YTickLabel', {''}, 'FontSize', 10);
% xlabel('Hazard Ratio', 'FontSize', 14);

% Set logarithmic scale
set(gca, 'XScale', 'log');
% Set the X-axis ticks and tick labels
set(gca, 'XTick', [ ]);
set(gca, 'XTickLabel', {''});

% Add a vertical line to represent HR=1
line([1, 1], [-5, height(resultsTable) + 7], 'Color', 'b', 'LineStyle', '--');
ylim([-5 65]);
xlim([0.2 5]);
% Reverse y-axis direction
set(gca, 'TickDir', 'out', 'YDir', 'reverse');
hold off;
%% Get fixed time ROC
% Define data columns for D, M, and DM
D = data{:,1:21};
M = data{:,[(22:35) (58:61) 63]};
DM = data{:,[16,18,19 (22:63)]};
E = ~C; % If C indicates censoring, then an event occurrence is the logical negation of C; if C indicates event occurrence, use C directly
timePoints = [1, 2, 3, 5, 15]; % Define the time points for ROC
T_year = cell(length(timePoints), 1);
E_year = cell(length(timePoints), 1);
for i = 1:length(timePoints)
    t = timePoints(i);
    T_year{i} = T;
    % Update event indicators for each time point, setting those still alive past the time point as censored
    T_year{i}(T > t) = t;
    E_year{i} = E;
    % Update event indicators for each time point, setting those still alive past the time point as censored
    E_year{i}(T > t) = 0;
end
n = 3; % Index for choosing the 3-year time point
% Fit Cox proportional hazards models
[bD, loglD, HD, statsD] = coxphfit(D, T_year{n}, 'Censoring', ~E_year{n}, 'Baseline', 0);
[bM, loglM, HM, statsM] = coxphfit(M, T_year{n}, 'Censoring', ~E_year{n}, 'Baseline', 0);
[bDM, loglDM, HDM, statsDM] = coxphfit(DM, T_year{n}, 'Censoring', ~E_year{n}, 'Baseline', 0);
% Calculate risk scores
riskScoresD = D * bD;
riskScoresM = M * bM;
riskScoresDM = DM * bDM;

% Plot ROC curves at each time point
for i = 1:length(timePoints)
    figure(i+2); clf; hold on; % Create a new figure and hold it for multiple plots
    % Calculate ROC and AUC for Model D
    [XrocD, YrocD, ~, AUCD] = perfcurve(E_year{i}, riskScoresD, 1);
    plot(XrocD, YrocD, 'b-', 'LineWidth', 2);
    % Calculate ROC and AUC for Model M
    [XrocM, YrocM, ~, AUCM] = perfcurve(E_year{i}, riskScoresM, 1);
    plot(XrocM, YrocM, 'r-', 'LineWidth', 2);
    % Calculate ROC and AUC for Model DM
    [XrocDM, YrocDM, ~, AUCDM] = perfcurve(E_year{i}, riskScoresDM, 1);
    plot(XrocDM, YrocDM, 'g-', 'LineWidth', 2);
    
    % Add legend and format chart
    c_indexD = compute_cindex(T_year{i}, ~E, riskScoresD);
    c_indexM = compute_cindex(T_year{i}, ~E, riskScoresM);
    c_indexDM = compute_cindex(T_year{i}, ~E, riskScoresDM);
    legend(sprintf('Data (AUC = %.3f)', AUCD), ...
        sprintf('Model (AUC = %.3f)', AUCM), ...
        sprintf('Data and Model (AUC = %.3f)', AUCDM), 'Location', 'southeast','Box', 'off', 'fontsize', 12);
    % legend(sprintf('Data (C-index = %.3f)', c_indexD), ...
    %        sprintf('Model (C-index = %.3f)', c_indexM), ...
    %        sprintf('Data and Model (C-index = %.3f)', c_indexDM), 'Location', 'southeast', 'Box', 'off', 'fontsize', 12);
    % 
    % Hide tick labels from x-axis and y-axis
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
    grid off; % Hide grid
    box on; % Show box
    pbaspect([1,1,1]); % Aspect ratio 1:1:1 for the plot
    hold off; % release figure hold
end
%% Get time-dependent ROC
% Create figures for plotting risk scores at each time point
for i = 1:length(timePoints)
    % Fit Cox proportional hazards model for each data set at current time point
    [bD, loglD, HD, statsD] = coxphfit(D, T_year{i}, 'Censoring', ~E_year{i}, 'Baseline', 0);
    [bM, loglM, HM, statsM] = coxphfit(M, T_year{i}, 'Censoring', ~E_year{i}, 'Baseline', 0);
    [bDM, loglDM, HDM, statsDM] = coxphfit(DM, T_year{i}, 'Censoring', ~E_year{i}, 'Baseline', 0);
    
    % Calculate risk scores
    riskScoresD = D * bD;
    riskScoresM = M * bM;
    riskScoresDM = DM * bDM;
    
    % % Calculate C-index for each model
    % c_indexD = compute_cindex(T_year{i}, ~E_year{i}, riskScoresD);
    % c_indexM = compute_cindex(T_year{i}, ~E_year{i}, riskScoresM);
    % c_indexDM = compute_cindex(T_year{i}, ~E_year{i}, riskScoresDM);
    
    % Create a new figure and hold for multiple ROC plots
    figure(i+10); hold on; 
    
    % Calculate and plot ROC curve and AUC for Model D
    [XrocD, YrocD, ~, AUCD] = perfcurve(E_year{i}, riskScoresD, 1);
    plot(XrocD, YrocD, 'b-', 'LineWidth', 2);

    % Calculate and plot ROC curve and AUC for Model M
    [XrocM, YrocM, ~, AUCM] = perfcurve(E_year{i}, riskScoresM, 1);
    plot(XrocM, YrocM, 'r-', 'LineWidth', 2);

    % Calculate and plot ROC curve and AUC for Model DM
    [XrocDM, YrocDM, ~, AUCDM] = perfcurve(E_year{i}, riskScoresDM, 1);
    plot(XrocDM, YrocDM, 'g-', 'LineWidth', 2);

    % Add legend and format the chart
     legend(sprintf('Data (AUC = %.3f)', AUCD), ...
        sprintf('Model (AUC = %.3f)', AUCM), ...
        sprintf('Data and Model (AUC = %.3f)', AUCDM), 'Location', 'southeast','box','off','fontsize',12);  % Hide x-axis and y-axis tick labels
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
    
    % Set aspect ratio and box
    pbaspect([1,1,1]);
    box on;
    
    % Turn off grid and release the figure hold
    grid off; hold off;
end