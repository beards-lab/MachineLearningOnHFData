clear
Stable = readtable('PredictorsFinal.xlsx','VariableNamingRule','preserve');
Stable.HTTime = datetime(Stable.HTTime);
% Stable = sortrows(Stable,"ModelTime_End","ascend");
Endtime = min(Stable(:,{'DeathTime','LVADTime','HTTime'}),[],2);
figure(1);
clf;hold on;
ms = 64;
for patId = 1:height(Stable)
    Win = plot([Stable.ModelTime_Start(patId),Stable.ModelTime_End(patId)], [patId, patId], LineWidth=2, Color='k');
    if isnat(Endtime.min(patId))
        plot([Stable.ModelTime_End(patId),Stable.LastTimeToVisit(patId)], [patId, patId], LineWidth=0.1, Color=[0.5 0.5 0.5])
        Censor = scatter(Stable.LastTimeToVisit(patId),patId,ms,"blue",'filled','o');
    else
        plot([Stable.ModelTime_End(patId),Endtime.min(patId)], [patId, patId], LineWidth=0.1, Color=[0.5 0.5 0.5])
        EndEvent = scatter(Endtime.min(patId),patId,ms,"red",'filled','o');
    end
end
xrefline = datetime('2024-07-31');
xline(xrefline, '--k', 'LineWidth', 2);
% for legend only
xlim([datetime('2009-03-01') datetime('2024-10-31')]);
ylim([-9 356]);
set(gca, 'YTick', []);
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
pbaspect([1.67 1 1]);

%% right sensored 5 year survival
% Define the 5-year threshold
FiveYearThreshold = years(5);

% Prepare the figure
figure(2);
clf; hold on;
ms = 64;
Duration = NaN(343,1);
EndEvent = NaN(343,1);
for patId = 1:height(Stable)
    % Calculate the duration from ModelTime_End to the event or last visit
    if isnat(Endtime.min(patId))
        DurationToEvent = Stable.LastTimeToVisit(patId) - Stable.ModelTime_End(patId);
        Color = 'blue'; % Censoring beyond 5 years
        EndEvent(patId) = 0;
        if DurationToEvent > FiveYearThreshold
            DurationToEvent = FiveYearThreshold;
        end
    else
        DurationToEvent = Endtime.min(patId) - Stable.ModelTime_End(patId);
        if DurationToEvent > FiveYearThreshold
            DurationToEvent = FiveYearThreshold;
            Color = 'blue'; % Censoring beyond 5 years
            EndEvent(patId) = 0;
        else
            Color = 'red'; % Event within 5 years
            EndEvent(patId) = 1;
        end
    end
    Duration(patId) = years(DurationToEvent);

    % Truncate at 5 years if the duration exceeds it

    % Plot the timeline starting from ModelTime_End
    plot([0, DurationToEvent], [patId, patId], 'LineWidth', 0.1,  Color=[0.5 0.5 0.5]);

    % Plot the event marker
    scatter(DurationToEvent, patId, ms, Color, 'filled', 'o');
end

% Adjust axes limits and appearance
xlim([0 FiveYearThreshold + years(1)]); % Add some space beyond 5 years
ylim([-9 356]);
set(gca, 'YTick', []);
set(gca, 'XTick', [0 years(1) years(2) years(3) years(4) years(5) years(6)]);
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
pbaspect([1.67 1 1]);
%% Compute Statistics

% Create a logical index for HFpEF and other groups
isHFpEF = strcmp(Stable.HFtype, 'HFpEF');
isOther = ~isHFpEF;

% Numerical variables to analyze
numVars = Stable.Properties.VariableNames(82:121);

% Categorical variables to analyze
catVars = {'First_diagnosis_of_heart_failure_in_the_past_18_months'};

% Initialize structures to hold the results
results.HFpEF = struct();
results.Other = struct();

% Loop through numerical variables
for i = 1:length(numVars)
    varName = numVars{i};
    
    % HFpEF group
    dataHFpEF = Stable{isHFpEF, varName};
    nHFpEF = sum(~isnan(dataHFpEF)); % sample size
    results.HFpEF.(varName).mean = mean(Stable{isHFpEF, varName}, 'omitnan');
    results.HFpEF.(varName).std = std(Stable{isHFpEF, varName}, 'omitnan')/sqrt(nHFpEF);
    
    % Other group
    dataOther = Stable{isOther, varName};
    nOther = sum(~isnan(dataOther)); % sample size
    results.Other.(varName).mean = mean(Stable{isOther, varName}, 'omitnan');
    results.Other.(varName).std = std(Stable{isOther, varName}, 'omitnan')/sqrt(nOther);
end

% Loop through categorical variables
for i = 1:length(catVars)
    varName = catVars{i};
    
    % HFpEF group
    [uniqueValsHFpEF, ~, idxHFpEF] = unique(Stable{isHFpEF, varName});
    countsHFpEF = accumarray(idxHFpEF, 1);
    percentagesHFpEF = 100 * countsHFpEF / sum(isHFpEF);
    
    % Other group
    [uniqueValsOther, ~, idxOther] = unique(Stable{isOther, varName});
    countsOther = accumarray(idxOther, 1);
    percentagesOther = 100 * countsOther / sum(isOther);
    
    % Convert unique values to cell array of character vectors
    uniqueValsHFpEF = cellstr(string(uniqueValsHFpEF));
    uniqueValsOther = cellstr(string(uniqueValsOther));
    
    % Create a table for the HFpEF group
    hfpefTable = array2table([countsHFpEF'; percentagesHFpEF'], ...
                             'VariableNames', uniqueValsHFpEF', ...
                             'RowNames', {'Counts', 'Percentages'});
    
    % Create a table for the Other group
    otherTable = array2table([countsOther'; percentagesOther'], ...
                             'VariableNames', uniqueValsOther', ...
                             'RowNames', {'Counts', 'Percentages'});
    
    % Store the tables in the structure
    results.HFpEF.(varName) = hfpefTable;
    results.Other.(varName) = otherTable;
end