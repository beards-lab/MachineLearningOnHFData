function [sortedVariables, sequence] = sortKeyDrivers(model12, model23, predictorNames, Threshold, meanRiskGroup1, meanRiskGroup3)
    % Collect coefficients and p-values for two models
    coeffs = [model12.Coefficients.Value(2:end), model23.Coefficients.Value(2:end)];
    pValues = [model12.Coefficients.pValue(2:end), model23.Coefficients.pValue(2:end)];
    
    % Calculate the average of the coefficients
    avgCoeffs = mean(coeffs, 2);
    
    % Define variable categories based on the sign of coefficients and p-values
    cat1Indices = all(coeffs < 0, 2) & any(pValues < Threshold, 2);
    cat3InitialIndices = (sign(coeffs(:, 1)) ~= sign(coeffs(:, 2))) | all(pValues > Threshold, 2);
    cat5Indices = all(coeffs > 0, 2) & any(pValues < Threshold, 2);

    % Calculate ratio for reclassification based on threshold
    ratio = meanRiskGroup1' ./ meanRiskGroup3';
    reclassCat1ToCat3 = (ratio <= 2) & (ratio > 0);
    reclassCat5ToCat3 =  ratio >= 0.5;
    
    % Combine indices for categories 1 and 5 except those that need reclassification
    cat1Indices = cat1Indices & ~reclassCat1ToCat3;
    cat5Indices = cat5Indices & ~reclassCat5ToCat3;
    % Take initial cat3 indices and add those that are reclassified from cat1 and cat5
    cat3Indices = cat3InitialIndices | reclassCat1ToCat3|reclassCat5ToCat3;

    % Assign categories to variables
    variableCategories = zeros(size(avgCoeffs));
    variableCategories(cat1Indices) = 1;
    variableCategories(cat3Indices) = 3;
    variableCategories(cat5Indices) = 5;

    % Create a table of the data
    tableData = table(predictorNames, avgCoeffs, meanRiskGroup1', meanRiskGroup3', variableCategories, ...
                      'VariableNames', {'Variable', 'AverageCoeff', 'MeanRiskGroup1', 'MeanRiskGroup3', 'Category'});
                  
    % Separate tableData into categories 1, 3, and 5
    cat1 = tableData(tableData.Category == 1, :);
    cat3 = tableData(tableData.Category == 3, :);
    cat5 = tableData(tableData.Category == 5, :);
    
    % Apply initial sort based on risk level for categories 1 and 5
    cat1 = sortrows(cat1, 'MeanRiskGroup1', 'descend');
    cat3 = sortrows(cat3, 'MeanRiskGroup1', 'descend');
    cat5 = sortrows(cat5, 'MeanRiskGroup3', 'descend');
    
    % Sort by coefficients (secondary sort) for categories 1 and 5
    cat1 = [cat1, table((1:height(cat1))', 'VariableNames', {'RiskLevelOrder'})];
    cat1 = sortrows(cat1, 'AverageCoeff', 'ascend');
    cat1 = [cat1, table((1:height(cat1))', 'VariableNames', {'CoeffOrder'})];
    
    cat5 = [cat5, table((1:height(cat5))', 'VariableNames', {'RiskLevelOrder'})];
    cat5 = sortrows(cat5, 'AverageCoeff', 'descend');
    cat5 = [cat5, table((1:height(cat5))', 'VariableNames', {'CoeffOrder'})];
    
    % Calculate the average of the ranking orders (for categories 1 and 5)
    cat1.AvgOrder = mean([cat1.RiskLevelOrder, cat1.CoeffOrder], 2);
    cat5.AvgOrder = mean([cat5.RiskLevelOrder, cat5.CoeffOrder], 2);

    % Final sort based on the average ranking orders (for categories 1 and 5)
    cat1 = sortrows(cat1, {'AvgOrder', 'MeanRiskGroup1'}, {'ascend', 'descend'});
    cat5 = sortrows(cat5, {'AvgOrder', 'MeanRiskGroup3'}, {'descend', 'descend'});

    % Category 3 remains the same (already sorted by MeanRiskGroup1)
    
    % Combine the sorted tables into one
    sortedVariables = [cat1(:,1:5); cat3; cat5(:,1:5)];

    % Find the matching sequence for original predictor names within the sortedVariables
    [~, sequence] = ismember(sortedVariables.Variable, predictorNames);
end