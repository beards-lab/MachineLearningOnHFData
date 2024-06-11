function plotModelCoefficients(model, predictorNames, threshold)
    % Extract model coefficients and their p-values
    coeffs = table2array(model.Coefficients(2:end, 'Value'));
    pValues = table2array(model.Coefficients(2:end, 'pValue'));
    coeffsCI = coefCI(model);
    coeffsCI = exp(coeffsCI(2:end,:));
    oddsRatios = exp(coeffs);

    % Combine variable names, Odds Ratios, p-values, and Confidence Intervals into a table
    resultsTable = table(predictorNames', oddsRatios, coeffsCI(:,1), coeffsCI(:,2), pValues, ...
        'VariableNames', {'Variable', 'OddsRatio', 'LowerCI', 'UpperCI', 'pValue'});

    % Sort table by Odds Ratio size
    resultsTable = sortrows(resultsTable, 'OddsRatio');

    % Create a forest plot
    hold on;
    for i = 1:height(resultsTable)
        if resultsTable.pValue(i) < threshold
            if resultsTable.OddsRatio(i) > 1
                scatterClr = 'r'; % Positive and significant coefficients in red
            else
                scatterClr = 'b'; % Negative and significant coefficients in blue
            end
        else
            scatterClr = 'k'; % Non-significant coefficients in black
        end
        scatter(resultsTable.OddsRatio(i), i, 'MarkerEdgeColor', scatterClr, 'MarkerFaceColor', scatterClr,'SizeData',36);
        line([resultsTable.LowerCI(i) resultsTable.UpperCI(i)], [i i], 'Color', scatterClr);
    end

    % Set y-axis and font size
    ax = gca; % Get the handle to the current axis
    set(ax, 'YTick', 1:height(resultsTable), 'YTickLabel', resultsTable.Variable, 'FontSize', 12,'FontWeight','bold', 'TickDir', 'out');  % Set y-axis tick labels
    ax.XAxis.FontSize = 14;
    % Calculate original x-axis range
    minCI = min(resultsTable.LowerCI);
    maxCI = max(resultsTable.UpperCI);

    % Round the range to the nearest power of 10
    minLimit = 10^(floor(log10(minCI))-1);
    maxLimit = 10^(ceil(log10(maxCI))+1);

    % Apply the range to xlim
    xlim([minLimit, maxLimit]);

    % Set x-axis to logarithmic scale
    set(gca, 'XScale', 'log');
    xline(1, 'LineStyle', '--'); % Add a reference line at OR=1
    hold off;
end