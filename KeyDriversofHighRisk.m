%% **********************************************************************************
%           Find key drivers of High Risk HF patients
% ***********************************************************************************
%
%   High Risk based on the heriarchial and kms clusters. Union of cluster 3
%
%   Script created on:  3 May  2024
%   Last modified on:   10 December 2019
%
%   Developed by        Feng Gu
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************

%close all
%clc
clear all
%% Load potential predictors values
PredictorsT = readtable("Predictors.xlsx",'VariableNamingRule','preserve');
% Predictslot = [(1:6) (11:21) (24:26) (35:38) (43:48) (50:56) (59:70)];
% Predictslot = [(1:6) (11:21) (24:26) (35:38) (43:48) (50:56) 70];
% Predictslot = [(1:38) (59:69)];
% PredictorsT = PredictorsT(:,Predictslot);
RiskLevel = PredictorsT{:,end};
% check out std and delete the repeat ones
[~,locs] = unique(round(table2array(mean(PredictorsT(:,3:end))),4));
locs = locs+2;
locs = sort(locs);
locs = [1;2;locs(1:end-1)];
PredictorsT = PredictorsT(:,locs);
rng(0);

% SIMULATING MISSING DATA
% Generate random creatinine levels based on a normal distribution with the
% mean set at 1.1 mg/dL and a standard deviation that represents a reasonable range
creatinineMean = 1.1;           % Mean creatinine level
creatinineStd = 0.2;            % Standard deviation for creatinine
creatinine = creatinineMean + creatinineStd * randn(346, 1); % Generate random values

% Ensure all generated creatinine values are positive by setting any negative
% values to a small positive number (e.g., 0.1)
creatinine(creatinine < 0) = 0.1;

% Add generated creatinine levels to the table of predictors
PredictorsT.Creatinine = creatinine;

% Simulate categorical and binary predictors using random integer generation
PredictorsT.("NYHA class") = randi([1, 4], 346, 1);  % NYHA class: values from 1 to 4
PredictorsT.("First diagnosis of heart failure in the past 18 months") = randi([0, 1], 346, 1);  % Binary: 0 or 1
PredictorsT.('Usage of beta-blocker') = randi([0, 1], 346, 1);  % Binary: 0 or 1
PredictorsT.('Usage of ACEI') = randi([0, 1], 346, 1);  % Binary: 0 or 1
% replace the out of bound values by 95% confidence lower bound and upper
% bound
data_normalized = zscore(PredictorsT{:,77:106});
z_threshold = 3;
outliers1 = data_normalized > z_threshold;
outliers2 = data_normalized < -z_threshold;
[~,PredictorsNumber] = size(data_normalized);
for i = 1:PredictorsNumber
    PredictorsT{:,i+76}(outliers1(:,i)) = max(PredictorsT{:,i+76}(~outliers1(:,i)));
    PredictorsT{:,i+76}(outliers2(:,i)) = min(PredictorsT{:,i+76}(~outliers2(:,i)));
end
A_Optp = PredictorsT{:,3:end};
ANorm_Optp = zscore(A_Optp,[],'omitnan');
PredictorName  = PredictorsT.Properties.VariableNames(:,3:end);
Lable = {'C_S_A','C_S_V','C_P_A','C_P_V','kact_L_V','kact_R_V','kpas_L_V','kpas_R_V','Vw_L_V','Vw_R_V','Vw_S_E_P','V_a_u','V_a_c',...
    'V4c','K_p_c','B_p_c','V_v_s','R_S_A','R_P_A','R_t_S_A','R_t_P_A','R_A_t','R_t_o','R_t_c','R_p_o','R_p_c','R_m_o','R_m_c','R_a_o','R_a_c'};
PredictorName(:,75:104) = Lable;
% Lable = {'C_S_A','C_S_V','C_P_A','C_P_V','G_S_A','G_P_A','G_t_S_A','G_t_P_A','kact_L_V','kact_R_V','kpas_L_V','kpas_R_V',...
%      'Vw_L_V','Vw_R_V','Vw_S_E_P','V_a_u','V_a_c','G_a','V4c','K_p_e_r_i','B_p_e_r_i','G_t_o','G_t_c','G_p_o','G_p_c','G_m_o','G_m_c','G_a_o','G_a_c','V_v_s'};
% PredictorName(:,1:30) = Lable;

%% Finding key drivers via logistic regression analysis for high-risk groups compared to low-risk group
%
% % Ordinal Logistic Regression
% % Convert ANorm_Optp into a data table dataTbl
% dataTbl = array2table(ANorm_Optp, 'VariableNames', PredictorName);
%
% % Add the categorical RiskLevel to the data table
% dataTbl.RiskLevel = categorical(RiskLevel, 'Ordinal', true);
% dataTbl = sortrows(dataTbl,"RiskLevel","descend");
% % Use the fitmnr function with the data table to fit an ordinal logistic regression model
% model = fitmnr(dataTbl, 'RiskLevel', 'ModelType', 'nominal');
% ci = coefCI(model);
% % ci = ci(3:end,:);
% ci = ci([(2:40) (42:end)],:);
% % coefftable = model.Coefficients(3:end,:);
% coefftable = model.Coefficients([(2:40) (42:end)],:);
% coefftable.OddsRatio = exp(coefftable.Value);
% coefftable.CI95Lower = exp(ci(:,1));
% coefftable.CI95Upper = exp(ci(:,2));
%
% coefftable.LogPValue = -log10(coefftable.pValue);
% coefftable.LogPValue(coefftable.Value<0) = -coefftable.LogPValue(coefftable.Value<0);
% coefftable = sortrows(coefftable, "OddsRatio", "ascend");
%
% % Convert coefficients and their standard errors to Odds Ratios (OR) and their 95% Confidence Intervals
% oddsRatios = coefftable.OddsRatio;
% CI95Lower = coefftable.CI95Lower;
% CI95Upper = coefftable.CI95Upper;
%
% % Create a forest plot
% figure;
% for i = 1:length(oddsRatios)
%     line([CI95Lower(i) CI95Upper(i)], [i i], 'Color', 'k'); % 95% CI line
%     hold on;
% end
% plot(oddsRatios, 1:length(oddsRatios), 'ro', 'MarkerFaceColor', 'r'); % OR dots
%
% % Set y-axis labels and title
% set(gca, 'YTick', 1:length(oddsRatios), 'YTickLabel', coefftable.Properties.RowNames, 'FontSize', 14);
% xlabel('Odds Ratio', 'FontSize', 14);
%
% % Set logarithmic scale
% set(gca, 'XScale', 'log');
% % Set the X-axis ticks and tick labels
% set(gca, 'XTick', [0.01 0.1 1 10 100]);
% set(gca, 'XTickLabel', {'0.01', '0.1', '1', '10', '100'});
%
% % Add a vertical line to represent OR=1
% line([1, 1], [0, length(oddsRatios) + 1], 'Color', 'b', 'LineStyle', '--');
% % ylim([-1 80]);
% % Reverse y-axis direction
% set(gca, 'YDir', 'reverse');
% hold off;

% % Separate analysis for low, middle, and high-risk groups
% % Perform logistic regression for RiskLevel 1 vs 2
% % Identify the indices of instances with RiskLevel 1 or 2
% index12 = RiskLevel == 1 | RiskLevel == 2;
% % Create a binary response variable where class 2 becomes 0, and class 1
% % becomes 1
% binaryResponse12 = double(RiskLevel(index12) == 1);
%
% % Logistic regression for 1 vs 2 risk levels comparison
% % Fit a multinomial logistic regression model using the relevant subset of data
% model12 = fitmnr(ANorm_Optp(index12,:), binaryResponse12,'Tolerance',1e-5);
%
% % Perform logistic regression for RiskLevel 1 vs 3
% % Identify the indices of instances with RiskLevel 1 or 3
% index13 = RiskLevel == 1 | RiskLevel == 3;
% % Create binary response variable where class 3 becomes 0, and class 1
% % becomes 1
% binaryResponse13 = double(RiskLevel(index13) == 1);
%
% % Logistic regression for 1 vs 3 risk levels comparison
% % Fit a multinomial logistic regression model using the relevant subset of data
% model13 = fitmnr(ANorm_Optp(index13,:), binaryResponse13,'Tolerance',1e-5);
%
% % Perform logistic regression for RiskLevel 2 vs 3
% % Identify the indices of instances with RiskLevel 2 or 3
% index23 = RiskLevel == 2 | RiskLevel == 3;
% % Create binary response variable where class 3 becomes 0, and class 2
% % becomes 1
% binaryResponse23 = double(RiskLevel(index23) == 2);
%
% % Logistic regression for 2 vs 3 risk levels comparison
% % Fit a multinomial logistic regression model using the relevant subset of data
% model23 = fitmnr(ANorm_Optp(index23,:), binaryResponse23,'Tolerance',1e-5);
%
%
% % Bonferroni correction
% Threshold = 0.05/39;
%
% % Visualize the coefficients of the models to identify key drivers
% % Each plotModelCoefficients call will create a forest plot for the respective logistic regression model
%
% % Visualize the coefficients for the 1 vs 2 model
% figure(1);clf;
% set(gcf, 'Position', [800, 50, 600, 1000]); % Adjust figure size as [left, bottom, width, height]
% plotModelCoefficients(model12, PredictorName,Threshold);
% % Visualize the coefficients for the 1 vs 3 model
% figure(2);clf;
% set(gcf, 'Position', [200, 50, 600, 1000]); % Adjust figure size as [left, bottom, width, height]
% plotModelCoefficients(model13, PredictorName,Threshold);
% % Visualize the coefficients for the 2 vs 3 model
% figure(3);clf;
% set(gcf, 'Position', [200, 50, 600, 1000]); % Adjust figure size as [left, bottom, width, height]
% plotModelCoefficients(model23, PredictorName,Threshold);
%
%
% % select key drivers
% % Call a custom function and sort the variables according to rules
% [sortedVariables,sequence] = sortKeyDrivers(model12, model13,PredictorName', Threshold,meanRiskGroup1,meanRiskGroup3);
%
% % Display the sorted variable names
% disp(sortedVariables);
%% Preparing data for 3D valcono and ternary plot
% Convert the ANorm_Optp array into a data table and name the columns with the PredictorName array
dataTbl = array2table(ANorm_Optp, 'VariableNames', PredictorName);
dataTbl.RiskLevel = categorical(RiskLevel, 'Ordinal', true);

% Retrieve the names of all variables in the table
variableNames = dataTbl.Properties.VariableNames;

% Iterate over each column to adjust negative values if present
for i = 1:length(variableNames)
    % Extract the data for the current column
    columnData = dataTbl.(variableNames{i});

    % Check if the column data is numeric
    if isnumeric(columnData)
        % Find the minimum value in the column
        minVal = min(columnData);

        % If the minimum value is less than 0, offset it to make the minimum 0
        if minVal < 0
            % Update the data in the table by adding the absolute value of the minimum
            dataTbl.(variableNames{i}) = columnData - minVal;
        end
    end
end

% Prepare an empty table for ternary plot data with three variables
ternaryplotData = array2table(zeros(149, 3), 'VariableNames', {'G1','G2','G3'});
% Assign row names based on the variable names of 'dataTbl', excluding the last variable
ternaryplotData.Properties.RowNames = dataTbl.Properties.VariableNames(1:end-1);

% Iterate over each variable (each column excluding the last one)
for i = 1:width(dataTbl)-1
    % Calculate the mean value for each group within the current variable
    Feature = dataTbl{:,i};
    meanG1 = mean(Feature(RiskLevel == 1),'omitnan');
    meanG2 = mean(Feature(RiskLevel == 2),'omitnan');
    meanG3 = mean(Feature(RiskLevel == 3),'omitnan');

    % Check for NaN and replace with the mean of the other two groups
    if isnan(meanG1)
        meanG1 = (meanG2 + meanG3) / 2;
    end

    if isnan(meanG2)
        meanG2 = (meanG1 + meanG3) / 2;
    end

    if isnan(meanG3)
        meanG3 = (meanG1 + meanG2) / 2;
    end

    % Calculate the sum of the means
    sumMeans = meanG1 + meanG2 + meanG3;

    % Calculate and assign the ratios for ternary plot data
    ternaryplotData.G1(i) = meanG1 / sumMeans;
    ternaryplotData.G2(i) = meanG2 / sumMeans;
    ternaryplotData.G3(i) = meanG3 / sumMeans;
end

% Perform ANOVA analysis for each predictor variable separately and collect the p-values
numPredictors = width(dataTbl) - 1; % Exclude the RiskLevel variable count
volcano3DData = array2table(zeros(149, 7), 'VariableNames', {'P','RD12','P12','RD23','P23','RD31','P31'});
volcano3DData.Properties.RowNames = dataTbl.Properties.VariableNames(1:end-1);

for i = 1:numPredictors
    % Select the i-th predictor variable
    predictor = dataTbl{:, i};

    % Find indices of rows that do not contain NaN values
    validRows = ~isnan(predictor);
    % Use index to exclude rows with NaN values
    validPredictor = predictor(validRows);
    validRiskLevel = dataTbl.RiskLevel(validRows);

    % Perform ANOVA on data without NaN values
    [p, tbl, stats] = anova1(validPredictor, validRiskLevel, 'off');
    % Multiple comparisons
    [results, ~, ~, gnames] = multcompare(stats, "CType", "bonferroni", "Display", "off");
    volcano3DData{dataTbl.Properties.VariableNames(i),"P"} = min([p*numPredictors, 1]);
    volcano3DData{dataTbl.Properties.VariableNames(i),"RD12"} = -results(1,4);
    volcano3DData{dataTbl.Properties.VariableNames(i),"P12"} = min([results(1,6)*numPredictors, 1]);
    volcano3DData{dataTbl.Properties.VariableNames(i),"RD23"} = -results(3,4);
    volcano3DData{dataTbl.Properties.VariableNames(i),"P23"} = min([results(3,6)*numPredictors, 1]);
    volcano3DData{dataTbl.Properties.VariableNames(i),"RD31"} = results(2,4);
    volcano3DData{dataTbl.Properties.VariableNames(i),"P31"} = min([results(2,6)*numPredictors, 1]);
end
%% define property for figures
maxRD = max([max(volcano3DData.RD12), max(volcano3DData.RD23), max(volcano3DData.RD31)]);
maxRadius = 1.8 * maxRD;  % This is used to determine the length of the polar axis.
angles = [pi/2, 2*pi/3 + pi/2, 4*pi/3 + pi/2]; % These are the angles for polar coordinates.
% Define colors for plotting
lightGray = [0.7 0.7 0.7];
blue = [0,0.45,0.74];
red = [0.84,0.08,0.18];
green = [0.4660 0.6740 0.1880];
lightblue = [0 0.75 0.75];
orange = [0.8500 0.3250 0.0980];
yellow = [0.93,0.69,0.13];
lightPurple = [0.7, 0.2, 0.7];
deepPink = [1.0, 0.08, 0.58];
deepPurple = [0.48, 0.19, 0.8];
colororder = [red;orange;yellow;[1 1 0]; green; lightblue;blue;[0 0 1];deepPurple;lightPurple;deepPink;[1 0 0];lightGray];
sectorangles = (0:pi/6:2*pi)+pi/12;
%% Visualization of Polar Coordinate System Classification Regions
% Create a polar scatter plot
figure(1); hold on;
clf;  % Clear current figure
polaraxes;  % Create a polar axes
hold on;  % Retain plots so that new plots added to the axes do not delete existing plots
pax = gca;  % Get the current polar axes
% Add labels at specified angles

for i = 1:length(sectorangles)-1
    theta = [sectorangles(i) sectorangles(i+1)];
    raddi1 = [1 maxRadius];
    raddi2 = [0 1];
    polarregion(theta,raddi2,FaceColor=lightGray,FaceAlpha=0.75)
    if i >= 6
        polarregion(theta,raddi1,FaceColor=colororder(12-i+6,:),FaceAlpha=0.75)
    else
        polarregion(theta,raddi1,FaceColor=colororder(1-i+5,:),FaceAlpha=0.75)
    end
end


% Remove angle labels
pax.ThetaTickLabel = {};  % Remove theta tick labels to declutter the plot
pax.RTickLabel = {};  % Remove Radias tick labels to declutter the plot

pax.FontSize = 28;
pax.LineWidth = 2;
% Set the polar axis limit
pax.RLim = [0 maxRadius];  % Set radius limits for the polar axis
set(pax, 'Color', 'none');  % Set the color of the axes to none (transparent)

% Release the hold on the current plot
hold off;
% Set up for 3D volcano plot.
figure(2); clf; % Clear the current figure.
hold on; % Retain plots so that new plots added to the figure do not delete existing plots.

% Prepare the z-coordinates.
alpha = 0.05; % Significance level for determining the z reference (threshold).
z_refer = -log10(alpha); % Compute the reference z-coordinate for the given alpha significance level.
color = zeros(height(volcano3DData), 3); % Initialize a matrix to hold RGB color values for each data point.
distance = zeros(height(volcano3DData), 1);% used for sorting
% Loop through each row, calculating polar coordinates and corresponding z-coordinates for each data point
for i = 1:height(volcano3DData)
    % Extract the RD and p-values for the current data point
    RD_values = table2array(volcano3DData(i, {'RD12','RD23','RD31'}));
    p_values = table2array(volcano3DData(i, {'P','P12','P23','P31'}));

    % Convert polar to Cartesian coordinates
    [x_vals, y_vals] = pol2cart(angles, RD_values);
    z_vals = -log10(p_values);
    % Sum the x and y vectors to get the final point
    x_final = sum(x_vals, 'omitnan');
    y_final = sum(y_vals, 'omitnan');
    z_coord = z_vals(1);
    % Convert Cartesian coordinates back to polar to determine angle and radius
    [theta, r] = cart2pol(x_final, y_final);
    theta = mod(theta, 2*pi);
    [~, ~, binIndex] = histcounts(theta, sectorangles);
    if binIndex >= 6
        color(i,:)=colororder(12-binIndex+6,:);
    else
        color(i,:)=colororder(1-binIndex+5,:);
    end
    % Determine which p-value to use based on the angle and set color and z_coord accordingly

    if z_coord <= z_refer || r <= 1
        color(i,:) = lightGray;
    end

    % Plot the 3D scatter point
    scatter3(x_final, y_final, z_coord, 380, color(i,:), 'filled','MarkerEdgeColor',[0.75,0.75,0.75],'AlphaData',0.5);
    distance(i)  = sqrt(x_final^2+y_final^2+z_coord^2);% used for sorting
end
% Draw polar coordinate axes through the origin
labels = {'G2-G1', 'G3-G2', 'G1-G3'}; % Labels corresponding to each angle
offset = 0.05 * maxRadius;  % Text offset from max radius, adjust as needed
% Draw lines parallel to the Z-axis and concentric circles
minZ = 0;      % Minimum cylinder height
maxZ = 30;     % Maximum cylinder height
for i = 1:length(angles)
    % Choose a long axis length, here assumed to be the max RD value from the data
    [x_end, y_end] = pol2cart(angles(i), maxRadius);
    line([0, x_end], [0, y_end], [0, 0], 'Color', 'k', 'LineWidth', 2.5);
    % Calculate position along the extended line outside the circular plot
    % this is used for ouput figure
    if i == 1
        xPos = (maxRadius + i*offset) * cos(angles(i))+1e-1;
        yPos = (maxRadius + i*offset) * sin(angles(i));
    elseif i ==2
        xPos = (maxRadius + 6*offset) * cos(angles(i));
        yPos = (maxRadius + 6*offset) * sin(angles(i));
    else
        xPos = (maxRadius + 5*offset) * cos(angles(i));
        yPos = (maxRadius + 5*offset) * sin(angles(i));
    end
    % This is used for movie
    xPos = (maxRadius + 4*offset) * cos(angles(i));
    yPos = (maxRadius + 4*offset) * sin(angles(i));
    % Add text annotations at the specified position

    text(xPos, yPos, minZ, labels{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 20, ...
        'FontWeight', 'normal', ...
        'Color', 'k');
end

thetaLines = 0 + pi/2 : pi/3 : 2*pi + pi/2;  % Angle intervals for vertical lines

% Draw vertical lines parallel to the Z-axis
for theta = thetaLines
    xLine = maxRadius * cos(theta);
    yLine = maxRadius * sin(theta);
    line([xLine, xLine], [yLine, yLine], [minZ, maxZ], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
end

% Draw concentric circles at the bottom and top
thetaCircle = linspace(0, 2*pi, 1000);  % For a smoother circle
radii = [1, 1.5, 2, maxRadius];  % Define radii for the concentric circles

for r = radii
    xCircle = r * cos(thetaCircle);
    yCircle = r * sin(thetaCircle);

    % Bottom circle
    plot3(xCircle, yCircle, repmat(minZ, size(xCircle)), 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5);
end

% Top circle
plot3(xCircle, yCircle, repmat(maxZ, size(xCircle)), 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5);

% Radius labels are to be placed at an angle of thetaLabel
thetaLabel = 0; % Assuming the desire to place radius labels at pi/4 radians

% Draw radius labels
for r = radii(1:end-1)
    xLabel = (r + 0.05) * cos(thetaLabel); % Position between the edge and diameter of the circle
    yLabel = (r + 0.05) * sin(thetaLabel); % Slightly outside on the circle to place the label
    text(xLabel, yLabel, minZ, num2str(r), ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 20, ...
        'Color', 'k');
end

% Custom Z-axis properties
thetaZAxis = pi;  % Position of Z-axis in polar coordinates
% X and Y coordinates of the intersection at the bottom of the Z-axis
xZAxis = maxRadius * cos(thetaZAxis);
yZAxis = maxRadius * sin(thetaZAxis);

% Draw custom Z-axis
line([xZAxis, xZAxis], [yZAxis, yZAxis], [minZ, maxZ], 'Color', 'k', 'LineWidth', 2.5);

% Simulating Z-axis tick marks
zTicks = minZ : 10 : maxZ; % For example, a tick mark every 10 units
lenTick = 0.05 * maxRadius;  % Length of tick marks, in the outward direction
% Add tick marks
for zTick = zTicks
    % Calculate the endpoint of the tick mark, along the polar axis direction
    xTick = xZAxis + lenTick * cos(thetaZAxis);
    yTick = yZAxis + lenTick * sin(thetaZAxis);
    % Draw small horizontal lines for each tick mark
    plot3([xZAxis, xTick], [yZAxis, yTick], [zTick, zTick], 'k-', 'LineWidth', 1.5);

    % Add numeric labels at the end of each tick mark
    text(xTick, yTick, zTick, num2str(zTick), ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', 'FontSize', 20);
end

% Hide X and Y axis labels
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
% Alternatively, keep the Z-axis visible if needed
ax.ZAxis.Visible = 'off';

% Ensure x and y axis scales are equal and symmetric
xlim([-1.1*maxRadius, 1.1*maxRadius]);
ylim([-1.1*maxRadius, 1.1*maxRadius]);

% Set aspect ratio to 1:1:1 for x, y, and z axes
pbaspect([1, 1, 1]);
% Set the view angle for the plot
view([5, 20]);
% Set the color of the axes to none (transparent background)
set(gca, 'Color', 'none');
% Release the hold on the current plot
hold off;
% Set up the video writer
videoName = 'my3DAnimation.avi';
vidObj = VideoWriter(videoName);
vidObj.Quality = 100;
vidObj.FrameRate = 24;
open(vidObj);  % Open video file

% Rotate the figure around the Z-axis and capture frames
for t = linspace(0, 360, 361)
    view(t, 20);  % Update view angle
    drawnow;      % Refresh figure
    frame = getframe(gcf);  % Capture frame
    writeVideo(vidObj, frame);  % Write frame to video
end

% Finish and close the video file
close(vidObj);  % Close video file


% Draw polar coordinate system
% Create a polar scatter plot
figure(3);
clf;  % Clear current figure
polaraxes;  % Create a polar axes
hold on;  % Retain plots so that new plots added to the axes do not delete existing plots
pax = gca;  % Get the current polar axes
labelOffset = maxRadius * 1.05;  % Label offset from the center
% Add labels at specified angles
for i = 1:length(angles)
    polarplot([angles(i) angles(i)], [0 maxRadius], 'k-', 'LineWidth', 3);  % Draw black axis lines
    % Add text labels
    if i == 1
        text(angles(i), labelOffset, labels{i},...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle',...
            'FontSize', 28,... % Adjust font size
            'FontWeight', 'normal',...
            'Parent',pax); %
    else
        text(angles(i), 1.08*labelOffset, labels{i},...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle',...
            'FontSize', 28,... % Adjust font size
            'FontWeight', 'normal',...
            'Parent',pax); %
    end
end

% Calculate the position of each point and plot
theta_final = zeros(height(volcano3DData), 1);
rho_final = zeros(height(volcano3DData), 1);

for i = 1:height(volcano3DData)
    RD_values = table2array(volcano3DData(i, {'RD12', 'RD23', 'RD31'}));

    [x_vals, y_vals] = pol2cart(angles, RD_values);
    x_final = sum(x_vals, 'omitnan');
    y_final = sum(y_vals, 'omitnan');
    [theta_final(i), rho_final(i)] = cart2pol(x_final, y_final);

    polarplot(theta_final(i), rho_final(i), 'o', 'MarkerSize', 28,...
        'MarkerFaceColor', color(i, :), 'MarkerEdgeColor', [0.85 0.85 0.85]);
end

% Remove angle labels
pax.ThetaTickLabel = {};  % Remove theta tick labels to declutter the plot
pax.RTickLabel = {};  % Remove Radias tick labels to declutter the plot
pax.FontSize = 28;
pax.LineWidth = 2;
% Set the polar axis limit
pax.RLim = [0 maxRadius];  % Set radius limits for the polar axis
set(pax, 'Color', 'none');  % Set the color of the axes to none (transparent)

% Enable interactive data cursor mode
dcm = datacursormode(gcf);  % Get the data cursor mode object for the current figure
set(dcm, 'UpdateFcn', {@myupdatefcn, volcano3DData, theta_final, rho_final});  % Customize the data cursor

% Release the hold on the current plot
hold off;

%% Visualization of Ternary Phase Diagram
% Draw the ternary phase diagram and the translated angle bisectors, ensuring that the lines do not exceed the boundaries of the triangle
% Define the coordinates of the three vertices (equilateral triangle)
vertex1 = [0, 0]; % Coordinates of the first point (x1, y1)
vertex2 = [1, 0]; % Coordinates of the second point (x2, y2)
sideLength = norm(vertex2 - vertex1); % Length of the triangle side
H = (sqrt(3) / 2) * sideLength; % Height of the equilateral triangle

vertex3 = [(vertex1(1) + vertex2(1)) / 2, H]; % Calculate the coordinates of the third point (x3, y3)
figure(4); clf;
hold on;

% Draw the three sides of the triangle
plot([vertex1(1), vertex2(1)], [vertex1(2), vertex2(2)], 'k-');
plot([vertex1(1), vertex3(1)], [vertex1(2), vertex3(2)], 'k-');
plot([vertex2(1), vertex3(1)], [vertex2(2), vertex3(2)], 'k-');

% Set the interval for the scale
interval = 0.1; % Define the interval for the scale, here it is 0.1

% Calculate the coordinates of the scale points and draw the scale lines and scale marks
for frac = interval:interval:1-interval
    % Calculate the scale points from vertex1 to vertex3
    x13 = vertex1(1) + frac * (vertex3(1) - vertex1(1));
    y13 = vertex1(2) + frac * (vertex3(2) - vertex1(2));

    % Calculate the scale points from vertex2 to vertex3
    x23 = vertex2(1) + frac * (vertex3(1) - vertex2(1));
    y23 = vertex2(2) + frac * (vertex3(2) - vertex2(2));

    x12 = vertex1(1) + frac;
    x21 = vertex1(1)+ 1 - frac;
    % Draw the scale lines parallel to the base
    plot([x13, x23], [y13, y23], 'k:');
    plot([x13, x12], [y13, 0], 'k:');
    plot([x23, x21], [y23, 0], 'k:');
end

% Calculate the midpoints of the three sides
midpoint12 = (vertex1 + vertex2) / 2;
midpoint23 = (vertex2 + vertex3) / 2;
midpoint31 = (vertex3 + vertex1) / 2;

% Calculate the boundary lines of the internal area
bound1 = polyfit([vertex1(1), vertex2(1)], [vertex1(2), vertex2(2)], 1);
bound2 = polyfit([vertex2(1), vertex3(1)], [vertex2(2), vertex3(2)], 1);
bound3 = polyfit([vertex3(1), vertex1(1)], [vertex3(2), vertex1(2)], 1);

% Calculate the direction vectors of the angle bisectors
dir1 = midpoint23 - vertex1;
dir2 = midpoint31 - vertex2;
dir3 = midpoint12 - vertex3;

% Calculate the perpendicular unit vectors in the direction of the angle bisectors
perpDir1 = [-dir1(2), dir1(1)] / norm(dir1);
perpDir2 = [-dir2(2), dir2(1)] / norm(dir2);
perpDir3 = [-dir3(2), dir3(1)] / norm(dir3);

% Translation distance
shiftDist = 0.05;
% Define the three boundary lines of the triangle
boundaries = [vertex1; vertex2; vertex3; vertex1]; % Line from vertex1 to vertex2

% Draw the translated angle bisectors, ensuring the lines are within the triangle
intersection1 = draw_clipped_line(midpoint23, vertex1, shiftDist, perpDir1, boundaries);
intersection2 = draw_clipped_line(midpoint23, vertex1, -shiftDist, perpDir1, boundaries);
intersection3 = draw_clipped_line(midpoint31, vertex2, shiftDist, perpDir2, boundaries);
intersection4 = draw_clipped_line(midpoint31, vertex2, -shiftDist, perpDir2, boundaries);
intersection5 = draw_clipped_line(midpoint12, vertex3, shiftDist, perpDir3, boundaries);
intersection6 = draw_clipped_line(midpoint12, vertex3, -shiftDist, perpDir3, boundaries);

% Collect the coordinates of all endpoints of the line segments
lines = {...
    [intersection1(1,:); intersection1(2,:)], ...
    [intersection2(1,:); intersection2(2,:)], ...
    [intersection3(1,:); intersection3(2,:)], ...
    [intersection4(1,:); intersection4(2,:)], ...
    [intersection5(1,:); intersection5(2,:)], ...
    [intersection6(1,:); intersection6(2,:)] ...
    };

% Initialize the set of intersection points
line_intersections = [];

% Compare each pair of line segments and find the intersection points
for i = 1:length(lines)
    for j = i+1:length(lines)
        % Get the i-th line segment
        seg1 = lines{i};
        % Get the j-th line segment
        seg2 = lines{j};

        % Calculate the intersection points of the two line segments using the modified lineSegmentIntersect function
        [xi, yi] = lineSegmentIntersect(seg1, seg2);

        % If there are intersection points, add them to the set
        if ~isempty(xi)
            line_intersections = [line_intersections; xi, yi];
        end
    end
end
sorted_intersections = sortrows(line_intersections, [1, 2]);
% Assume x and y are the coordinates of the vertices of your polygon
x1 = [vertex1(1) intersection1(1,1) sorted_intersections(1,1) sorted_intersections(3,1) sorted_intersections(2,1) intersection2(1,1)];
y1 = [vertex1(2) intersection1(1,2) sorted_intersections(1,2) sorted_intersections(3,2) sorted_intersections(2,2) intersection2(1,2)];
x2 = [intersection6(1,1)  sorted_intersections(2,1) intersection2(1,1)];
y2 = [intersection6(1,2) sorted_intersections(2,2) intersection2(1,2)];
x3 = [intersection6(1,1)  sorted_intersections(2,1) sorted_intersections(7,1) sorted_intersections(8,1) intersection5(1,1)];
y3 = [intersection6(1,2) sorted_intersections(2,2) sorted_intersections(7,2) sorted_intersections(8,2) intersection5(1,2)];
x4 = [sorted_intersections(8,1) intersection5(1,1) intersection3(2,1)];
y4 = [sorted_intersections(8,2) intersection5(1,2) intersection3(2,2)];
x5 = [vertex2(1) intersection4(2,1) sorted_intersections(12,1) sorted_intersections(9,1) sorted_intersections(8,1) intersection3(2,1)];
y5 = [vertex2(2) intersection4(2,2) sorted_intersections(12,2) sorted_intersections(9,2) sorted_intersections(8,2) intersection3(2,2)];
x6 = [intersection4(2,1)  sorted_intersections(12,1) intersection2(2,1)];
y6 = [intersection4(2,2) sorted_intersections(12,2) intersection2(2,2)];
x7 = [intersection2(2,1)  sorted_intersections(12,1) sorted_intersections(10,1) sorted_intersections(11,1) intersection1(2,1)];
y7 = [intersection2(2,2) sorted_intersections(12,2) sorted_intersections(10,2) sorted_intersections(11,2) intersection1(2,2)];
x8 = [sorted_intersections(11,1) intersection1(2,1) intersection5(2,1)];
y8 = [sorted_intersections(11,2) intersection1(2,2) intersection5(2,2)];
x9 = [vertex3(1) intersection5(2,1) sorted_intersections(11,1) sorted_intersections(6,1) sorted_intersections(4,1) intersection6(2,1)];
y9 = [vertex3(2) intersection5(2,2) sorted_intersections(11,2) sorted_intersections(6,2) sorted_intersections(4,2) intersection6(2,2)];
x10 = [intersection6(2,1)  sorted_intersections(4,1) intersection4(1,1)];
y10 = [intersection6(2,2) sorted_intersections(4,2) intersection4(1,2)];
x11 = [intersection4(1,1)  sorted_intersections(4,1) sorted_intersections(5,1) sorted_intersections(1,1) intersection3(1,1)];
y11 = [intersection4(1,2) sorted_intersections(4,2) sorted_intersections(5,2) sorted_intersections(1,2) intersection3(1,2)];
x12 = [sorted_intersections(1,1) intersection3(1,1) intersection1(1,1)];
y12 = [sorted_intersections(1,2) intersection3(1,2) intersection1(1,2)];
x = [sorted_intersections(1,1) sorted_intersections(3,1) sorted_intersections(2,1) sorted_intersections(7,1) sorted_intersections(8,1) sorted_intersections(9,1) sorted_intersections(12,1) sorted_intersections(10,1) sorted_intersections(11,1) sorted_intersections(6,1) sorted_intersections(4,1) sorted_intersections(5,1)];
y = [sorted_intersections(1,2) sorted_intersections(3,2) sorted_intersections(2,2) sorted_intersections(7,2) sorted_intersections(8,2) sorted_intersections(9,2) sorted_intersections(12,2) sorted_intersections(10,2) sorted_intersections(11,2) sorted_intersections(6,2) sorted_intersections(4,2) sorted_intersections(5,2)];

% Fill the triangles with colors
fill(x,y,lightGray,'FaceAlpha',0.75)
fill(x1, y1, 'b','FaceAlpha',0.75);
fill(x2, y2, blue,'FaceAlpha',0.75);
fill(x3, y3, lightblue,'FaceAlpha',0.75);
fill(x4, y4, green,'FaceAlpha',0.75);
fill(x5, y5, 'y','FaceAlpha',0.75);
fill(x6, y6, yellow,'FaceAlpha',0.75);
fill(x7, y7, orange,'FaceAlpha',0.75);
fill(x8, y8, red,'FaceAlpha',0.75);
fill(x9, y9, 'r','FaceAlpha',0.75);
fill(x10, y10, deepPink,'FaceAlpha',0.75);
fill(x11, y11, lightPurple,'FaceAlpha',0.75);
fill(x12, y12, deepPurple,'FaceAlpha',0.75);

axis equal;
axis off;
xlim([min([vertex1(1), vertex2(1), vertex3(1)])-0.1, max([vertex1(1), vertex2(1), vertex3(1)])+0.1]);
ylim([min([vertex1(2), vertex2(2), vertex3(2)])-0.1, max([vertex1(2), vertex2(2), vertex3(2)])+0.1]);

hold off;



%% Extract three ratios from ternary phase diagram
colors = zeros(size(ternaryplotData, 1), 3);

% Iterate through data points
for i = 1:size(ternaryplotData, 1)
    G1 = ternaryplotData{i, 1};
    G2 = ternaryplotData{i, 2};
    G3 = ternaryplotData{i, 3};

    % Calculate differences
    diff12 = abs(G1 - G2);
    diff23 = abs(G2 - G3);
    diff13 = abs(G1 - G3);

    % Determine color
    if (diff12 < 0.1 && diff23 < 0.1) || (diff12 < 0.1 && diff13 < 0.1) || (diff23 < 0.1 && diff13 < 0.1)
        % Rule 1: At least two group differences less than 0.1
        colors(i, :) = lightGray; % Gray
    elseif diff12 < 0.1
        % Rule 2: Difference between G1 and G2 less than 0.1
        if G3 > max(G1, G2)
            colors(i, :) = [1, 0, 0]; % Red
        elseif G3 < min(G1, G2)
            colors(i, :) = lightblue; % Light blue
        end
    elseif diff23 < 0.1
        % Rule 3: Difference between G2 and G3 less than 0.1
        if G1 > max(G2, G3)
            colors(i, :) = [0, 0, 1]; % Blue
        elseif G1 < min(G2, G3)
            colors(i, :) = orange; % Orange
        end
    elseif diff13 < 0.1
        % Rule 4: Difference between G1 and G3 less than 0.1
        if G2 > max(G1, G3)
            colors(i, :) = [1, 1, 0]; % Yellow
        elseif G2 < min(G1, G3)
            colors(i, :) = lightPurple; % Light purple
        end
    else
        % Rule 5: Other cases
        if G3 > G2 && G2 > G1
            colors(i, :) = red; % Deep red
        elseif G3 > G1 && G1 > G2
            colors(i, :) = deepPink; % Deep pink
        elseif G1 > G3 && G3 > G2
            colors(i, :) = deepPurple; % Purple
        elseif G1 > G2 && G2 > G3
            colors(i, :) = blue; % Blue
        elseif G2 > G1 && G1 > G3
            colors(i, :) = green; % Green
        elseif G2 > G3 && G3 > G1
            colors(i, :) = yellow; % Deep yellow
        end
    end
end
figure(5); clf;
colors(volcano3DData.P > 0.05, 1) = 0.7;
colors(volcano3DData.P > 0.05, 2) = 0.7;
colors(volcano3DData.P > 0.05, 3) = 0.7;
ternaryplotData.c = colors;
ternaryplotData.P = volcano3DData.P;

% Initialize the index array of same length as the number of rows in colors
colorIndex = zeros(size(colors, 1), 1);

% Find the index for each color in colors within the colororder array
for i = 1:size(colors, 1)
    for j = 1:size(colororder, 1)
        if isequal(colors(i, :), colororder(j, :))
            colorIndex(i) = j;
            break;
        end
    end
end

% Sort the colorIndex while keeping track of the original indices
[~, sortedIndices] = sort(colorIndex);

% Create a sorted version of ternaryplotData
sortedTernaryplotData = ternaryplotData(sortedIndices, :);

% Create a new empty table to save sorted results
sortedDataGroups = cell(size(colororder, 1), 1);

for i = 1:size(colororder, 1)
    % Find all rows of the current color group
    colorGroupIndices = ismember(sortedTernaryplotData.c, colororder(i,:), 'rows');
    colorGroupData = sortedTernaryplotData(colorGroupIndices, :);

    % Sort by group and sorting rules
    if any(i == [1 2 11 12 13]) % First 2 and last 3 groups
        colorGroupData = sortrows(colorGroupData, 'G3', 'descend');
    elseif any(i == [3 4 5]) % Groups 3 to 5
        colorGroupData = sortrows(colorGroupData, 'G2', 'descend');
    elseif any(i == [6 7 8 9]) % Groups 6 to 9
        colorGroupData = sortrows(colorGroupData, 'G1', 'descend');
    elseif i == 10 % Tenth group
        colorGroupData = sortrows(colorGroupData, 'G2', 'ascend');
    end

    % Save the sorted data segment
    sortedDataGroups{i} = colorGroupData;
end

% Combine the sorted data segments into the final sorted table
finalSortedTernaryplotData = vertcat(sortedDataGroups{:});
addpath("github_repo\");

% Loop through each data point for plotting
% for i = height(sortedTernaryplotData):-1:1
%     ternplot(sortedTernaryplotData.G1(i), sortedTernaryplotData.G2(i), sortedTernaryplotData.G3(i),-log10(ternaryplotData.P(i)), 'o', 'Color', sortedTernaryplotData.c(i, :), 'MarkerFaceColor', sortedTernaryplotData.c(i, :),'MarkerSize',16);
%     hold on; % Keep the image to plot the next data point
% end

% Predefine the minimum and maximum values for P
p_min = 0; % When P is 1
p_max = -log10(0.05) * 2; % When P is 1e-30

% Create a colormap from deep blue to yellow
c_map = colormap("default");

% Calculate color indices corresponding to -log10(P)
p_indices = round( ( -log10(finalSortedTernaryplotData.P) - p_min ) / (p_max - p_min) * (size(c_map, 1) - 1) ) + 1;

% Ensure indices are within the correct range
p_indices(p_indices < 1) = 1;
p_indices(p_indices > size(c_map, 1)) = size(c_map, 1);

% Loop to plot points
for i = height(sortedTernaryplotData):-1:1
    P_value = ternaryplotData.P(i);
    P_index = p_indices(i);
    P_color = c_map(P_index, :); % Get corresponding color from colormap

    ternplot(sortedTernaryplotData.G1(i), ...
        sortedTernaryplotData.G2(i), ...
        sortedTernaryplotData.G3(i), ...
        'o', ...
        'Color', P_color, ...
        'MarkerFaceColor', P_color, ...
        'MarkerSize', 16);
    hold on;
end

cbr = colorbar;
% Set the range of colorbar (should match the range of data mapped to colormap)
clim([p_min p_max]); % If the range of P values is from 1e-0 to 1e-30

% Set the position of the colorbar tick mark at 0.5, the specific value needs to be located according to the setting of caxis
cbr.Ticks = [-log10(0.05)]; 

% Remove tick labels from the colorbar, but keep the tick mark corresponding to 0.5
cbr.TickLabels{cbr.Ticks == -log10(0.05)} = 'P = 0.05';
cbr.TickLabels(cbr.Ticks ~= -log10(0.05)) = {''};

% If only the tick mark is to be displayed without labels
set(cbr, 'TickLabels', []); 

% Get the current Axes object
ax = gca;

% Remove ticks
set(ax, 'XTick', [], 'YTick', [], 'ZTick', []);
% ternlabel('G1 ratio', 'G2 ratio', 'G3 ratio');


%%
% Define priority for sorting by color
indices = (1:149).';  % Create column vector of indices

% Create the 'Source' column
source = [
    repmat({'data'}, 74, 1);
    repmat({'parameters'}, 30, 1);
    repmat({'prediction'}, 45, 1)
    ];
T = table(source, color, distance, indices, 'VariableNames', {'Source', 'Color', 'Distance', 'Index'});
T.Properties.RowNames = PredictorName;

% Define priorities for 'Color' and 'Source' sorting
colorPriority = [red; [1,0,0];orange; deepPink; yellow; [1,1,0]; lightGray;lightPurple; green;deepPurple ;lightblue; [0 0 1];blue];
sourcePriority = {'data', 'parameters', 'prediction'};

% Create sort order columns
[~, sourceOrder] = ismember(T.Source, sourcePriority);

% Convert RGB array to color name strings
strColors = cellfun(@(c) ...
    strcat(num2str(c(1)), '_', num2str(c(2)), '_', num2str(c(3))), ...
    num2cell(T.Color, 2), 'UniformOutput', false);
strColorPriority = cellfun(@(c) ...
    strcat(num2str(c(1)), '_', num2str(c(2)), '_', num2str(c(3))), ...
    num2cell(colorPriority, 2), 'UniformOutput', false);

% Convert color name strings to categorical variables
catColors = categorical(strColors);
catColorPriority = categorical(strColorPriority);

% Get row index of 'Color' within 'colorPriority'
[~, colorOrder] = ismember(catColors, catColorPriority);
T.SourceOrder = sourceOrder;
T.ColorOrder = colorOrder;

% Create a sorting order for 'Distance'
% Rule: red-orange-yellow (descending), lightGray (unaltered), green-cyan-blue (ascending)
for i = 1:height(T)
    if T.ColorOrder(i) <= 6 % For red-orange-yellow
        T.DistanceOrder(i) = T.Distance(i); % Negate to sort in descending order
    elseif T.ColorOrder(i) == 7 % For lightGray
        T.DistanceOrder(i) = NaN; % Set NaN to exclude from sorting
    else % For green-cyan-blue
        T.DistanceOrder(i) = T.Distance(i); % Keep as is to sort in ascending order
    end
end

% Sort the table by 'Source', 'Color', and 'Distance'
[sortedT, sortOrder] = sortrows(T, {'SourceOrder', 'DistanceOrder'});

% Retrieve the index post sorting
sequence = sortedT.Index; % Final sorted sequence of the original index
%%
% % Initialize forest plot data structure
% forestData = {};
% forestIdx = 1;
%
% % Correction threshold for multiple comparisons
% alpha_corrected = 0.05 / numPredictors;
% pValues = zeros(1, numPredictors);
% for i = 1:numPredictors
%     % Select the i-th predictor variable
%     predictor = dataTbl{:, i};
%
%     % Find row indices without NaN values
%     validRows = ~isnan(predictor);
%     % Exclude rows with NaN to obtain valid data
%     validPredictor = predictor(validRows);
%     validRiskLevel = dataTbl.RiskLevel(validRows);
%
%     % Perform ANOVA using data without NaN values
%     [p, tbl, stats] = anova1(validPredictor, validRiskLevel, 'off');
%     pValues(i) = p;
%     % Multiple comparisons
%     [results, ~, ~, gnames] = multcompare(stats, "CType", "bonferroni", "Display", "off");
%
%     % Extract relevant data and store in the forest plot data structure
%     for j = 1:size(results, 1)
%         comparisonName = sprintf('%s%d%d', PredictorName{i}, results(j,1), results(j,2));
%         estimate = results(j,4);
%         lowerCI = results(j,3);
%         upperCI = results(j,5);
%
%         % Set color based on significance
%         if results(j,6) < alpha_corrected
%             color = estimate < 0 ? 'r' : 'b';
%         else
%             color = [0.5, 0.5, 0.5]; % Grey for non-significant results
%         end
%
%         % Store in forestData
%         forestData{forestIdx, 1} = comparisonName; % Name
%         forestData{forestIdx, 2} = estimate;       % Estimate
%         forestData{forestIdx, 3} = lowerCI;        % Lower confidence interval
%         forestData{forestIdx, 4} = upperCI;        % Upper confidence interval
%         forestData{forestIdx, 5} = color;          % Color
%
%         forestIdx = forestIdx + 1;
%     end
% end
%
% % Sort forestData by the estimate values
% forestData = sortrows(forestData, 2);
%
% % Plot the forest plot
% figure(1); clf;
% set(gcf, 'Position', [0, 100, 500, 1000]); % Adjust figure size as necessary
% hold on;
%
% % Define the axes range
% yLimits = [-0.5, size(forestData, 1) + 0.5];
%
% % Draw a reference line at x = 0 (vertical line)
% line([0 0], yLimits, 'Color', 'k', 'LineStyle', '--');
%
% % Loop to draw horizontal confidence interval lines
% for i = 1:size(forestData, 1)
%     % Horizontal confidence interval line
%     line([forestData{i, 3}, forestData{i, 4}], [i, i], 'Color', forestData{i, 5}, 'LineWidth', 2);
%     % Draw a point at the midpoint of the confidence interval line (vertically)
%     plot(forestData{i, 2}, i, 'o', 'MarkerFaceColor', forestData{i, 5}, 'MarkerEdgeColor', 'none');
% end
%
% % Remove axis text labels while keeping tick marks
% set(gca, 'XTickLabel', {}, 'YTick', 1:size(forestData, 1), 'YTickLabel', {}, 'TickLength', [0.01, 0.035]);
%
% % Set the Y-axis tick marks to the outside
% ax = gca; % Get current axes
% ax.YAxis.TickDirection = 'out';
%
% % Set figure axis ranges
% xlim([-1.75, 1.75]);
% ylim(yLimits);
%
% % Release the plot hold
% hold off;
%
% % Initialize arrays for grouping
% groupNum = zeros(length(PredictorName), 1);
% sortIdx = zeros(length(PredictorName), 1);
% for i = 1:length(PredictorName)
%     % Find rows in forestData that match the current PredictorName
%     matches = strncmp(PredictorName{i}, forestData(:, 1), length(PredictorName{i}));
%     matchRows = find(matches);
%
%     % Retrieve the color data for matching rows
%     matchColors = forestData(matches, 5);
%
%     % Count the number of blues and reds
%     numBlues = sum(cellfun(@(c) isequal(c, 'b'), matchColors));
%     numReds = sum(cellfun(@(c) isequal(c, 'r'), matchColors));
%
%     % Grouping criteria
%     if numBlues == 3
%         groupNum(i) = 7; % Group 1: 3 blues
%     elseif numBlues == 2
%         groupNum(i) = 6; % Group 2: 2 blues
%     elseif numBlues == 0 && numReds == 0
%         groupNum(i) = 4; % Group 4: no blues, no reds
%     elseif numReds == 1
%         groupNum(i) = 3; % Group 5: 1 red
%     elseif numReds == 2
%         groupNum(i) = 2; % Group 6: 2 reds
%     elseif numReds == 3
%         groupNum(i) = 1; % Group 7: 3 reds
%     else
%         groupNum(i) = 5; % Group 3: 1 blue
%     end
%
%     % Set sorting index according to the rule above
%     % Now we need to consider both effect direction and group number
%     if mean([forestData{matchRows, 2}]) < 0
%         sortIdx(i) = min(matchRows); % For negative effects, take the smallest row
%     else
%         sortIdx(i) = max(matchRows); % For positive effects, take the largest row
%     end
% end
%
% % Combine group numbers and row numbers for sorting
% % Sort first by group number then within the same group by sortIdx
% [~, sequence] = sortrows([groupNum, sortIdx]);
%
% % Reorder PredictorName according to the sorted sequence
% sortedPredictorName = PredictorName(sequence);

%% Making heatmap on all predictors
% Extract data based on the processed sequence
Data = ANorm_Optp(:,sequence);% for data
% Data = ANorm_Optp(:,sequence(75:104));% for parameters
% Data = ANorm_Optp(:,sequence(105:end));% for predictions
% Create and clear figure
figure(8); clf;
set(gcf, 'Position', [500, 50, 300, 1000]); % Adjust figure size as [left, bottom, width, height]

% Extract data for the different risk groups
group1Data = Data(RiskLevel == 1, :);  % Data for risk group 1
group3Data = Data(RiskLevel == 3, :);  % Data for risk group 3

% Calculate the row average values for risk groups 1 and 3
group1Avg = mean(group1Data(:, 1:20), 2);    % Row averages for the first 8 columns (group 1)
group3Avg = mean(group3Data(:, end-16:end), 2); % Row averages for the last 9 columns (group 3)

% Sort the rows within risk groups 1 and 3 based on their average values
[~, group1SortOrder] = sort(group1Avg, 'ascend');  % Order group 1 by ascending mean value
[~, group3SortOrder] = sort(group3Avg, 'descend'); % Order group 3 by descending mean value

% Reorder the rows within each risk group
group1Data = group1Data(group1SortOrder, :);
group3Data = group3Data(group3SortOrder, :);

% Place reordered group data back into the main dataset
Data(RiskLevel == 1, :) = group1Data;
Data(RiskLevel == 3, :) = group3Data;

% Sort data by risk level in descending order
[~,classifier] = sort(RiskLevel,'ascend');
Data = Data(classifier,:);
% Data = flipud(Data);
% Data = fliplr(Data);
% Create a heatmap with empty labels for readability
% make sure the row height is the same for each subplot
totalUnits = 74 + 30 + 45;
heightPerUnit = 0.984 / totalUnits;
subplot('Position', [0.1, 9e-3 + heightPerUnit*45+heightPerUnit*30, 0.8, heightPerUnit*74]);
hm1 = heatmap(Data(:,1:74)','YDisplayLabels', repmat({''}, size(Data(:,1:74), 2), 1),...
    'XDisplayLabels', repmat({''}, size(Data, 1), 1),'MissingDataColor', [0.8 0.8 0.8]);

% Interpolate colors for the colormap to have smooth gradient
originalSize = size(redbluecmap, 1);
newColorMap = interp1(1:originalSize, redbluecmap, linspace(1, originalSize, 256), 'linear');
hm1.Colormap = newColorMap;
% boldLabels = cellfun(@(x) ['\bf{' x '}'], PredictorName(sequence), 'UniformOutput', false);
% hm.YDisplayLabels = boldLabels; % Apply predictor names for y-axis labels
hm1.FontSize = 18;
hm1.ColorLimits = [-3 3];  % Set the color axis limits
hm1.GridVisible = 'off'; % Disable the grid for a cleaner look
hm1.ColorbarVisible = 'off';
subplot('Position', [0.1, 7e-3 + heightPerUnit*45, 0.8, heightPerUnit*30]);
hm2 = heatmap(Data(:,75:104)','YDisplayLabels', repmat({''}, size(Data(:,75:104), 2), 1),...
    'XDisplayLabels', repmat({''}, size(Data, 1), 1),'MissingDataColor', [0.8 0.8 0.8]);

% Interpolate colors for the colormap to have smooth gradient
originalSize = size(redbluecmap, 1);
newColorMap = interp1(1:originalSize, redbluecmap, linspace(1, originalSize, 256), 'linear');
hm2.Colormap = newColorMap;
% boldLabels = cellfun(@(x) ['\bf{' x '}'], PredictorName(sequence), 'UniformOutput', false);
% hm.YDisplayLabels = boldLabels; % Apply predictor names for y-axis labels
hm2.FontSize = 18;
hm2.ColorLimits = [-3 3];  % Set the color axis limits
hm2.GridVisible = 'off'; % Disable the grid for a cleaner look
hm2.ColorbarVisible = 'off';
subplot('Position', [0.1, 5e-3, 0.8, heightPerUnit*45]);
hm3 = heatmap(Data(:,105:end)','YDisplayLabels', repmat({''}, size(Data(:,105:end), 2), 1),...
    'XDisplayLabels', repmat({''}, size(Data, 1), 1),'MissingDataColor', [0.8 0.8 0.8]);

% Interpolate colors for the colormap to have smooth gradient
originalSize = size(redbluecmap, 1);
newColorMap = interp1(1:originalSize, redbluecmap, linspace(1, originalSize, 256), 'linear');
hm3.Colormap = newColorMap;
% boldLabels = cellfun(@(x) ['\bf{' x '}'], PredictorName(sequence), 'UniformOutput', false);
% hm.YDisplayLabels = boldLabels; % Apply predictor names for y-axis labels
hm3.FontSize = 18;
hm3.ColorLimits = [-3 3];  % Set the color axis limits
hm3.GridVisible = 'off'; % Disable the grid for a cleaner look
hm3.ColorbarVisible = 'off';
%% Making heatmap on selective  predictors
Newsequence = sequence(~isnan(sortedT.DistanceOrder));
% Extract data based on the processed sequence
SelectedData = ANorm_Optp(:,Newsequence);% for data
% Create and clear figure
figure(9); clf;
set(gcf, 'Position', [525, 50, 300, 1000]); % Adjust figure size as [left, bottom, width, height]

% Extract data for the different risk groups
group1Data = SelectedData(RiskLevel == 1, :);  % Data for risk group 1
group3Data = SelectedData(RiskLevel == 3, :);  % Data for risk group 3

% Calculate the row average values for risk groups 1 and 3
group1Avg = mean(group1Data(:, 1:20), 2);    % Row averages for the first 8 columns (group 1)
group3Avg = mean(group3Data(:, end-16:end), 2); % Row averages for the last 9 columns (group 3)

% Sort the rows within risk groups 1 and 3 based on their average values
[~, group1SortOrder] = sort(group1Avg, 'ascend');  % Order group 1 by ascending mean value
[~, group3SortOrder] = sort(group3Avg, 'descend'); % Order group 3 by descending mean value

% Reorder the rows within each risk group
group1Data = group1Data(group1SortOrder, :);
group3Data = group3Data(group3SortOrder, :);

% Place reordered group data back into the main dataset
SelectedData(RiskLevel == 1, :) = group1Data;
SelectedData(RiskLevel == 3, :) = group3Data;

% Sort data by risk level in descending order
[~,classifier] = sort(RiskLevel,'ascend');
SelectedData = SelectedData(classifier,:);
% make sure every subplot row height is the same
totalUnits = 21 + 14 + 28;
heightPerUnit = 0.984 / totalUnits;
subplot('Position', [0.1, 9e-3 + heightPerUnit*28+heightPerUnit*14, 0.64, heightPerUnit*21]);
hm1 = heatmap(SelectedData(:,1:21)',...
    'XDisplayLabels', repmat({''}, size(SelectedData, 1), 1),'MissingDataColor', [0.8 0.8 0.8]);

% Interpolate colors for the colormap to have smooth gradient
originalSize = size(redbluecmap, 1);
newColorMap = interp1(1:originalSize, redbluecmap, linspace(1, originalSize, 256), 'linear');
hm1.Colormap = newColorMap;
% boldLabels1 = cellfun(@(x) ['\bf{' x '}'], PredictorName(Newsequence(1:21)), 'UniformOutput', false);
% hm1.YDisplayLabels = boldLabels1; % Apply predictor names for y-axis labels
hm1.YDisplayLabels = PredictorName(Newsequence(1:21));
hm1.FontSize = 10;
innerPosition = get(hm1, 'InnerPosition');
innerPosition(1) = innerPosition(1) + 0.25; % shift
set(hm1, 'InnerPosition', innerPosition);
hm1.ColorLimits = [-3 3];  % Set the color axis limits
hm1.GridVisible = 'off'; % Disable the grid for a cleaner look
hm1.ColorbarVisible = 'off';
subplot('Position', [0.1, 7e-3 + heightPerUnit*28, 0.64, heightPerUnit*14]);
hm2 = heatmap(SelectedData(:,22:35)','YDisplayLabels', repmat({''}, size(SelectedData(:,22:35), 2), 1),...
    'XDisplayLabels', repmat({''}, size(SelectedData, 1), 1),'MissingDataColor', [0.8 0.8 0.8]);

% Interpolate colors for the colormap to have smooth gradient
originalSize = size(redbluecmap, 1);
newColorMap = interp1(1:originalSize, redbluecmap, linspace(1, originalSize, 256), 'linear');
hm2.Colormap = newColorMap;
% boldLabels2 = cellfun(@(x) ['\bf{' x '}'], PredictorName(Newsequence(22:35)), 'UniformOutput', false);
% hm2.YDisplayLabels = boldLabels2; % Apply predictor names for y-axis labels
hm2.YDisplayLabels =  PredictorName(Newsequence(22:35));
innerPosition = get(hm2, 'InnerPosition');
innerPosition(1) = innerPosition(1) + 0.25; % shift
set(hm2, 'InnerPosition', innerPosition);
hm2.FontSize = 10;
hm2.ColorLimits = [-3 3];  % Set the color axis limits
hm2.GridVisible = 'off'; % Disable the grid for a cleaner look
hm2.ColorbarVisible = 'off';
subplot('Position', [0.1, 5e-3, 0.64, heightPerUnit*28]);
hm3 = heatmap(SelectedData(:,36:end)','YDisplayLabels', repmat({''}, size(SelectedData(:,36:end), 2), 1),...
    'XDisplayLabels', repmat({''}, size(SelectedData, 1), 1),'MissingDataColor', [0.8 0.8 0.8]);

% Interpolate colors for the colormap to have smooth gradient
originalSize = size(redbluecmap, 1);
newColorMap = interp1(1:originalSize, redbluecmap, linspace(1, originalSize, 256), 'linear');
hm3.Colormap = newColorMap;
% boldLabels3 = cellfun(@(x) ['\bf{' x '}'], PredictorName(Newsequence(36:end)), 'UniformOutput', false);
% hm3.YDisplayLabels = boldLabels3; % Apply predictor names for y-axis labels
hm3.YDisplayLabels = PredictorName(Newsequence(36:end));
hm3.FontSize = 10;
innerPosition = get(hm3, 'InnerPosition');
innerPosition(1) = innerPosition(1) + 0.25; % shift
set(hm3, 'InnerPosition', innerPosition);
hm3.ColorLimits = [-3 3];  % Set the color axis limits
hm3.GridVisible = 'off'; % Disable the grid for a cleaner look
hm3.ColorbarVisible = 'off';
% 这个函数输入线段起始点、结束点、偏移距离和三角形的三条边界
function intersections = draw_clipped_line(startPoint, endPoint, shiftDist, perpDir, boundaries)
% 计算平移的起始和终止点
shiftedStart = startPoint + perpDir * shiftDist;
shiftedEnd = endPoint + perpDir * shiftDist;
direction = shiftedEnd-shiftedStart;
entension_factor = 2;
shiftedStart = shiftedStart-entension_factor*direction;
shiftedEnd  = shiftedEnd + entension_factor*direction;
% 计算交点
intersections = find_intersection(shiftedStart, shiftedEnd, boundaries);

% 用计算得到的交点绘制线段
% 此处逻辑确保我们有至少两个交点
if size(intersections, 1) == 2
    plot(intersections(:, 1), intersections(:, 2), 'k-');
end
end

function intersections = find_intersection(shiftedStart, shiftedEnd, boundaries)
% 初始化交点数组
intersections = [];

% 使用 polyxpoly 查找交点
[xi, yi] = polyxpoly([shiftedStart(1), shiftedEnd(1)], [shiftedStart(2), shiftedEnd(2)], ...
    boundaries(:,1), boundaries(:,2));
% 收集所有交点
intersections = [intersections; xi(:), yi(:)];

% 如果有超过两个交点，说明线段经过了三角形的一个顶点
% 这里我们仅保留前两个点作为线段的起止点
if size(intersections, 1) > 2
    intersections = intersections(1:2, :);
end

% 按 x 坐标排序交点（如果是垂直线，则按 y 坐标排序)
if shiftedStart(1) == shiftedEnd(1)
    [~, order] = sort(intersections(:,2));  % 垂直线段，按 y 排序
else
    [~, order] = sort(intersections(:,1));  % 非垂直线段，按 x 排序
end
intersections = intersections(order, :);
end
function [x, y] = lineSegmentIntersect(seg1, seg2)
% 将线段表示为参数方程 p = p0 + t*(p1 - p0)
p0 = seg1(1,:);
p1 = seg1(2,:);
p2 = seg2(1,:);
p3 = seg2(2,:);

s10_x = p1(1) - p0(1);  % 线段1的x方向增量
s10_y = p1(2) - p0(2);  % 线段1的y方向增量
s32_x = p3(1) - p2(1);  % 线段2的x方向增量
s32_y = p3(2) - p2(2);  % 线段2的y方向增量

denom = s10_x * s32_y - s32_x * s10_y;
if denom == 0
    % denom为0时线段平行或共线，不相交
    x = [];
    y = [];
    return;
end

denomPositive = denom > 0;
s02_x = p0(1) - p2(1);
s02_y = p0(2) - p2(2);
s_numer = s10_x * s02_y - s10_y * s02_x;
if (s_numer < 0) == denomPositive
    % 不相交
    x = [];
    y = [];
    return;
end

t_numer = s32_x * s02_y - s32_y * s02_x;
if (t_numer < 0) == denomPositive
    % 不相交
    x = [];
    y = [];
    return;
end

if ((s_numer > denom) == denomPositive) || ((t_numer > denom) == denomPositive)
    % 不相交
    x = [];
    y = [];
    return;
end

% 相交
t = t_numer / denom;
x = p0(1) + (t * s10_x);
y = p0(2) + (t * s10_y);
end