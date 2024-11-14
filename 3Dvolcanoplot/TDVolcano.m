%% Script Summary:
% This script generates a 3D Volcano plot.

% Created and designed by Feng Gu
% Last modified: 10/29/2024

clear
%% Load input features
PredictorsT = readtable("BaselineFeaturesNoPHI.xlsx",'VariableNamingRule','preserve');
RiskLevel = PredictorsT.Risk;

% Dummy variable encoding
PredictorsT.HFtype = strcmp(PredictorsT.HFtype, 'HFpEF');
PredictorsT.HFtype = double(PredictorsT.HFtype);  

PredictorsT.Race = strcmp(PredictorsT.Race, 'Caucasian');
PredictorsT.Race = double(PredictorsT.Race);  

PredictorsT.Smoking = strcmp(PredictorsT.Smoking, 'Current');
PredictorsT.Smoking = double(PredictorsT.Smoking);  

PredictorsT.Gender(PredictorsT.Gender == 2) = 0;  

PredictorsT.Alcohol = strcmp(PredictorsT.Alcohol, 'Yes');
PredictorsT.Alcohol = double(PredictorsT.Alcohol); 

PredictorsT.Drug = strcmp(PredictorsT.Drug, 'Yes');
PredictorsT.Drug = double(PredictorsT.Drug);  

columnsToConvert = PredictorsT.Properties.VariableNames;
for i = 4:length(columnsToConvert)
    colName = columnsToConvert{i};
    if iscell(PredictorsT.(colName))
    PredictorsT.(colName) = str2double(PredictorsT.(colName));
    end
end

% Verify features and remove any duplicates, all-zero, and all-NaN columns if present
[~,locs] = unique(round(table2array(mean(PredictorsT(:,4:end-10))),4));
locs = locs+3;
locs = sort(locs);
locs = [1;locs];
PredictorsT = PredictorsT(:,locs);
colsToRemove_NaN = varfun(@(x) all(isnan(x)), PredictorsT, 'OutputFormat', 'uniform');
colsToRemove_Zero = varfun(@(x) all(x == 0), PredictorsT, 'OutputFormat', 'uniform');
colsToRemove = colsToRemove_NaN | colsToRemove_Zero;
PredictorsT = PredictorsT(:, ~colsToRemove);
rng(0);

%% replace the out of bound values by 99% confidence lower bound and upper bound (just for parameters)
startCol = find(strcmp(PredictorsT.Properties.VariableNames, 'C_SA'));
endCol = find(strcmp(PredictorsT.Properties.VariableNames, 'R_a_c'));
data_normalized = zscore(PredictorsT{:,startCol:endCol});
z_threshold = 3;
outliers1 = data_normalized > z_threshold;
outliers2 = data_normalized < -z_threshold;
[~,PredictorsNumber] = size(data_normalized);
for i = 1:PredictorsNumber
    PredictorsT{:,i+startCol}(outliers1(:,i)) = max(PredictorsT{:,i+startCol}(~outliers1(:,i)));
    PredictorsT{:,i+startCol}(outliers2(:,i)) = min(PredictorsT{:,i+startCol}(~outliers2(:,i)));
end
binaryCols = varfun(@(x) all(ismember(x, [0 1])), PredictorsT, 'OutputFormat', 'uniform');
nonBinaryCols = ~binaryCols;
A_NonBinary = PredictorsT{:, nonBinaryCols};
ANorm_NonBinary = zscore(A_NonBinary(:,2:end), [], 'omitnan');
ANorm_NonBinary = [A_NonBinary(:,1) ANorm_NonBinary];
PredictorsT{:, nonBinaryCols} = ANorm_NonBinary;
ANorm_Optp = PredictorsT{:, 2:end};
PredictorName  = PredictorsT.Properties.VariableNames(:,2:end);
Lable = {'C_S_A','C_S_V','C_P_A','C_P_V','kact_L_V','kact_R_V','kpas_L_V','kpas_R_V','Vw_L_V','Vw_R_V','Vw_S_E_P','Amref_L_V','Amref_R_V','Amref_S_E_P',...
    'V4c','K_p_c','B_p_c','R_S_A','R_P_A','R_t_S_A','R_t_P_A','R_A_t','R_t_c','R_p_o','R_p_c','R_m_o','R_m_c','R_a_o','R_a_c'};
PredictorName(:,startCol-1:endCol-1) = Lable;
PredictorsT.Properties.VariableNames(:,2:end) = PredictorName;

%% Preparing data for 3D valcono plot
% Convert the ANorm_Optp array into a data table and name the columns with the PredictorName array
dataTbl = array2table(ANorm_Optp, 'VariableNames', PredictorName);
dataTbl.RiskLevel = categorical(RiskLevel, 'Ordinal', true);

% Perform ANOVA analysis for each predictor variable separately and collect the p-values
isContinuous = varfun(@(x) isnumeric(x) && length(unique(x)) > 2, dataTbl, 'OutputFormat', 'uniform');
numContinuous = sum(isContinuous)-1;
isBinary = varfun(@(x) isnumeric(x) && all(ismember(unique(x), [0, 1])), dataTbl, 'OutputFormat', 'uniform');
numBinary = sum(isBinary);
volcano3DData = array2table(zeros(width(dataTbl)-1, 7), 'VariableNames', {'P','RD12','P12','RD23','P23','RD31','P31'});
volcano3DData.Properties.RowNames = dataTbl.Properties.VariableNames(1:end-1);
numPredictors = width(dataTbl) - 1; % Exclude the RiskLevel variable

for i = 1:numPredictors
    % Select the i-th predictor variable
    predictor = dataTbl{:, i};

    % Find indices of rows that do not contain NaN values
    validRows = ~isnan(predictor);
    % Use index to exclude rows with NaN values
    validPredictor = predictor(validRows);
    validRiskLevel = dataTbl.RiskLevel(validRows);

    % Check if the predictor is a binary variable (0 or 1)
    if all(ismember(unique(validPredictor), [0, 1]))
        % Perform Chi-Square test for binary variables
        [tbl, ~, p] = crosstab(validPredictor, validRiskLevel);
        
        % Store the adjusted p-value in volcano3DData
        volcano3DData{dataTbl.Properties.VariableNames{i}, "P"} = min(p * (numContinuous+numBinary), 1);

        % Post-hoc pairwise comparisons
        gnames = categories(validRiskLevel);
        idx = [1 2; 2 3; 3 1];

        for k = 1:size(idx, 1)
            subTbl = tbl(:, idx(k, :));
            [~, ~, subP] = crosstab(validPredictor, validRiskLevel == gnames{idx(k, 1)} | validRiskLevel == gnames{idx(k, 2)});
            volcano3DData{dataTbl.Properties.VariableNames{i}, sprintf("RD%d%d", idx(k, 1), idx(k, 2))} = 5*(-subTbl(2, 1)/(subTbl(1, 1)+subTbl(2, 1)) + subTbl(2, 2)/(subTbl(1, 2)+subTbl(2, 2)));  % Difference in proportions
            volcano3DData{dataTbl.Properties.VariableNames{i}, sprintf("P%d%d", idx(k, 1), idx(k, 2))} = min(subP * nchoosek(length(gnames), 2) * numBinary, 1);
        
        end
    else
        % Perform ANOVA on data without NaN values for continuous variables
        [p, tbl, stats] = anova1(validPredictor, validRiskLevel, 'off');
        % Multiple comparisons
        [results, ~, ~, gnames] = multcompare(stats, "CType", "bonferroni", "Display", "off");
        volcano3DData{dataTbl.Properties.VariableNames{i}, "P"} = min([p * (numContinuous+numBinary), 1]);
        volcano3DData{dataTbl.Properties.VariableNames{i}, "RD12"} = -results(1, 4);
        volcano3DData{dataTbl.Properties.VariableNames{i}, "P12"} = min([results(1, 6) * numContinuous, 1]);
        volcano3DData{dataTbl.Properties.VariableNames{i}, "RD23"} = -results(3, 4);
        volcano3DData{dataTbl.Properties.VariableNames{i}, "P23"} = min([results(3, 6) * numContinuous, 1]);
        volcano3DData{dataTbl.Properties.VariableNames{i}, "RD31"} = results(2, 4);
        volcano3DData{dataTbl.Properties.VariableNames{i}, "P31"} = min([results(2, 6) * numContinuous, 1]);
    end  
end

%% Define property for figures

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

%% Set up for 3D volcano plot.
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
labels = {'PG2-PG1', 'PG3-PG2', 'PG1-PG3'}; % Labels corresponding to each angle
offset = 0.05 * maxRadius;  % Text offset from max radius, adjust as needed
% Draw lines parallel to the Z-axis and concentric circles
minZ = 0;      % Minimum cylinder height
maxZ = 60;     % Maximum cylinder height
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
radii = [1, 2, 3, maxRadius];  % Define radii for the concentric circles

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

zLabelPosition = maxZ - 25; 

text(xZAxis-0.5, yZAxis, zLabelPosition, '-Log_{10}P', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 20, ...
    'FontName', 'Times New Roman', ...
    'Rotation',90,...
    'Interpreter', 'tex');

set(gca, 'FontName', 'Times New Roman');

% Hide X and Y axis labels
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
% Alternatively, keep the Z-axis visible if needed
ax.ZAxis.Visible = 'off';

% Ensure x and y axis scales are equal and symmetric
xlim([-1.1*maxRadius, 1.1*maxRadius]);
ylim([-1.1*maxRadius, 1.1*maxRadius]);

% Set aspect ratio to 1:1:1.33 for x, y, and z axes
pbaspect([1, 1, 1.11]);
% Set the view angle for the plot
view([-10, 20]);
% Set the color of the axes to none (transparent background)
set(gca, 'Color', 'none');
dcm = datacursormode(gcf);  % Get the data cursor mode object for the current figure
set(dcm, 'UpdateFcn', {@myupdatefcn, volcano3DData,angles});  % Customize the data cursor
% Release the hold on the current plot
hold off;

%% Set up the video writer
% Set up the video writer
videoName = 'my3DAnimation.avi';
vidObj = VideoWriter(videoName);
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);  % Open video file

% Set the aspect ratio to 16:9 (e.g., 1280x720)

set(gcf, 'Position', [100, 100, 1280, 720]);  % Set figure size

% Rotate the figure around the Z-axis and capture frames
for t = linspace(0, 360, 361)
    view(t, 20);  % Update view angle
    drawnow;      % Refresh figure
    frame = getframe(gcf);  % Capture frame
    writeVideo(vidObj, frame);  % Write frame to video
end

% Finish and close the video file
close(vidObj);  % Close video file

%% Create figure for the base of the 3D volcano plot

figure(3);
clf;  % Clear current figure

% Create a polar axes for the data points and labels
polarAxes = polaraxes;
hold(polarAxes, 'on');
labelOffset = maxRadius * 1.05;  % Label offset from the center
sectorcolors = jet(length(sectorangles)-1); % Colors

% Draw labels and data points
for i = 1:length(angles)
    polarplot(polarAxes, [angles(i), angles(i)], [0, maxRadius], 'k-', 'LineWidth', 3);  % Draw black axis lines
    % Add text labels
    if i == 1
        text(angles(i), 1.2 * labelOffset, labels{i},...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle',...
            'FontSize', 28,... % Adjust font size
            'FontWeight', 'normal',...
            'Parent', polarAxes); %
    else
        text(angles(i), 1.38 * labelOffset, labels{i},...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle',...
            'FontSize', 28,... % Adjust font size
            'FontWeight', 'normal',...
            'Parent', polarAxes); %
    end
end

theta_final = zeros(height(volcano3DData), 1);
rho_final = zeros(height(volcano3DData), 1);

for i = 1:height(volcano3DData)
    RD_values = table2array(volcano3DData(i, {'RD12', 'RD23', 'RD31'}));

    [x_vals, y_vals] = pol2cart(angles, RD_values);
    x_final = sum(x_vals, 'omitnan');
    y_final = sum(y_vals, 'omitnan');
    [theta_final(i), rho_final(i)] = cart2pol(x_final, y_final);

    polarplot(polarAxes, theta_final(i), rho_final(i), 'o', 'MarkerSize', 28,...
        'MarkerFaceColor', color(i, :), 'MarkerEdgeColor', [0.85, 0.85, 0.85]);
end

polarAxes.ThetaTickLabel = {};  % Remove theta tick labels to declutter the plot
polarAxes.RTickLabel = {};  % Remove Radias tick labels to declutter the plot
polarAxes.FontSize = 28;
polarAxes.LineWidth = 2;
polarAxes.RLim = [0 1.2*maxRadius];  % Set radius limits for the polar axis
set(polarAxes, 'Color', 'none');  % Set the color of the axes to none (transparent)

% Create a Cartesian axes for the color ring
cartesianAxes = axes;
set(cartesianAxes, 'Position', get(polarAxes, 'Position'));  % Match the position with the polar axes
set(cartesianAxes, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');  % Hide axes
hold(cartesianAxes, 'on');

% Draw the colored ring with adjusted width
drawColoredRing(cartesianAxes, sectorangles, colororder, maxRadius, 0.5);  % Increase width to 0.2

% Link axes properties for synchronizing zoom/pan
linkprop([polarAxes, cartesianAxes], 'Position');
axis(cartesianAxes, 'equal');
cartesianAxes.XTick = [];  % Remove X ticks
cartesianAxes.YTick = [];  % Remove Y ticks
cartesianAxes.XColor = 'none';  % Hide X axis
cartesianAxes.YColor = 'none';  % Hide Y axis

dcm = datacursormode(gcf);  % Get the data cursor mode object for the current figure
set(dcm, 'UpdateFcn', {@myupdatefcn, volcano3DData, theta_final, rho_final});  % Customize data cursor

hold(polarAxes, 'off');
hold(cartesianAxes, 'off');

%% Define heatmap sorting priority based on 3D valcano plot group color

indices = (1:width(dataTbl)-1).';  % Create column vector of indices

% Create the 'Source' column
source = [
    repmat({'data'}, startCol-1, 1);
    repmat({'parameters'}, endCol - startCol, 1);
    repmat({'prediction'}, width(PredictorsT)-endCol, 1)
    ];
T = table(source, color, distance, indices, 'VariableNames', {'Source', 'Color', 'Distance', 'Index'});
T.Properties.RowNames = PredictorName;

% Define priorities for 'Color' and 'Source' sorting
colorPriority = [red; orange;[1,0,0]; yellow;deepPink;  [1,1,0]; lightGray;lightPurple; green;deepPurple ;lightblue; [0 0 1];blue];
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
        T.DistanceOrder(i) = -T.Distance(i); % Negate to sort in descending order
    elseif T.ColorOrder(i) == 7 % For lightGray
        T.DistanceOrder(i) = NaN; % Set NaN to exclude from sorting
    else % For green-cyan-blue
        T.DistanceOrder(i) = T.Distance(i); % Keep as is to sort in ascending order
    end
end

% Sort the table by 'Source', 'Color', and 'Distance'
[sortedT, sortOrder] = sortrows(T, {'SourceOrder', 'ColorOrder','DistanceOrder'});

% Retrieve the index post sorting
sequence = sortedT.Index; % Final sorted sequence of the original index

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
group3Avg = mean(group3Data(:, end-6:end), 2); % Row averages for the last 9 columns (group 3)

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

% Create a heatmap with empty labels for readability
% make sure the row height is the same for each subplot
% Identify columns with specific suffixes and other columns
D_Columns = find(endsWith(PredictorName, '_D'));
D_Columns = (1:D_Columns(end));
S_Columns = find(endsWith(PredictorName, '_S'));
totalColumns = size(Data, 2);
otherColumns = (D_Columns(end)+1:S_Columns(1)-1);

% Calculate the number of units for each group
numD = length(D_Columns);
numS = length(S_Columns);
numOther = length(otherColumns);
totalUnits = numD + numS + numOther;

heightPerUnit = 0.984 / totalUnits;
subplot('Position', [0.1, 9e-3 + heightPerUnit*numS+heightPerUnit*numOther, 0.8, heightPerUnit*numD]);
hm1 = heatmap(Data(:,D_Columns)','YDisplayLabels', repmat({''}, size(Data(:,D_Columns), 2), 1),...
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
subplot('Position', [0.1, 7e-3 + heightPerUnit*numS, 0.8, heightPerUnit*numOther]);
hm2 = heatmap(Data(:,otherColumns)','YDisplayLabels', repmat({''}, size(Data(:,otherColumns), 2), 1),...
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
subplot('Position', [0.1, 5e-3, 0.8, heightPerUnit*numS]);
hm3 = heatmap(Data(:,S_Columns)','YDisplayLabels', repmat({''}, size(Data(:,S_Columns), 2), 1),...
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

%% Making heatmap on selective predictors

% Remove NaNs from sortedT.DistanceOrder and create a new sequence
Newsequence = sequence(~isnan(sortedT.DistanceOrder));

% Extract data based on the processed sequence
SelectedData = ANorm_Optp(:, Newsequence); % for data

% Create and clear figure
figure(9); clf;
set(gcf, 'Position', [525, 50, 400, 1050]); % Adjust figure size as [left, bottom, width, height]

% Extract data for different risk groups
group1Data = SelectedData(RiskLevel == 1, :); % Data for risk group 1
group3Data = SelectedData(RiskLevel == 3, :); % Data for risk group 3

% Calculate row averages for risk groups 1 and 3
group1Avg = mean(group1Data(:, 1:20), 2); % Row averages for the first 8 columns (group 1)
group3Avg = mean(group3Data(:, end-16:end), 2); % Row averages for the last 9 columns (group 3)

% Sort rows within risk groups 1 and 3 based on their average values
[~, group1SortOrder] = sort(group1Avg, 'ascend'); % Order group 1 by ascending mean value
[~, group3SortOrder] = sort(group3Avg, 'descend'); % Order group 3 by descending mean value

% Reorder rows within each risk group
group1Data = group1Data(group1SortOrder, :);
group3Data = group3Data(group3SortOrder, :);

% Place reordered group data back into the main dataset
SelectedData(RiskLevel == 1, :) = group1Data;
SelectedData(RiskLevel == 3, :) = group3Data;

% Sort data by risk level in ascending order
[~, classifier] = sort(RiskLevel, 'ascend');
SelectedData = SelectedData(classifier, :);
% Identify columns with specific suffixes and other columns
D_Columns = find(endsWith(PredictorName(Newsequence), '_D'));
D_Columns = (1:D_Columns(end));
S_Columns = find(endsWith(PredictorName(Newsequence), '_S'));
totalColumns = size(SelectedData, 2);
otherColumns = (D_Columns(end)+1:S_Columns(1)-1);

% Calculate the number of units for each group
numD = length(D_Columns);
numS = length(S_Columns);
numOther = length(otherColumns);
totalUnits = numD + numS + numOther;

% Calculate height per unit for subplot positioning
heightPerUnit = 0.984 / totalUnits;

% Plot heatmaps for each group of columns
% Heatmap for _D columns
subplot('Position', [0.1, 9e-3 + heightPerUnit*numOther + heightPerUnit*numS, 0.64, heightPerUnit*numD]);
hm1 = heatmap(SelectedData(:, D_Columns)', ...
    'XDisplayLabels', repmat({''}, size(SelectedData, 1), 1), 'MissingDataColor', [0.8 0.8 0.8]);

% Interpolate colors for the colormap to have a smooth gradient
originalSize = size(redbluecmap, 1);
newColorMap = interp1(1:originalSize, redbluecmap, linspace(1, originalSize, 256), 'linear');
hm1.Colormap = newColorMap;
hm1.YDisplayLabels = PredictorName(Newsequence(D_Columns));
hm1.FontSize = 10;
innerPosition = get(hm1, 'InnerPosition');
innerPosition(1) = innerPosition(1) + 0.25; % shift
set(hm1, 'InnerPosition', innerPosition);
hm1.ColorLimits = [-3 3]; % Set the color axis limits
hm1.GridVisible = 'off'; % Disable the grid for a cleaner look
hm1.ColorbarVisible = 'off';

% Heatmap for other columns
subplot('Position', [0.1, 7e-3 + heightPerUnit*numS, 0.64, heightPerUnit*numOther]);
hm3 = heatmap(SelectedData(:, otherColumns)', ...
    'YDisplayLabels', repmat({''}, size(SelectedData(:, otherColumns), 2), 1), ...
    'XDisplayLabels', repmat({''}, size(SelectedData, 1), 1), 'MissingDataColor', [0.8 0.8 0.8]);

hm3.Colormap = newColorMap;
hm3.YDisplayLabels = PredictorName(Newsequence(otherColumns));
hm3.FontSize = 10;
innerPosition = get(hm3, 'InnerPosition');
innerPosition(1) = innerPosition(1) + 0.25; % shift
set(hm3, 'InnerPosition', innerPosition);
hm3.ColorLimits = [-3 3]; % Set the color axis limits
hm3.GridVisible = 'off'; % Disable the grid for a cleaner look
hm3.ColorbarVisible = 'off';

% Heatmap for _S columns
subplot('Position', [0.1, 5e-3, 0.64, heightPerUnit*numS]);
hm2 = heatmap(SelectedData(:, S_Columns)', ...
    'YDisplayLabels', repmat({''}, size(SelectedData(:, S_Columns), 2), 1), ...
    'XDisplayLabels', repmat({''}, size(SelectedData, 1), 1), 'MissingDataColor', [0.8 0.8 0.8]);

hm2.Colormap = newColorMap;
hm2.YDisplayLabels = PredictorName(Newsequence(S_Columns));
innerPosition = get(hm2, 'InnerPosition');
innerPosition(1) = innerPosition(1) + 0.25; % shift
set(hm2, 'InnerPosition', innerPosition);
hm2.FontSize = 10;
hm2.ColorLimits = [-3 3]; % Set the color axis limits
hm2.GridVisible = 'off'; % Disable the grid for a cleaner look
hm2.ColorbarVisible = 'off';
