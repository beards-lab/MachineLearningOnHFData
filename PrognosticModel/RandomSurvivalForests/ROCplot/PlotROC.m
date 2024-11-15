clear;
clc;

% List of file names
file_names = {'UmichDM_roc_results.xlsx', 'UmichD_roc_results.xlsx', 'UmichM_roc_results.xlsx', 'UmichMAGGIC_roc_results.xlsx', 'UmichMAGGIC_M_roc_results.xlsx'};
model_names = {'DM', 'D', 'M', 'MAGGIC', 'MAGGIC_with_M'};

% Initialize a structure to store all data
roc_data_struct = struct();

% Read data from each file and sheet
for file_idx = 1:length(file_names)
    file_name = file_names{file_idx};
    model_name = model_names{file_idx};

    % Get sheet names
    [~, sheet_names] = xlsfinfo(file_name);

    for sheet_idx = 1:length(sheet_names)
        sheet_name = sheet_names{sheet_idx};
        % Extract the prefix before the last underscore
        last_underscore_idx = find(sheet_name == '_', 1, 'last');
        sheet_name_prefix = sheet_name(1:last_underscore_idx-1);

        % Match the time point in the sheet name using a regex
        tokens = regexp(sheet_name, sprintf('%s_(\\d+\\.?\\d*)', sheet_name_prefix), 'tokens');

        time_point = str2double(tokens{1}{1});

        % Read data from the current sheet
        data = readtable(file_name, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');

        % Ensure field names are valid
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_'); % Replace dots

        % Save data to the structure
        roc_data_struct.(field_name).FPR = data.FPR;
        roc_data_struct.(field_name).TPR = data.TPR;
        roc_data_struct.(field_name).Model = model_name;
        roc_data_struct.(field_name).TimePoint = time_point;
    end
end

% Define unique time points and model names
timePoints = [0.5, 1];
models = model_names;

% Number of repetitions for each model at each time point
num_repeats = 25;

% Calculate the number of original data points
num_points = length(roc_data_struct.(['Model_', models{1}, '_Time_0_5']).FPR) / num_repeats;

% Define the number of new data points for interpolation
new_points = 350;

% Generate a new FPR vector for interpolation
new_FPR = linspace(0, 1, new_points);

% Initialize a structure for storing processed data
data_struct = struct();

% Process and interpolate data for each time point and model
for t = 1:length(timePoints)
    time_point = timePoints(t);

    for m = 1:length(models)
        model_name = models{m};

        % Generate a valid field name
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_'); % Replace dots

        % Retrieve FPR and TPR data
        currentModelDataFPR = roc_data_struct.(field_name).FPR;
        currentModelDataTPR = roc_data_struct.(field_name).TPR;

        % Handle special case for MAGGIC model (no repetitions)
        if strcmp(model_name, 'MAGGIC')
            num_points = length(currentModelDataFPR); % Use all data points
            num_repeats_corrected = 1; % No repetitions
        else
            num_repeats_corrected = num_repeats; % Use specified repetitions
            num_points = length(currentModelDataFPR) / num_repeats_corrected;
            currentModelDataFPR = reshape(currentModelDataFPR, num_points, num_repeats_corrected);
            currentModelDataTPR = reshape(currentModelDataTPR, num_points, num_repeats_corrected);
        end

        % Initialize a matrix to store interpolated TPR values
        new_TPR_matrix = zeros(new_points, num_repeats_corrected);

        % Perform interpolation for each repetition
        for i = 1:num_repeats_corrected
            if num_repeats_corrected == 1
                FPR_col = currentModelDataFPR;
                TPR_col = currentModelDataTPR;
            else
                FPR_col = currentModelDataFPR(:, i);
                TPR_col = currentModelDataTPR(:, i);
            end

            % Remove duplicate FPR values
            [FPR_unique, unique_idx] = unique(FPR_col);
            TPR_unique = TPR_col(unique_idx);

            % Ensure the range [0, 1] is included
            if FPR_unique(1) ~= 0
                FPR_unique = [0; FPR_unique];
                TPR_unique = [0; TPR_unique];
            end
            if FPR_unique(end) ~= 1
                FPR_unique = [FPR_unique; 1];
                TPR_unique = [TPR_unique; 1];
            end

            % Interpolate TPR values
            new_TPR_matrix(:, i) = interp1(FPR_unique, TPR_unique, new_FPR, 'linear');
            new_TPR_matrix(1, i) = 0; % Ensure the first value is 0
        end

        % Create a matrix for the new FPR vector
        new_FPR_matrix = repmat(new_FPR', 1, num_repeats_corrected);

        % Save interpolated data to the structure
        data_struct.(field_name).FPR = new_FPR_matrix;
        data_struct.(field_name).TPR = new_TPR_matrix;
        data_struct.(field_name).Model = model_name;
        data_struct.(field_name).TimePoint = time_point;
    end
end



%% Visualize data for each time point
colors = [
    119, 172, 48;
    237, 176, 33;
    160, 0, 0;
    0, 0, 0;
    0, 155, 189;
] / 255;

for t = 1:length(timePoints)
    time_point = timePoints(t);

    % Create a new figure window
    figure(t); clf;
    hold on;

    % Iterate through all models
    for m = length(models)-1:-1:1
        model_name = models{m};

        % Ensure the field name is valid
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_');

        % Retrieve FPR and TPR matrices
        FPR_matrix = data_struct.(field_name).FPR;
        TPR_matrix = data_struct.(field_name).TPR;

        % Calculate mean and standard error
        if strcmp(model_name, 'MAGGIC')
            TPR_mean = TPR_matrix;
            TPR_se = zeros(length(TPR_mean), 1); % No repetitions, standard error is zero
        else
            TPR_mean = mean(TPR_matrix, 2);
            TPR_se = std(TPR_matrix, 0, 2) / sqrt(num_repeats);
        end
        FPR_mean = mean(FPR_matrix, 2);

        % Plot the mean curve
        if strcmp(model_name, 'MAGGIC')
            plot(FPR_mean, TPR_mean, 'LineWidth', 2, 'Color', colors(m,:), 'LineStyle', '--');
        else
            plot(FPR_mean, TPR_mean, 'LineWidth', 2, 'Color', colors(m,:), 'LineStyle', '-');
        end

        % Calculate area
        area_mean = trapz(FPR_mean, TPR_mean);
        if strcmp(model_name, 'MAGGIC')
            area_results.(field_name).MeanArea = area_mean;
        else
            area_lower = trapz(FPR_mean, TPR_mean - 1.96 * TPR_se);
            area_upper = trapz(FPR_mean, TPR_mean + 1.96 * TPR_se);
            area_results.(field_name).MeanArea = area_mean;
            area_results.(field_name).LowerArea = area_lower;
            area_results.(field_name).UpperArea = area_upper;

            % Plot the standard error region
            fill([FPR_mean; flipud(FPR_mean)], [TPR_mean - 1.96 * TPR_se; flipud(TPR_mean + 1.96 * TPR_se)], ...
                colors(m,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end

    % Set axis labels
    xlabel('False Positive Rate (1 - Specificity)');
    ylabel('True Positive Rate (Sensitivity)');

    % Set axis limits and diagonal reference line
    axis([0 1 0 1]);
    plot([0 1], [0 1], 'k--'); % Add diagonal reference line
    legend off
    hold off;
    pbaspect([1,1,1]);

    % Remove background and box
    set(gca, 'Color', 'none');

    % Remove title
    title('');

    % Remove axis labels
    xlabel('');
    ylabel('');

    % Adjust y-axis ticks to intervals of 0.2
    yticks(0:0.2:1);

    % Remove grid lines
    grid on;
    box on;

    % Remove x-axis and y-axis tick labels
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
end

% Calculate the AUC of each repetition
results = struct();

for t = 1:length(timePoints)
    time_point = timePoints(t);
    time_point_str = sprintf('Time_%0.1f', time_point);
    time_point_str = strrep(time_point_str, '.', '_'); % Replace dot with underscore

    % Initialize the results table for this time point
    results.(time_point_str) = zeros(num_repeats, length(models));
    col_idx = 1;

    for m = 1:length(models)
        model_name = models{m};

        % Ensure the field name is valid
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_');

        % Retrieve FPR and TPR matrices
        FPR_matrix = data_struct.(field_name).FPR;
        TPR_matrix = data_struct.(field_name).TPR;

        % Calculate the AUC for each column (repetition)
        if strcmp(model_name, 'MAGGIC')
            auc = trapz(FPR_matrix, TPR_matrix);  % Calculate AUC for a single case
            results.(time_point_str)(1, col_idx) = auc;
        else
            for i = 1:num_repeats
                FPR_col = FPR_matrix(:, i);
                TPR_col = TPR_matrix(:, i);
                auc = trapz(FPR_col, TPR_col);  % Calculate area
                results.(time_point_str)(i, col_idx) = auc;
            end
        end

        col_idx = col_idx + 1;
    end
end

% Create result tables and write to an Excel file
output_filepath = 'AUC_results.xlsx';
for t = 1:length(timePoints)
    time_point = timePoints(t);
    time_point_str = sprintf('Time_%0.1f', time_point);
    time_point_str = strrep(time_point_str, '.', '_');

    % Retrieve the result matrix
    auc_matrix = results.(time_point_str);

    % Create a table
    T = array2table(auc_matrix, 'VariableNames', models);

    % Write to an Excel file
    writetable(T, output_filepath, 'Sheet', time_point_str);
end

disp('AUC results have been written to Excel file.');
