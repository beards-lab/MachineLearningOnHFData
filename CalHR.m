clear;
clc;

% 文件名列表
file_names = {'DM_roc_results.xlsx', 'D_roc_results.xlsx', 'M_roc_results.xlsx', 'MAGGIC_roc_results.xlsx', 'MAGGIC_M_roc_results.xlsx'};
model_names = {'DM', 'D', 'M', 'MAGGIC', 'MAGGIC_with_M'};

% 初始化结构体存储所有数据
roc_data_struct = struct();

% 读取数据
for file_idx = 1:length(file_names)
    file_name = file_names{file_idx};
    model_name = model_names{file_idx};

    % 读取附表名称
    [~, sheet_names] = xlsfinfo(file_name);

    for sheet_idx = 1:length(sheet_names)
        sheet_name = sheet_names{sheet_idx};
        % 先获取 sheet_name 最后一个下划线之前的部分
        last_underscore_idx = find(sheet_name == '_', 1, 'last');
        sheet_name_prefix = sheet_name(1:last_underscore_idx-1);

        % 生成正则表达式来匹配 sheet_name 中的时间点
        tokens = regexp(sheet_name, sprintf('%s_(\\d+\\.?\\d*)', sheet_name_prefix), 'tokens');

        time_point = str2double(tokens{1}{1});

        % 读取当前附表的数据
        data = readtable(file_name, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');

        % 确保字段名称合法
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_'); % 替换点号

        % 将数据保存到结构体中
        roc_data_struct.(field_name).FPR = data.FPR;
        roc_data_struct.(field_name).TPR = data.TPR;
        roc_data_struct.(field_name).Model = model_name;
        roc_data_struct.(field_name).TimePoint = time_point;
    end
end
%%
% 获取唯一的时间点和模型名称 
timePoints = [0.5, 1];
models = model_names;

% 设置每个时间点每个模型的重复次数
num_repeats = 25;

% 获取当前数据点数量
num_points = length(roc_data_struct.(['Model_', models{1}, '_Time_0_5']).FPR) / num_repeats;

% 设定新的数据点数量
new_points = 350;

% 生成新的 FPR 向量
new_FPR = linspace(0, 1, new_points);

% 初始化结构体
data_struct = struct();

% 按时间点和模型分离数据，并进行插值
for t = 1:length(timePoints)
    time_point = timePoints(t);

    for m = 1:length(models)
        model_name = models{m};

        % 确保字段名称合法
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_'); % 替换点号

        % 获取 FPR 和 TPR 数据并重塑为矩阵
        currentModelDataFPR = roc_data_struct.(field_name).FPR;
        currentModelDataTPR = roc_data_struct.(field_name).TPR;

        if strcmp(model_name, 'MAGGIC')
            num_points = length(currentModelDataFPR); % 对于MAGGIC，使用全部数据点
            num_repeats_corrected = 1; % MAGGIC没有重复
        else
            % 确认重复次数对于其他模型
            num_repeats_corrected = num_repeats;
            num_points = length(currentModelDataFPR) / num_repeats_corrected;
            currentModelDataFPR = reshape(currentModelDataFPR, num_points, num_repeats_corrected);
            currentModelDataTPR = reshape(currentModelDataTPR, num_points, num_repeats_corrected);
        end

        % 初始化新的 TPR 矩阵
        new_TPR_matrix = zeros(new_points, num_repeats_corrected);

        % 对每一列进行插值
        for i = 1:num_repeats_corrected
            if num_repeats_corrected == 1
                FPR_col = currentModelDataFPR;
                TPR_col = currentModelDataTPR;
            else
                FPR_col = currentModelDataFPR(:, i);
                TPR_col = currentModelDataTPR(:, i);
            end

            % 检查并去除重复的 FPR 值
            [FPR_unique, unique_idx] = unique(FPR_col);
            TPR_unique = TPR_col(unique_idx);

            % 在插值前强制添加起点 (0, 0) 和终点 (1, 1)
            if FPR_unique(1) ~= 0
                FPR_unique = [0; FPR_unique];
                TPR_unique = [0; TPR_unique];
            end
            if FPR_unique(end) ~= 1
                FPR_unique = [FPR_unique; 1];
                TPR_unique = [TPR_unique; 1];
            end

            % 插值
            new_TPR_matrix(:, i) = interp1(FPR_unique, TPR_unique, new_FPR, 'linear');
            new_TPR_matrix(1, i) = 0;
        end

        % 将新的 FPR 向量重复 num_repeats_corrected 次
        new_FPR_matrix = repmat(new_FPR', 1, num_repeats_corrected);

        % 将新的 FPR 向量和插值后的 TPR 矩阵保存到结构体中
        data_struct.(field_name).FPR = new_FPR_matrix;
        data_struct.(field_name).TPR = new_TPR_matrix;
        data_struct.(field_name).Model = model_name;
        data_struct.(field_name).TimePoint = time_point;
    end
end


%% 可视化每个时间点的数据
colors = [
    119, 172, 48;
    237, 176, 33;
    160, 0, 0;
    0, 0, 0;
     0, 155, 189;
    ] / 255;

for t = 1:length(timePoints)
    time_point = timePoints(t);

    % 创建新的图形窗口
    figure(t); clf;
    hold on;

    % 遍历所有模型
    for m = length(models)-2:-1:1
        model_name = models{m};

        % 确保字段名称合法
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_');

        % 获取 FPR 和 TPR 矩阵
        FPR_matrix = data_struct.(field_name).FPR;
        TPR_matrix = data_struct.(field_name).TPR;

        % 计算均值和标准误
        if strcmp(model_name, 'MAGGIC')
            TPR_mean = TPR_matrix;
            TPR_se = zeros(length(TPR_mean), 1); % 没有重复，标准误为零
        else
            TPR_mean = mean(TPR_matrix, 2);
            TPR_se = std(TPR_matrix, 0, 2) / sqrt(num_repeats);
        end
        FPR_mean = mean(FPR_matrix, 2);
        % 绘制均值线
        if strcmp(model_name, 'MAGGIC')
            plot(FPR_mean, TPR_mean, 'LineWidth', 2, 'Color', colors(m,:),'LineStyle','--');
        else
            plot(FPR_mean, TPR_mean, 'LineWidth', 2, 'Color', colors(m,:),'LineStyle','-');
        end

        % 计算面积
        area_mean = trapz(FPR_mean, TPR_mean);
        if strcmp(model_name, 'MAGGIC')
            area_results.(field_name).MeanArea = area_mean;
        else
            area_lower = trapz(FPR_mean, TPR_mean - 1.96 * TPR_se);
            area_upper = trapz(FPR_mean, TPR_mean + 1.96 * TPR_se);
            area_results.(field_name).MeanArea = area_mean;
            area_results.(field_name).LowerArea = area_lower;
            area_results.(field_name).UpperArea = area_upper;

            % 绘制标准误区域
            fill([FPR_mean; flipud(FPR_mean)], [TPR_mean - 1.96 * TPR_se; flipud(TPR_mean + 1.96 * TPR_se)], ...
                colors(m,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end

    % 设置轴标签
    xlabel('False Positive Rate (1 - Specificity)');
    ylabel('True Positive Rate (Sensitivity)');

    % 设置轴比例和对角线
    axis([0 1 0 1]);
    plot([0 1], [0 1], 'k--'); % 添加对角线作为参考
    legend off
    hold off;
    pbaspect([1,1,1])

    % 去掉背景和框线
    set(gca, 'Color', 'none');

    % 移除标题
    title('');

    % 移除坐标轴标签
    xlabel('');
    ylabel('');

    % 调整y轴刻度，每0.2一个刻度
    yticks(0:0.2:1);

    % 移除网格线
    grid on;
    box on;

    % 移除X轴和Y轴的刻度标签
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
end

% 计算每次重复的曲线面积（AUC）
results = struct();

for t = 1:length(timePoints)
    time_point = timePoints(t);
    time_point_str = sprintf('Time_%0.1f', time_point);
    time_point_str = strrep(time_point_str, '.', '_'); % 替换点号

    % 初始化该时间点的结果表格
    results.(time_point_str) = zeros(num_repeats, length(models));
    col_idx = 1;

    for m = 1:length(models)
        model_name = models{m};

        % 确保字段名称合法
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_');

        % 获取 FPR 和 TPR 矩阵
        FPR_matrix = data_struct.(field_name).FPR;
        TPR_matrix = data_struct.(field_name).TPR;

        % 对每一列（每次重复）计算面积（AUC）
        if strcmp(model_name, 'MAGGIC')
            auc = trapz(FPR_matrix, TPR_matrix);  % 计算单次AUC
            results.(time_point_str)(1, col_idx) = auc;
        else
            for i = 1:num_repeats
                FPR_col = FPR_matrix(:, i);
                TPR_col = TPR_matrix(:, i);
                auc = trapz(FPR_col, TPR_col);  % 计算面积
                results.(time_point_str)(i, col_idx) = auc;
            end
        end

        col_idx = col_idx + 1;
    end
end

% 创建结果表格并写入 Excel 文件
output_filepath = 'AUC_results.xlsx';
for t = 1:length(timePoints)
    time_point = timePoints(t);
    time_point_str = sprintf('Time_%0.1f', time_point);
    time_point_str = strrep(time_point_str, '.', '_');

    % 获取结果矩阵
    auc_matrix = results.(time_point_str);

    % 创建表格
    T = array2table(auc_matrix, 'VariableNames', models);

    % 写入 Excel 文件
    writetable(T, output_filepath, 'Sheet', time_point_str);
end

disp('AUC results have been written to Excel file.');
%%
for t = 1:length(timePoints)
    time_point = timePoints(t);

    % 创建新的图形窗口
    figure(t); clf;
    hold on;

    % 遍历所有模型
    for m = 1:3
        model_name = models{m};

        % 确保字段名称合法
        field_name = sprintf('Model_%s_Time_%0.1f', model_name, time_point);
        field_name = strrep(field_name, '.', '_');

        % 获取 FPR 和 TPR 矩阵
        FPR_matrix = data_struct.(field_name).FPR;
        TPR_matrix = data_struct.(field_name).TPR;

        % 计算均值和标准误
        if strcmp(model_name, 'MAGGIC')
            TPR_mean = TPR_matrix;
            TPR_se = zeros(length(TPR_mean), 1); % 没有重复，标准误为零
        else
            TPR_mean = mean(TPR_matrix, 2);
            TPR_se = std(TPR_matrix, 0, 2) / sqrt(num_repeats);
        end
        FPR_mean = mean(FPR_matrix, 2);
        % 绘制均值线
        if strcmp(model_name, 'MAGGIC')
            plot(FPR_mean, TPR_mean, 'LineWidth', 2, 'Color', colors(m,:),'LineStyle','--');
        else
            plot(FPR_mean, TPR_mean, 'LineWidth', 2, 'Color', colors(m,:),'LineStyle','-');
        end

        % 计算面积
        area_mean = trapz(FPR_mean, TPR_mean);
        if strcmp(model_name, 'MAGGIC')
            area_results.(field_name).MeanArea = area_mean;
        else
            area_lower = trapz(FPR_mean, TPR_mean - 1.96 * TPR_se);
            area_upper = trapz(FPR_mean, TPR_mean + 1.96 * TPR_se);
            area_results.(field_name).MeanArea = area_mean;
            area_results.(field_name).LowerArea = area_lower;
            area_results.(field_name).UpperArea = area_upper;

            % 绘制标准误区域
            fill([FPR_mean; flipud(FPR_mean)], [TPR_mean - 1.96 * TPR_se; flipud(TPR_mean + 1.96 * TPR_se)], ...
                colors(m,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end

    % 设置轴标签
    xlabel('False Positive Rate (1 - Specificity)');
    ylabel('True Positive Rate (Sensitivity)');

    % 设置轴比例和对角线
    axis([0 1 0 1]);
    plot([0 1], [0 1], 'k--'); % 添加对角线作为参考
    legend off
    hold off;
    pbaspect([1,1,1])

    % 去掉背景和框线
    set(gca, 'Color', 'none');

    % 移除标题
    title('');

    % 移除坐标轴标签
    xlabel('');
    ylabel('');

    % 调整y轴刻度，每0.2一个刻度
    yticks(0:0.2:1);

    % 移除网格线
    grid on;
    box on;

    % 移除X轴和Y轴的刻度标签
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
end
