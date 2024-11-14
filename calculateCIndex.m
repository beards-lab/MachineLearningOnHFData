function [cindex, paired_matrix] = calculateCIndex(predicted_scores, survival_times, status)
    n = length(survival_times);
    nconc = 0;  % 协调对数量
    nrel = 0;   % 有效对数量
    nuncer = 0; % 不确定对数量

    % 初始化配对矩阵
    paired_matrix = zeros(n, n);

    for i = 1:n
        for j = 1:n
            % 避免自比较
            if i ~= j
                dx = predicted_scores(i) - predicted_scores(j);
                dy = survival_times(i) - survival_times(j);

                if (status(i) && (dy < 0)) || (status(i) && ~status(j) && (dy == 0))
                    nrel = nrel + 1;
                    paired_matrix(i, j) = 1; % 计入配对
                    nconc = nconc + (dx > 0) + 0.5 * (dx == 0);
                elseif (status(j) && (dy > 0)) || (status(j) && ~status(i) && (dy == 0))
                    nrel = nrel + 1;
                    paired_matrix(i, j) = 1; % 计入配对
                    nconc = nconc + (dx < 0) + 0.5 * (dx == 0);
                else
                    if ~(status(i) && status(j))
                        nuncer = nuncer + 1;
                    end
                end
            end
        end
    end

    % 计算C-index
    cindex = nconc / nrel;

    % 确保C-index在0到1之间
    if isnan(cindex)
        cindex = 0;  % 无可用对时返回0
    end
end