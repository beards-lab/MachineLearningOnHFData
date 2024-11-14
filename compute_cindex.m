function [c_index, paired_matrix] = compute_cindex(times, censoring, risk_scores)
    n = length(times);
    num_correct = 0; % 初始化协调数
    num_tied = 0; % 初始化平局数
    num_possible = 0; % 初始化可能配对数

    % 初始化配对矩阵
    paired_matrix = zeros(n, n);

    for i = 1:n
        for j = 1:n
            % 避免自比较
            if i ~= j
                % 检查至少一个患者有事件
                if censoring(i) == 0 || censoring(j) == 0
                    % 两者都发生事件
                    if times(i) == times(j) && censoring(i) == 0 && censoring(j) == 0
                        continue;
                    end
                    if censoring(i) == 0 && censoring(j) == 0
                        num_possible = num_possible + 1;
                        paired_matrix(i, j) = 1; % 计入配对
                        if (risk_scores(i) < risk_scores(j) && times(i) > times(j)) || ...
                                (risk_scores(i) > risk_scores(j) && times(i) < times(j))
                            num_correct = num_correct + 1;
                        elseif risk_scores(i) == risk_scores(j)
                            num_tied = num_tied + 1;
                        end
                    % 仅一方事件另一方删失
                    elseif censoring(i) == 0 && censoring(j) == 1 && times(j) >= times(i)
                        num_possible = num_possible + 1;
                        paired_matrix(i, j) = 1; % 计入配对
                        if risk_scores(i) > risk_scores(j)
                            num_correct = num_correct + 1;
                        elseif risk_scores(i) == risk_scores(j)
                            num_tied = num_tied + 1;
                        end
                    elseif censoring(i) == 1 && censoring(j) == 0 && times(i) >= times(j)
                        num_possible = num_possible + 1;
                        paired_matrix(i, j) = 1; % 计入配对
                        if risk_scores(i) < risk_scores(j)
                            num_correct = num_correct + 1;
                        elseif risk_scores(i) == risk_scores(j)
                            num_tied = num_tied + 1;
                        end
                    end
                end
            end
        end
    end

    c_index = (num_correct + 0.5 * num_tied) / num_possible; % 计算C指数
end