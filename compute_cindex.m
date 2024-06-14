function c_index = compute_cindex(times, censoring, risk_scores)
    n = length(times);
    num_correct = 0; % 初始化协调数
    num_tied = 0; % 初始化平局数
    num_possible = 0; % 初始化可能配对数
      
    for i = 1:n
        for j = i+1:n
            % 检查两个患者至少有一个发生了事件
            if censoring(i) == 0 || censoring(j) == 0
                if censoring(i) == 0 && censoring(j) == 0
                    % 若两者都发生了事件
                    num_possible = num_possible + 1;
                    if risk_scores(i) < risk_scores(j) && times(i) > times(j) || ...
                       risk_scores(i) > risk_scores(j) && times(i) < times(j)
                        % 若风险得分与生存时间协调
                        num_correct = num_correct + 1;
                    elseif risk_scores(i) == risk_scores(j)
                        % 若风险得分一致
                        num_tied = num_tied + 1;
                    end
                elseif censoring(i) == 0 && censoring(j) == 1 && times(j) >= times(i) || ...
                       censoring(i) == 1 && censoring(j) == 0 && times(i) >= times(j)
                    % 若一者发生事件（假设为患者i），另一者删失（假设为患者j），且患者j的删失时间大于等于患者i的事件时间
                    num_possible = num_possible + 1;
                    if risk_scores(i) > risk_scores(j) && censoring(i) == 0 || ...
                       risk_scores(i) < risk_scores(j) && censoring(j) == 0
                        % 若风险得分高的患者发生了事件
                        num_correct = num_correct + 1;
                    elseif risk_scores(i) == risk_scores(j)
                        % 若风险得分一致
                        num_tied = num_tied + 1;
                    end
                end
            end
        end
    end
      
    c_index = (num_correct + 0.5 * num_tied) / num_possible; % 计算C指数
end