clear
T = readtable("PredictorsFinal.xlsx",'VariableNamingRule','preserve');
x = "EAr_D";  % 你可以更改这个变量来选择你想要的列

% 找到Risk列的位置
riskColumn = T.Risk;

% 创建新的列名
columnNames = {'x_Risk1', 'x_Risk2', 'x_Risk3'};

% 初始化新的表格
T_new = table();

% 获取每个Risk值对应的数据
for i = 1:3
    % 对应Risk值等于i时的x值，其他位置填入NaN
    x_values = T.(x)(riskColumn == i);
    
    % 如果长度不够，则填充NaN
    if length(x_values) < height(T)
        x_values(end+1:height(T)) = NaN;
    end
    
    % 添加数据到新表格
    T_new.(columnNames{i}) = x_values;
end

% 显示生成的新表格
disp(T_new);

