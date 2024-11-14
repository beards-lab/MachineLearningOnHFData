%% CALCULATE MAGGIC SCORES
clear
T = readtable('Predictors.xlsx',VariableNamingRule='preserve');
% Assuming T is your table with the necessary columns
T.MAGGIC_Score = zeros(height(T), 1); % Initialize the MAGGIC score column

for i = 1:height(T)
    % Calculate individual parameters
    EF = (T.LVEDV_D(i) - T.LVESV_D(i)) / T.LVEDV_D(i);
    BMI = T.BMI(i);
    Age = T.Age(i);
    SBP = T.SBP_D(i);
    Creat = T.Creatinine(i);
    Sex = T.Gender(i);
    Smoking = strcmp(T.Smoking(i), 'Current');
    Diabetic = T.DM(i);
    COPD = T.COPD(i);
    BB = T.("Usage of beta-blocker")(i);
    ACE = T.("Usage of ACEI")(i);
    ARB = T.("Usage of ARB")(i);
    ARNi = T.("Usage of ARNi")(i);
    NYHA = T.("NYHA class")(i);
    First_Diagnosis = T.("First diagnosis of heart failure in the past 18 months")(i);

    % Initialize score
    score = 0;
    
    % Ejection Fraction
    if EF < 0.2
        score = score + 7;
    elseif EF < 0.25
        score = score + 6;
    elseif EF < 0.30
        score = score + 5;
    elseif EF < 0.35
        score = score + 3;
    elseif EF < 0.40
        score = score + 2;
    end
    
    % Age
    if Age >= 55
        if EF < 0.30
            if Age < 60
                score = score + 1;
            elseif Age < 65
                score = score + 2;
            elseif Age < 70
                score = score + 4;
            elseif Age < 75
                score = score + 6;
            elseif Age < 80
                score = score + 8;
            else
                score = score + 10;
            end
        elseif EF < 0.40
            if Age < 60
                score = score + 2;
            elseif Age < 65
                score = score + 4;
            elseif Age < 70
                score = score + 6;
            elseif Age < 75
                score = score + 8;
            elseif Age < 80
                score = score + 10;
            else
                score = score + 12;
            end
        else
            if Age < 60
                score = score + 3;
            elseif Age < 65
                score = score + 5;
            elseif Age < 70
                score = score + 7;
            elseif Age < 75
                score = score + 9;
            elseif Age < 80
                score = score + 12;
            else
                score = score + 15;
            end
        end
    end
    % Systolic Blood Pressure (SBP)
    if SBP < 110
        if EF < 0.30
            score = score + 5;
        elseif EF < 0.40
            score = score + 3;
        else
            score = score + 2;
        end
    elseif SBP < 120
        if EF < 0.30
            score = score + 4;
        elseif EF < 0.40
            score = score + 1;
        else
            score = score + 1;
        end
    elseif SBP < 130
        if EF < 0.30
            score = score + 3;
        elseif EF < 0.40
            score = score + 1;
        else
            score = score + 1;
        end
    elseif SBP < 140
        if EF < 0.30
            score = score + 2;
        end
    elseif SBP < 150
        if EF < 0.30
            score = score + 1;
        end
    end
    
    % BMI
    if BMI < 15
        score = score + 6;
    elseif BMI < 20
        score = score + 5;
    elseif BMI < 25
        score = score + 3;
    elseif BMI < 30
        score = score + 2;
    end
    
    % Creatinine
    if Creat >= 90
        if Creat < 110
            score = score + 1;
        elseif Creat < 130
            score = score + 2;
        elseif Creat < 150
            score = score + 3;
        elseif Creat < 170
            score = score + 4;
        elseif Creat < 210
            score = score + 5;
        elseif Creat < 250
            score = score + 6;
        else
            score = score + 8;
        end
    end
    
    % Gender
    if  Sex == 1
        score = score + 1;
    end
    
    % Smoking
    if Smoking
        score = score + 1;
    end
    
    % Diabetes
    if Diabetic
        score = score + 3;
    end
    
    % COPD
    if COPD
        score = score + 2;
    end
    
    % Beta-blocker
    if ~BB
        score = score + 3;
    end
    
    % ACEI/ARB
    if ~ACE&&~ARB&&~ARNi
        score = score + 1;
    end

    % NYHA Class
    if NYHA == 1
        score = score+0;
    elseif NYHA == 2
        score = score + 2;
    elseif NYHA == 3
        score = score + 6;
    elseif NYHA == 4
        score = score + 8;
    else
        score = score+4;
    end
    
    % First Diagnosis of Heart Failure in the Past 18 Months
    if First_Diagnosis
        score = score + 2;
    end
    
    % Assign the score to the table
    T.MAGGIC_Score(i) = score;
end
