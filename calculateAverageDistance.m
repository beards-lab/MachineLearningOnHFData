function avgDist = calculateAverageDistance(inputData, somWeights)
    
    totalDist = 0;  % Variable to accumulate distances for all samples
    
    % inputData needs to be transposed because the previous code assumes samples by rows, now by columns
    for i = 1:size(inputData, 2)  % Note we use 2, because inputs are by columns
        inputVec = inputData(:, i);  % Retrieve a single sample column vector
        bmuIndex = findBMU(inputVec, somWeights);  % Note the transpose is included
        
        % Calculate and accumulate the distance
        totalDist = totalDist + norm(inputVec' - somWeights(bmuIndex, :));
    end
    
    % Calculate the average distance from all input samples to their respective BMUs
    avgDist = totalDist / size(inputData, 2);
end