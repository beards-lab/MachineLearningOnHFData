
function bmuIndex = findBMU(inputVec, somWeights)
% inputVec - Input vector, must be a column vector
% somWeights - SOM neuron weights matrix (each row represents the weight vector of a neuron)

% Initialize the minimum distance and neuron index
minDist = inf;
bmuIndex = -1;

% Iterate over all neurons to compute the distance
for i = 1:size(somWeights, 1)
    % Compute the Euclidean distance between the current neuron weights and the input vector
    dist = norm(somWeights(i, :) - inputVec');

    % If a smaller distance is found, update the minimum distance and neuron index
    if dist < minDist
        minDist = dist;
        bmuIndex = i;
    end
end
end
