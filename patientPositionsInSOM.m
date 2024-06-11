function [col, row] = patientPositionsInSOM(net, data)
    % Apply the SOM network to the dataset to find the BMU index for each patient
    bmuIndices = vec2ind(net(data));

    % Calculate the grid dimensions of the SOM
    gridDimensions = net.layers{1}.dimensions;
    numRows = gridDimensions(1);
    numCols = gridDimensions(2);

    % Initialize row and column vectors
    row = zeros(1, size(data, 2));
    col = zeros(1, size(data, 2));
    
    % Convert BMU index to row and column positions
    for i = 1:length(bmuIndices)
        [col(i), row(i)] = ind2sub([numRows, numCols], bmuIndices(i));
    end
end