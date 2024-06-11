function visualizeSOMGrid(clusterIndexMap, dimension1, dimension2, K, titleStr,n)
    figure(n); clf; hold on;

    % Define the hexagon radius and spacing
    hexRadius = 0.5;
    hexWidth = hexRadius * sqrt(3);
    horizSpacing = hexWidth;
    vertSpacing = hexRadius * 1.5;
    colorMap = [0.75, 0.75, 0.75; hsv(K)]; % Creating a colormap with gray for noise

    % Drawing the grid
    for r = 1:dimension1
        for c = 1:dimension2
            if mod(r, 2) == 0
                centerPos = [(c + 0.5) * horizSpacing, r * vertSpacing];
            else
                centerPos = [c * horizSpacing, r * vertSpacing];
            end

            angles = (-pi/6) + (0:5) * (2*pi/6);
            hexagonX = centerPos(1) + hexRadius * cos(angles);
            hexagonY = centerPos(2) + hexRadius * sin(angles);
            
            % Get the cluster index; if it's noise (-1), use the first color (gray)
            clusterIdx = clusterIndexMap(r, c);
            if clusterIdx == -1  % For noise points
                colorIdx = 1;
            else
                colorIdx = clusterIdx + 1; % Because we put gray at the first position
            end
            
            % Fill hexagons with the appropriate cluster color
            fill(hexagonX, hexagonY, colorMap(colorIdx, :), 'EdgeColor', 'white');
            text(centerPos(1), centerPos(2), num2str(clusterIdx), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
        end
    end

    axis equal;
    axis off;
    hold off;
    title(titleStr);
end