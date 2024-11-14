function drawColoredRing(cartesianAxes, sectorangles, colororder, maxRadius, width)
    for i = 1:length(sectorangles)-1
        theta_span = linspace(sectorangles(i), sectorangles(i+1), 100); % More points for a smooth arc
        inner_radius = maxRadius + 0.1; % Adjusted inner radius of the ring
        outer_radius = maxRadius + width; % Outer radius of the ring

        % Convert polar coordinates to Cartesian coordinates
        [x_inner, y_inner] = pol2cart(theta_span, inner_radius * ones(size(theta_span)));
        [x_outer, y_outer] = pol2cart(theta_span, outer_radius * ones(size(theta_span)));

        % Merge into polygon coordinates
        x_coords = [x_inner, fliplr(x_outer)];
        y_coords = [y_inner, fliplr(y_outer)];

        % Determine the correct color index as per your sector angle
        if i >= 6
            colorIdx = 12 - i + 6;
        else
            colorIdx = 1 - i + 5;
        end

        % Plot a filled polygon in the Cartesian coordinate system
        patch(x_coords, y_coords, colororder(colorIdx, :), 'FaceAlpha', 0.75, 'EdgeColor', 'none', 'Parent', cartesianAxes);
    end
end