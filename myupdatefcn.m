function txt = myupdatefcn(~, event_obj, dataTbl, angles)
    % Get the current cursor position
    pos = get(event_obj, 'Position');
    x = pos(1);
    y = pos(2);
    z = pos(3);

    % Data point coordinates
    xData = sum(cos(angles).*table2array(dataTbl(:, {'RD12', 'RD23', 'RD31'})), 2, 'omitnan');
    yData = sum(sin(angles).*table2array(dataTbl(:, {'RD12', 'RD23', 'RD31'})), 2, 'omitnan');
    zData = -log10(table2array(dataTbl(:, 'P')));

    % Calculate the distance to the clicked point
    distances = sqrt((xData - x).^2 + (yData - y).^2 + (zData - z).^2);

    % Find the closest point
    [minDistance, index] = min(distances);
    threshold = 0.01; % Distance threshold for determining match

    % If the minimum distance is less than the threshold, consider it a match
    if minDistance < threshold
        txt = dataTbl.Properties.RowNames{index};
    else
        txt = 'No match found';
    end
end