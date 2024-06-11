function txt = myupdatefcn(~, event_obj, dataTbl, theta_final, rho_final)
    % Get the cursor position on the plot
    pos = get(event_obj, 'Position');
    [~, thetaIndx] = min(abs(theta_final - pos(1)));
    [~, rhoIndx] = min(abs(rho_final - pos(2)));

    % Check if the indices match to verify it's the correct point
    if thetaIndx == rhoIndx
        index = thetaIndx;
    else
        index = [];  % No exact match found
    end

    % Generate the text to display
    if ~isempty(index)
        txt = {[dataTbl.Properties.RowNames{index}]};
    else
        txt = 'No match found';
    end
end