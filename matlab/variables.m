% ... (Keep the setup and data loading sections from the previous code) ...

% --- 3. PLOTTING ---
figure('Position',[100 50 1100 800], 'Color', 'w');
t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nVars
    ax = nexttile; hold(ax, 'on');
    vName = varNames{i};
    
    % Get Data
    try vecFU = Data.Forest.Urban.(vName); catch, vecFU = []; end
    try vecFW = Data.Forest.Wild.(vName);  catch, vecFW = []; end
    try vecSU = Data.Shrub.Urban.(vName);  catch, vecSU = []; end
    try vecSW = Data.Shrub.Wild.(vName);   catch, vecSW = []; end
    
    % Clean NaNs
    vecFU = vecFU(~isnan(vecFU)); vecFW = vecFW(~isnan(vecFW));
    vecSU = vecSU(~isnan(vecSU)); vecSW = vecSW(~isnan(vecSW));

    % Unit conversions
    if strcmp(vName, 'vpdmax')
        vecFU=vecFU/10; vecFW=vecFW/10; vecSU=vecSU/10; vecSW=vecSW/10;
    end
    
    % --- PLOT BOXES (Whiskers REMOVED) ---
    b1 = boxchart(ax, 1.0*ones(size(vecFW)), vecFW, 'BoxFaceColor', cWild,'Notch', 'on');
    b2 = boxchart(ax, 1.8*ones(size(vecFU)), vecFU, 'BoxFaceColor', cUrban,'Notch', 'on');
    b3 = boxchart(ax, 3.5*ones(size(vecSW)), vecSW, 'BoxFaceColor', cWild,'Notch', 'on');
    b4 = boxchart(ax, 4.3*ones(size(vecSU)), vecSU, 'BoxFaceColor', cUrban,'Notch', 'on');
    
    set([b1 b2 b3 b4], 'MarkerStyle', 'none', 'BoxFaceAlpha', 0.7, ...
        'WhiskerLineColor', 'none', 'BoxWidth', 0.5);

    % --- STATISTICS ---
    % Forest Stats
    dF = mean(vecFU) - mean(vecFW);
    [~, pF] = ttest2(vecFU, vecFW, 'Vartype','unequal');
    if pF < 0.001, pStrF = 'p<0.001'; else, pStrF = sprintf('p=%.2f', pF); end
    
    % Shrub Stats
    dS = mean(vecSU) - mean(vecSW);
    [~, pS] = ttest2(vecSU, vecSW, 'Vartype','unequal');
    if pS < 0.001, pStrS = 'p<0.001'; else, pStrS = sprintf('p=%.2f', pS); end
    
    title(ax, sprintf('F:%+.1f (\\it%s\\rm)  S:%+.1f (\\it%s\\rm)', ...
          dF, pStrF, dS, pStrS), 'FontSize', 8, 'FontWeight', 'normal');

    % --- DYNAMIC Y-LIMITS (TIGHT ZOOM) ---
    % Calculate the 25th (bottom of box) and 75th (top of box) for all 4 groups
    q1_vals = [prctile(vecFW, 25), prctile(vecFU, 25), prctile(vecSW, 25), prctile(vecSU, 25)];
    q3_vals = [prctile(vecFW, 75), prctile(vecFU, 75), prctile(vecSW, 75), prctile(vecSU, 75)];
    
    % Filter out NaNs in case some data is missing
    q1_vals = q1_vals(~isnan(q1_vals));
    q3_vals = q3_vals(~isnan(q3_vals));
    
    if ~isempty(q1_vals) && ~isempty(q3_vals)
        global_min = min(q1_vals);
        global_max = max(q3_vals);
        height_range = global_max - global_min;
        
        if height_range == 0, height_range = 1; end % Prevent flat line error
        
        % Add 10% padding above the highest box and below the lowest box
        padding = 0.15 * height_range; 
        ylim(ax, [global_min - padding, global_max + padding]);
    end
    
    ylabel(ax, yLabels{i}, 'FontSize', 10);
    
    % X-Axis Settings
    xlim(ax, [0 5.5]);
    xticks(ax, [1.4, 3.9]); 
    xticklabels(ax, {'Forest', 'Shrub'});
    if i <= 6, xticklabels(ax, {}); end
    
    % Annotation
    xl = xlim(ax); yl = ylim(ax);
    text(ax, xl(1)+0.02*(xl(2)-xl(1)), yl(2)-0.02*(yl(2)-yl(1)), ...
         labels{i}, 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
     
    box(ax, 'on'); ax.TickDir = 'out';
end

% --- LEGEND ---
lg = legend([b1, b2], {'Wildland', 'Urban-edge'}, 'Orientation', 'horizontal');
lg.Layout.Tile = 'North';

exportgraphics(gcf, 'Fig/variables.png');