%% NDVI calculation, ShrubGrass only
clear;clc;

decades  = {'1990s', '2000s', '2010s', '2020s'};  nDec  = numel(decades);
fireType = {'WUI','Wildland'};             nFire = numel(fireType);

data     = shaperead('dataPrc/firePrmt/CalFire.shp');
ndviAll  = [data.ndvi]';
lc       = {data.lc}';
Dec      = {data.decade}';
Fire     = {data.FireType}';

colors   = {[216, 118, 89]/255; [41, 157, 143]/255};

figure('Position', [100, 100, 400, 500]); 
t = tiledlayout(nDec-1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
binW = 0.02;
edge = 0:binW:1;
xFit = 0:0.005:1; 
for d = 2:nDec
    ax = nexttile;
    hold(ax, 'on');      
    maskLC = strcmp(lc, 'ShrubGrass');
    
    urbanVals = []; wildVals = [];        
    for i = 1:nFire
        mask = maskLC & strcmp(Fire, fireType{i}) & strcmp(Dec, decades{d}) & ~isnan(ndviAll);
        ndvi = ndviAll(mask);
        if i == 1, urbanVals = ndvi; else, wildVals = ndvi; end
        
        mu = mean(ndvi); sigma = std(ndvi);
        histogram(ax, ndvi, 'BinEdges', edge, 'Normalization', 'probability', ...
            'FaceColor', 0.8*colors{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
            'HandleVisibility', 'off');
        
        yFit = normpdf(xFit, mu, sigma) * binW;
        plot(ax, xFit, yFit, 'Color', colors{i}, 'LineWidth', 1.5, ...
            'DisplayName', fireType{i});            
        xline(ax, mu, 'LineWidth', 1.5, 'LineStyle', ':', ...
            'Color', colors{i}, 'HandleVisibility', 'off');
    end        
    
    [~, p_val] = ttest2(urbanVals, wildVals, 'Vartype', 'unequal');
    if p_val < 0.01, pStr = '{\it p} < 0.01'; 
    else, pStr = sprintf('{\\it p} = %.2f', p_val); end
    text(ax, -0.02, 0.5, pStr, 'Units', 'normalized');
    
    xlim(ax, [0.1 0.9]); 
    ylim(ax, [0 0.15]); 
    
    if d == 2, title(ax, 'Shrub & Grassland', 'FontWeight','normal'); end        
    
    text(ax, 0.95, 0.80, decades{d}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right'); 
   
    if d == nDec, xlabel(ax, 'NDVI'); end
    set(ax, 'FontSize', 10, 'TickDir', 'out', 'Box', 'off');
    
    if d == 2
        lgd = legend(ax, 'Location', 'east', 'Box', 'off');
        lgd.ItemTokenSize = [15, 15];
    end
end
ylabel(t, 'Probability');

exportgraphics(gcf, 'figs/Fig5.pdf');




%% NDVI calculation, ShrubGrass only
clear;clc;
decades  = {'1990s', '2000s', '2010s', '2020s'};  nDec  = numel(decades);
fireType = {'WUI','Wildland'};             nFire = numel(fireType);
data     = shaperead('dataPrc/firePrmt/CalFire.shp');
ndviAll  = [data.ndvi]';
lc       = {data.lc}';
Dec      = {data.decade}';
Fire     = {data.FireType}';
colors   = {[216, 118, 89]/255; [41, 157, 143]/255};

maskLC_all = strcmp(lc, 'ShrubGrass') & ~isnan(ndviAll);
threshold_95 = prctile(ndviAll(maskLC_all), 95);

figure('Position', [100, 100, 400, 500]);
t = tiledlayout(nDec-1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
binW = 0.02;
edge = 0:binW:1;
xFit = 0:0.005:1;

for d = 2:nDec
    ax = nexttile;
    hold(ax, 'on');
    maskLC = strcmp(lc, 'ShrubGrass');
    urbanVals = []; wildVals = [];
    mu_vals = zeros(1, nFire);
    pct95_vals = zeros(1, nFire);
    
    for i = 1:nFire
        mask = maskLC & strcmp(Fire, fireType{i}) & strcmp(Dec, decades{d}) & ~isnan(ndviAll);
        ndvi = ndviAll(mask);
        if i == 1, urbanVals = ndvi; else, wildVals = ndvi; end
        mu = mean(ndvi); sigma = std(ndvi);
        mu_vals(i) = mu;
        pct95_vals(i) = prctile(ndvi, 95);
        
        histogram(ax, ndvi, 'BinEdges', edge, 'Normalization', 'probability', ...
            'FaceColor', 0.8*colors{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
            'HandleVisibility', 'off');
        yFit = normpdf(xFit, mu, sigma) * binW;
        plot(ax, xFit, yFit, 'Color', colors{i}, 'LineWidth', 1.5, 'DisplayName', fireType{i});
        xline(ax, mu, 'LineWidth', 1.5, 'LineStyle', ':', ...
            'Color', colors{i}, 'HandleVisibility', 'off');
        xline(ax, pct95_vals(i), 'LineWidth', 1.2, 'LineStyle', '--', ...
            'Color', colors{i}, 'HandleVisibility', 'off');
    end
    
    % Mean gap annotation with arrows
    y_mean = 0.10;
    mean_gap = mu_vals(1) - mu_vals(2);
    x_mid_mean = mean(mu_vals);
    % Draw line between mean lines
    plot(ax, [mu_vals(2), mu_vals(1)], [y_mean, y_mean], '-', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(ax, [mu_vals(2), mu_vals(2)], [y_mean-0.004, y_mean+0.004], '-', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(ax, [mu_vals(1), mu_vals(1)], [y_mean-0.004, y_mean+0.004], '-', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    text(ax, x_mid_mean, y_mean + 0.015, sprintf('%.2f', mean_gap), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'k');
    
    % 95th percentile gap annotation with arrows
    y_pct = 0.10;
    pct95_gap = pct95_vals(1) - pct95_vals(2);
    x_mid_pct = mean(pct95_vals);
    plot(ax, [pct95_vals(2), pct95_vals(1)], [y_pct, y_pct], '-', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(ax, [pct95_vals(2), pct95_vals(2)], [y_pct-0.004, y_pct+0.004], '-', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(ax, [pct95_vals(1), pct95_vals(1)], [y_pct-0.004, y_pct+0.004], '-', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    text(ax, x_mid_pct, y_pct + 0.015, sprintf('%.2f', pct95_gap), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'k');
    
    [~, p_val] = ttest2(urbanVals, wildVals, 'Vartype', 'unequal');
    if p_val < 0.01, pStr = '{\it p} < 0.01';
    else, pStr = sprintf('{\\it p} = %.2f', p_val); end
    text(ax, -0.02, 0.5, pStr, 'Units', 'normalized', 'Color', 'k');
    
    xlim(ax, [0.1 0.9]);
    ylim(ax, [0 0.15]);
    if d == 2, title(ax, 'Shrub & Grassland', 'FontWeight', 'normal'); end
    text(ax, 0.95, 0.80, decades{d}, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'Color', 'k');
    if d == nDec, xlabel(ax, 'NDVI'); end
    set(ax, 'FontSize', 10, 'TickDir', 'out', 'Box', 'off');
    if d == 2
        lgd = legend(ax, {'WUI', 'Wildland'}, 'Location', 'east', 'Box', 'off');
        lgd.ItemTokenSize = [15, 15];
    end
end

ylabel(t, 'Probability');
exportgraphics(gcf, 'figs/Fig5.pdf');
exportgraphics(gcf, 'figs/Fig5.png');