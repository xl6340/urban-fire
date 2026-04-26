%% Fig. 3B: z-score dot plot, split into Forest and ShrubGrass subplots
clear; clc;

varNames  = {'ignition','elevation','slope','ndviM','vs','ppt','tmax','tmin','tmean','vpdmax'};
xLabels   = {'%Human','Elevation','Slope','NDVI','Wind','P','T_{max}','T_{min}','T_{mean}','VPD_{max}'};

catColors = containers.Map(...
    {'vpdmax','tmean','tmin','tmax','ppt','vs','ndviM','slope','elevation','ignition'}, ...
    {[246,198,175]/255, [246,198,175]/255, [246,198,175]/255, ...
     [246,198,175]/255, [246,198,175]/255, [246,198,175]/255, ...
     [181,212,190]/255, ...
     [175,212,227]/255, [175,212,227]/255, ...
     [184,185,210]/255});

nVars    = numel(varNames);
lc       = {'Forest','ShrubGrass'};
fireType = {'Urban-edge','Wildland'};

% Load data
Data.Forest.Urban = table(); Data.Forest.Wild = table();
Data.Shrub.Urban  = table(); Data.Shrub.Wild  = table();
for i = 1:numel(lc)
    for j = 1:numel(fireType)
        tmp = readtable(sprintf('dataFig/variable/%s_%s.csv', lc{i}, fireType{j}));
        if strcmp(lc{i},'Forest') && strcmp(fireType{j},'Urban-edge')
            Data.Forest.Urban = [Data.Forest.Urban; tmp];
        elseif strcmp(lc{i},'Forest') && strcmp(fireType{j},'Wildland')
            Data.Forest.Wild  = [Data.Forest.Wild;  tmp];
        elseif strcmp(lc{i},'ShrubGrass') && strcmp(fireType{j},'Urban-edge')
            Data.Shrub.Urban  = [Data.Shrub.Urban;  tmp];
        elseif strcmp(lc{i},'ShrubGrass') && strcmp(fireType{j},'Wildland')
            Data.Shrub.Wild   = [Data.Shrub.Wild;   tmp];
        end
    end
end

% Compute z-scores and p-values
stats = zeros(nVars, 4);
pVals = zeros(nVars, 2);
for i = 1:nVars
    vName = varNames{i};
    for k = 1:2
        if k == 1
            uVec = Data.Forest.Urban.(vName); uVec = uVec(~isnan(uVec));
            wVec = Data.Forest.Wild.(vName);  wVec = wVec(~isnan(wVec));
        else
            uVec = Data.Shrub.Urban.(vName);  uVec = uVec(~isnan(uVec));
            wVec = Data.Shrub.Wild.(vName);   wVec = wVec(~isnan(wVec));
        end
        mu = mean(wVec); sigma = std(wVec);
        if sigma == 0, sigma = 1; end
        zScores = (uVec - mu) / sigma;
        stats(i, (k-1)*2+1) = mean(zScores);
        stats(i, (k-1)*2+2) = std(zScores) / sqrt(length(zScores));
        [~, p] = ttest2(uVec, wVec, 'Vartype','unequal');
        pVals(i, k) = p;
    end
end

% Sort by mean z-score descending
[~, idxForest] = sort(stats(:,1), 'descend');
[~, idxShrub]  = sort(stats(:,3), 'descend');

% Plotting
figure('Position', [100, 100, 800, 420]);
ax(1) = axes('Position', [0.08, 0.32, 0.40, 0.61]);
ax(2) = axes('Position', [0.50, 0.32, 0.40, 0.61]);

titles   = {'Forest', 'Shrub & Grassland'};
sortIdx  = {idxForest, idxShrub};
ecoIdx   = {1, 2};

for panel = 1:2
    axes(ax(panel));
    hold on;

    idx      = sortIdx{panel};
    eco      = ecoIdx{panel};
    labs     = xLabels(idx);
    vars     = varNames(idx);
    meanVals = stats(idx, (eco-1)*2+1);
    seVals   = stats(idx, (eco-1)*2+2);
    pV       = pVals(idx, eco);

    % Light background alternating columns for readability
    for i = 1:2:nVars
        patch([i-0.5, i+0.5, i+0.5, i-0.5], [-4, -4, 4, 4], ...
            [0.96 0.96 0.96], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    end

    % Zero line
    yline(0, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);

    % Dot + CI for each variable
    for i = 1:nVars
        clr = catColors(vars{i});
        ci  = 1.96 * seVals(i);

        % Vertical CI line
        plot([i i], [meanVals(i)-ci, meanVals(i)+ci], '-', ...
            'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);

        % Cap lines
        plot([i-0.15, i+0.15], [meanVals(i)+ci, meanVals(i)+ci], '-', ...
            'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);
        plot([i-0.15, i+0.15], [meanVals(i)-ci, meanVals(i)-ci], '-', ...
            'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);

        % Dot
        scatter(i, meanVals(i), 60, clr, 'filled', ...
            'MarkerEdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
    end

    % Significance stars
    for j = 1:nVars
        starStr = '';
        if     pV(j) < 0.001, starStr = '***';
        elseif pV(j) < 0.01,  starStr = '**';
        elseif pV(j) < 0.05,  starStr = '*';
        end
        if ~isempty(starStr)
            ci = 1.96 * seVals(j);
            if meanVals(j) >= 0
                y_pos = meanVals(j) + ci + 0.06;
                va = 'bottom';
            else
                y_pos = meanVals(j) - ci - 0.06;
                va = 'top';
            end
            text(j, y_pos, starStr, 'FontSize', 10, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', va);
        end
    end

    % Axes formatting
    set(gca, 'XTick', 1:nVars, 'XTickLabel', labs, ...
        'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
    xtickangle(45);
    xlim([0.5, nVars + 0.5]);
    ylim([-1.1, 1.1]);
    title(titles{panel}, 'FontSize', 12, 'FontWeight', 'bold');

    if panel == 1
        yticks(-0.9:0.3:0.9);
        ylabel('Standardized deviations (Z-score)', 'FontSize', 11);
        catNames = {'Climate','Vegetation','Terrain','Ignition'};
        catClrs  = {[246,198,175]/255, [181,212,190]/255, ...
                    [175,212,227]/255, [184,185,210]/255};
        for ci = 1:4
            scatter(nan, nan, 60, catClrs{ci}, 'filled', ...
                'MarkerEdgeColor', [0.3 0.3 0.3], ...
                'DisplayName', catNames{ci});
        end
    else
        set(gca, 'YTick', [], 'YColor', 'none');
    end

    hold off;
end
set(gcf, 'Color', 'w');


catNames = {'Climate','Vegetation','Terrain','Ignition'};
catClrs  = {[246,198,175]/255, [181,212,190]/255, ...
            [175,212,227]/255, [184,185,210]/255};

ax_leg = axes('Position', [0.15, 0.01, 0.70, 0.08]);
hold on;

nCat  = numel(catNames);
xPos  = linspace(0.15, 0.85, nCat); 

set(ax_leg, 'XLim', [0 1], 'YLim', [0 1], ...
    'Visible', 'off');  
rectangle('Position', [0, 0, 1, 1], ...
    'EdgeColor', [0.75 0.75 0.75], 'LineWidth', 0.8, 'FaceColor', 'none');

for ci = 1:nCat
    scatter(ax_leg, xPos(ci), 0.5, 60, catClrs{ci}, 'filled', ...
        'MarkerEdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.8);

    % Text label
    text(ax_leg, xPos(ci) + 0.04, 0.5, catNames{ci}, ...
        'FontSize', 11, ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left');
end

hold off;

set(gcf, 'Color', 'w');
exportgraphics(gcf, 'Fig/Fig3B.png');
exportgraphics(gcf, 'Fig/Fig3B.pdf');