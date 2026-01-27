%% Fig. 5A: z-score
clear; clc;

varNames  = {'ignition','elevation','slope', 'ndviM', 'vs', 'ppt', 'tmax', 'tmin', 'tmean', 'vpdmax'};
xLabels   = {'Human','Elevation', 'Slope', 'NDVI', 'Wind', 'P', 'T_{max}', 'T_{min}', 'T_{mean}', 'VPD_{max}'};

nVars     = numel(varNames);
lc        = {'Forest', 'ShrubGrass'};
fireType  = {'Urban-edge','Wildland'};

cForest = [50, 100, 50]/255;    % Dark Green (Forest)
cShrub  = [100, 180, 100]/255;  % Light Green (Shrub)

Data.Forest.Urban = table(); Data.Forest.Wild  = table();
Data.Shrub.Urban  = table(); Data.Shrub.Wild   = table();
for i = 1:numel(lc)
    for j = 1:numel(fireType)
        tmp = readtable(sprintf('dataFig/variable/%s_%s.csv', lc{i}, fireType{j}));
        if strcmp(lc{i}, 'Forest') && strcmp(fireType{j}, 'Urban-edge')
            Data.Forest.Urban = [Data.Forest.Urban; tmp];
        elseif strcmp(lc{i}, 'Forest') && strcmp(fireType{j}, 'Wildland')
            Data.Forest.Wild = [Data.Forest.Wild; tmp];
        elseif strcmp(lc{i}, 'ShrubGrass') && strcmp(fireType{j}, 'Urban-edge')
            Data.Shrub.Urban = [Data.Shrub.Urban; tmp];
        elseif strcmp(lc{i}, 'ShrubGrass') && strcmp(fireType{j}, 'Wildland')
            Data.Shrub.Wild = [Data.Shrub.Wild; tmp];
        end
    end
end

stats = zeros(nVars, 4); 
pVals = zeros(nVars, 2); 
for i = 1:nVars
    vName = varNames{i};
    % Forest
    uVec = Data.Forest.Urban.(vName); uVec = uVec(~isnan(uVec));
    wVec = Data.Forest.Wild.(vName);  wVec = wVec(~isnan(wVec));
    mu = mean(wVec); sigma = std(wVec);
    % Safety check for constant variables (avoid div by zero)
    if sigma == 0, sigma = 1; end 
    zScores = (uVec - mu) / sigma;    
    stats(i, 1) = mean(zScores);
    stats(i, 2) = std(zScores) / sqrt(length(zScores));
    [~, p] = ttest2(uVec, wVec, 'Vartype','unequal');
    pVals(i, 1) = p;
    
    % Shrub
    uVec = Data.Shrub.Urban.(vName); uVec = uVec(~isnan(uVec));
    wVec = Data.Shrub.Wild.(vName);  wVec = wVec(~isnan(wVec));    
    mu = mean(wVec); sigma = std(wVec);    
    if sigma == 0, sigma = 1; end
    zScores = (uVec - mu) / sigma;    
    stats(i, 3) = mean(zScores);
    stats(i, 4) = std(zScores) / sqrt(length(zScores));    
    [~, p] = ttest2(uVec, wVec, 'Vartype','unequal');
    pVals(i, 2) = p;
end

% Plotting 
figure('Position', [100, 100, 500, 500]); 
hold on;
for i = 0:nVars
    yline(i + 0.5, '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 1);
end
b = barh(stats(:, [1, 3]), 1, 'grouped');
b(1).FaceColor = cForest; b(1).EdgeColor = 'none';
b(2).FaceColor = cShrub;  b(2).EdgeColor = 'none';
ngroups = nVars;
nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    y = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    meanVal = stats(:, (i-1)*2 + 1);
    seVal   = stats(:, (i-1)*2 + 2);
    pVal    = pVals(:, i);
    
    errorbar(meanVal, y, seVal, 'horizontal', 'k', ...
        'linestyle', 'none', 'LineWidth', 1, 'CapSize', 6);
    
    for j = 1:ngroups
        starStr = '';
        if pVal(j) < 0.001, starStr = '***';
        elseif pVal(j) < 0.01, starStr = '**';
        elseif pVal(j) < 0.05, starStr = '*'; 
        end
        
        if ~isempty(starStr)     
            offset = 0.05;             
            if meanVal(j) >= 0
                x_pos = meanVal(j) + seVal(j) + offset;
                ha = 'left'; 
            else
                x_pos = meanVal(j) - seVal(j) - offset;
                ha = 'right'; 
            end
            text(x_pos, y(j)-0.04, starStr, 'FontSize', 12, ...
                'HorizontalAlignment', ha, 'VerticalAlignment', 'middle'); 
        end
    end
end
xline(0, 'k-', 'LineWidth', 1);
xline(-0.5, ':', 'Color', [0.6 0.6 0.6]); 
xlabel('Deviation of urban-edge from wildland', 'FontSize', 12);
set(gca, 'Box', 'off', 'YColor', 'none', 'TickDir', 'out');
xticks(-0.9:0.3:1);
xlim([-1.3, 1]); 
for i = 1:nVars
    text(-1.3, i, xLabels{i}, 'HorizontalAlignment', 'left', ... 
        'VerticalAlignment', 'middle', 'FontSize', 12);
end
lgd = legend([b(2), b(1)], {'ShrubGrass', 'Forest'}, ...
    'Location', 'southeast', 'Box', 'off', 'FontSize', 11);
lgd.ItemTokenSize = [15, 18];
lgd.Position(2) = lgd.Position(2) + 0.08;
hold off;
exportgraphics(gcf, 'Fig/Fig3B.png');
%% NDVI calculation, Forest and ShrubGrass
lcs      = {'Forest', 'ShrubGrass'};              nlc  = numel(lcs);
decades  = {'1990s', '2000s', '2010s', '2020s'};  nDec  = numel(decades);
fireType = {'Urban-edge','Wildland'};             nFire = numel(fireType);
colors   = {[216, 118, 89]/255; [41, 157, 143]/255};
data     = shaperead('dataPrc/firePrmt/CalFire.shp');
ndviAll  = [data.ndviM]';
lc       = {data.lndcvr}';
Dec      = {data.decade}';
Fire     = {data.FireType}';

figure('Position', [100, 100, 500, 500]); 
t = tiledlayout(nDec-1, nlc, 'TileSpacing', 'compact', 'Padding', 'compact');
binW = 0.02;
edge = 0:binW:1;
xFit = 0:0.005:1; 
for d = 2:nDec
    for l = 1:nlc
        ax = nexttile;
        hold(ax, 'on');        
        if l == 1, maskLC = strcmp(lc, 'forest');
        else, maskLC = strcmp(lc, 'shrub') | strcmp(lc, 'grass'); end
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
        xlabel('NDVI');
        xlim(ax, [0.1 1]); 
        ylim(ax, [0 0.15]); 
        if d == 2, title(ax, lcs{l}, 'FontWeight','normal'); end        
        if l ~= 1, ax.YAxis.Visible = 'off';
            text(ax, 0.95, 0.70, decades{d}, 'Units', 'normalized', ...
                'HorizontalAlignment', 'right'); 
        end        
        if d ~= nDec, xticks(ax, []); end
        set(ax, 'FontSize', 10, 'TickDir', 'out', 'Box', 'off');

        if d == 3 && l == 1
            lgd = legend(ax, 'Location', 'northwest', 'Box', 'off');
            lgd.ItemTokenSize = [15, 15];
        end
    end
end
ylabel(t, 'Probability');
exportgraphics(gcf, 'Fig/Fig5B.png');
