clear; clc;

data = shaperead('dataPrc/firePrmt/CalFire.shp');

lcs        = {data.lc}';
Fires      = {data.FireType}';
elevations = [data.elevation]';
FireSizes  = [data.size]'; 

Tmeans     = [data.tmean]'; 
Tmaxs      = [data.tmax]';
Tmins      = [data.tmin]';
VPDmaxs    = [data.vpdmax]'; 

lc        = {'Forest', 'Shrub', 'Grass'};                                         
fireType  = {'Urban-edge','Wildland'};              
nFire     = numel(fireType);
nBins     = 8; 

eleEdges = unique(quantile(elevations(strcmp(lcs, lc{1})), linspace(0, 1, nBins+1)));
eleCenters = (eleEdges(1:end-1) + eleEdges(2:end)) / 2;

vars = {'Tmean', 'Tmax', 'Tmin', 'VPDmax'};
nVars = length(vars);

results = struct();
for v = 1:nVars
    results.(vars{v}).Val = NaN(nFire, nBins);
    results.(vars{v}).Err = NaN(nFire, nBins);
end

allVarData = {Tmeans, Tmaxs, Tmins, VPDmaxs};
for i = 1:nFire
    for j = 1:nBins
        minEle = eleEdges(j);
        maxEle = eleEdges(j+1);
        
        idx = strcmp(Fires, fireType{i}) & ...
              strcmp(lcs, lc{1}) & ...
              elevations >= minEle & ...
              elevations < maxEle;    
        if sum(idx) < 10, continue; end
        for v = 1:nVars
            currentData = allVarData{v};
            binData = currentData(idx);
            
            m = mean(binData, 'omitnan');
            e = std(binData, 'omitnan') / sqrt(numel(binData));
            
            results.(vars{v}).Val(i,j) = m;
            results.(vars{v}).Err(i,j) = e;
        end
    end
end

figure; clf; 
set(gcf, 'Color', 'w', 'Position', [50, 50, 1000, 600]); 
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Elevation Histogram 
ax1 = nexttile; 
hold(ax1, 'on');
histogram(ax1, elevations(strcmp(lcs, lc{1})), 40, 'Normalization', 'pdf', ...
    'FaceColor', '#7785ac', 'EdgeColor', [0.9 0.9 0.9], 'FaceAlpha', 0.4);
for k = 1:length(eleCenters)
    xline(ax1, eleCenters(k),'LineStyle', '--','Color', '#7209b7','LineWidth', 1.5);
end
xlim([0 3000]);
xlabel(ax1, 'Elevation (m)', 'FontSize', 11);
ylabel(ax1, 'PDF', 'FontSize', 11);
title(ax1, 'Elevation Dist.', 'FontWeight', 'bold');
set(ax1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);

% Subplot 2-5: other climate variables 
yLabs    = {'T_{mean} (°C)', 'T_{max} (°C)', 'T_{min} (°C)', 'VPD_{max} (kPa)'};
color    = {[216, 118, 89]/255; [41, 157, 143]/255}; % Orange / Teal
for v = 1:nVars
    ax = nexttile;
    hold(ax, 'on');
    varName = vars{v};    
    val = results.(varName).Val;
    err = results.(varName).Err;
    
    hUrban = errorbar(ax, eleCenters, val(1,:), err(1,:), '-o', ...
        'Color', color{1}, 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'MarkerFaceColor', 'w', 'CapSize', 8);    
    hWild = errorbar(ax, eleCenters, val(2,:), err(2,:), '-^', ...
        'Color', color{2}, 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'MarkerFaceColor', 'w', 'CapSize', 8);    
    xlabel(ax, 'Elevation (m)', 'FontSize', 11);
    ylabel(ax, yLabs{v}, 'FontSize', 11);
    set(ax, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);    
    if v == 1
        legend([hUrban, hWild], {'Urban-edge', 'Wildland'}, ...
            'Location', 'best', 'Box', 'off', 'FontSize', 9);
    end
end
exportgraphics(gcf, 'Fig/FigSI_eleBin_EnvVars.png');