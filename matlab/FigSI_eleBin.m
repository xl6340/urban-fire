clear; clc;
% forest
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
fireType  = {'Urban-edge','Wildland'};              nFire = numel(fireType);

nBins = 8; 
eleEdges = unique(quantile(elevations(strcmp(lcs, lc{1})), linspace(0, 1, nBins+1)));

count    = NaN(nFire, nBins);
area     = NaN(nFire, nBins);
beta     = NaN(nFire, nBins);
betaErr  = NaN(nFire, nBins);
pBeta    = NaN(nFire, nBins);
R2       = NaN(nFire, nBins);
sizeEdges = 10 .^ (0:0.05:5); 
for i = 1:nFire
    for j = 1:nBins
        minEle = eleEdges(j);
        maxEle = eleEdges(j+1);
        
        idx = strcmp(Fires, fireType{i}) & ...
              strcmp(lcs, lc{1}) & ...
              elevations >= minEle & ...
              elevations < maxEle;          
        
        sz = FireSizes(idx);
        if numel(sz) < 10 
            continue; 
        end
        
        [counts, binEdges] = histcounts(sz, sizeEdges);
        binWidth = diff(binEdges);
        
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
    
        validIdx = probDensity > 0; 
        x = log10(binCenters(validIdx))'; 
        y = log10(probDensity(validIdx))';
        
        if length(x) < 10, continue; end        
        mdl = fitlm(x, y);       
        
        count(i,j)   = numel(sz);
        area(i,j)    = sum(sz);        
        beta(i,j)    = mdl.Coefficients.Estimate(2);
        betaErr(i,j) = mdl.Coefficients.SE(2);
        pBeta(i,j)   = mdl.Coefficients.pValue(2);
        R2(i,j)      = mdl.Rsquared.Ordinary;            
    end
end
%% plotting
figure; clf; 
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 350]); % Wider window
t = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

eleCenters = (eleEdges(1:end-1) + eleEdges(2:end)) / 2;

% --- Subplot 1: Elevation Histogram ---
ax1 = nexttile; 
hold(ax1, 'on');
histogram(ax1, elevations(strcmp(lcs, lc{1})), 40, 'Normalization', 'pdf', ...
    'FaceColor', '#7785ac', 'EdgeColor', [0.9 0.9 0.9], 'FaceAlpha', 0.4);

% CHANGE: Plot vertical lines at bin CENTERS to match Figure B
for k = 1:length(eleCenters)
    xline(ax1, eleCenters(k),'LineStyle', '--','Color', '#7209b7','LineWidth', 1.5);
end
xlim([0 3000]);
xlabel(ax1, 'Elevation (m)', 'FontSize', 12);
ylabel(ax1, 'Probability density', 'FontSize', 12);
set(ax1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);

% --- Subplot 2: Scatter of Beta Values ---
ax2 = nexttile; 
hold(ax2, 'on');
color = {[216, 118, 89]/255; [41, 157, 143]/255}; % Orange / Teal

hUrban = errorbar(ax2, eleCenters, beta(1,:), betaErr(1,:), '-o', ...
    'Color', color{1}, 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', 'CapSize', 10);
hWild = errorbar(ax2, eleCenters, beta(2,:), betaErr(2,:), '-^', ...
    'Color', color{2}, 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', 'CapSize', 10);

xlabel(ax2, 'Elevation (m)', 'FontSize', 12);
ylabel(ax2, '\beta value', 'FontSize', 12);
ylim(ax2, [-1.6 -0.8]); 
yticks(ax2, -1.6:0.2:-0.8);
set(ax2, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);

lgd = legend([hUrban, hWild], {'Urban-edge', 'Wildland'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
exportgraphics(gcf, 'Fig/FigSI_eleBin.png');