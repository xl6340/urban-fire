clear; clc;

% --- Data Loading ---
data = shaperead('dataPrc/firePrmt/CalFire.shp');

% Extract Variables
lcs        = {data.lc}';
Fires      = {data.FireType}';
elevations = [data.elevation]';
FireSizes  = [data.size]'; 
Tmeans     = [data.tmean]'; 
Tmaxs      = [data.tmax]';
Tmins      = [data.tmin]';

% Convert VPD from hPa to kPa (divide by 10)
VPDmaxs    = [data.vpdmax]' / 10; 

% --- Configurations ---
lc        = {'Forest', 'Shrub', 'Grass'};                                         
fireType  = {'Urban-edge','Wildland'};              
nFire     = numel(fireType);
nBins     = 8; 

% Define Elevation Bins
eleEdges = unique(quantile(elevations(strcmp(lcs, lc{1})), linspace(0, 1, nBins+1)));

% Initialize Storage
beta     = NaN(nFire, nBins);
betaErr  = NaN(nFire, nBins);
% Storage for Env Variables: 1=Tmean, 2=Tmax, 3=Tmin, 4=VPDmax
envMeans = NaN(nFire, nBins, 4); 
envErrs  = NaN(nFire, nBins, 4);

sizeEdges = 10 .^ (0:0.05:5); 

% --- Processing Loop ---
for i = 1:nFire
    for j = 1:nBins
        minEle = eleEdges(j);
        maxEle = eleEdges(j+1);
        
        idx = strcmp(Fires, fireType{i}) & ...
              strcmp(lcs, lc{1}) & ...
              elevations >= minEle & ...
              elevations < maxEle;          
        
        sz = FireSizes(idx);
        if numel(sz) < 10, continue; end
        
        % 1. Calculate Beta
        [counts, binEdges] = histcounts(sz, sizeEdges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCentersIdx = sqrt(binEdges(1:end-1) .* binEdges(2:end));
        validIdx = probDensity > 0; 
        x = log10(binCentersIdx(validIdx))'; 
        y = log10(probDensity(validIdx))';
        
        if length(x) >= 10
            mdl = fitlm(x, y);       
            beta(i,j)    = mdl.Coefficients.Estimate(2);
            betaErr(i,j) = mdl.Coefficients.SE(2);
        end
        
        % 2. Calculate Environmental Stats
        vars = {Tmeans, Tmaxs, Tmins, VPDmaxs};
        for v = 1:4
            vals = vars{v}(idx);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                envMeans(i,j,v) = mean(vals);
                envErrs(i,j,v)  = std(vals) / sqrt(numel(vals));
            end
        end
    end
end

%% --- Plotting ---
figure; clf; 
set(gcf, 'Color', 'w', 'Position', [100, 50, 600, 600]); 

% 3 Rows, 2 Columns
t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

eleCenters = (eleEdges(1:end-1) + eleEdges(2:end)) / 2;
colors = {[216, 118, 89]/255; [41, 157, 143]/255}; % Orange, Teal

% --- 1. Beta Plot (Row 1, Col 1) ---
ax1 = nexttile; 
hold(ax1, 'on');
hUrban = errorbar(ax1, eleCenters, beta(1,:), betaErr(1,:), '-o', ...
    'Color', colors{1}, 'MarkerSize', 5, 'MarkerFaceColor', colors{1});
plot(ax1, eleCenters, beta(2,:), '-', 'Color', colors{2}, 'LineWidth', 1);
hWild = errorbar(ax1, eleCenters, beta(2,:), betaErr(2,:), '-^', ...
    'Color', colors{2}, 'MarkerSize', 5, 'MarkerFaceColor',colors{2});

% xlabel(ax1, 'Elevation (m)', 'FontSize', 11);
ylabel(ax1, '\beta value', 'FontSize', 11);
set(ax1, 'Box', 'off', 'TickDir', 'out');
xlim(ax1, [0 3000]);
ylim(ax1, [-1.6 -0.8]);

% --- 2. Legend (Row 1, Col 2) ---
lgd = legend([hUrban, hWild], {'Urban-edge', 'Wildland'}, ...
    'Box', 'off', 'FontSize', 12);
lgd.Layout.Tile = 2;

% --- 3. Env Plots (Rows 2 & 3) ---
% Order: VPD, Tmean, Tmin, Tmax
plotOrder = [4, 1, 3, 2]; 
yLabels   = {'VPD_{max} (kPa)', 'T_{mean} (°C)', 'T_{min} (°C)', 'T_{max} (°C)'};
ylims     = {[2 5], [15 25], [8 16], [20 35]};
yTickVals = {2:1:5, 15:5:25, 8:4:16, 20:5:35};

for k = 1:4
    idx = plotOrder(k);
    ax = nexttile; % This will automatically fill tiles 3, 4, 5, 6
    hold(ax, 'on');
    
    % Urban
    errorbar(ax, eleCenters, envMeans(1,:,idx), envErrs(1,:,idx), '-o', ...
        'Color', colors{1},'MarkerSize', 5, 'MarkerFaceColor', colors{1});
    
    % Wildland
    errorbar(ax, eleCenters, envMeans(2,:,idx), envErrs(2,:,idx), '-^', ...
        'Color', colors{2},'MarkerSize', 5, 'MarkerFaceColor', colors{2});
    
    if k >2, xlabel(ax, 'Elevation (m)', 'FontSize', 11);end
    ylabel(ax, yLabels{k}, 'FontSize', 11);
    set(ax, 'Box', 'off', 'TickDir', 'out');
    xlim(ax, [0 3000]);
    ylim(ax, ylims{k});
    yticks(ax, yTickVals{k});
end

exportgraphics(gcf, 'Fig/Fig6_eleBins.png');