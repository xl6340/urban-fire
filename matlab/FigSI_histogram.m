clear; clc;

varNames  = {'vpdmax', 'tmean', 'tmin', 'tmax', 'ppt', 'vs', 'ndviM', 'slope', ...
             'elevation', 'FFMC', 'DMC', 'DC'};
yLabels   = {'VPD_{max} (kPa)', 'T_{mean} (°C)', 'T_{min} (°C)', 'T_{max} (°C)', ...
             'Precipitation (mm)', 'Wind speed (m/s)', 'NDVI', 'Slope (°)', ...
             'Elevation (m)', 'FFMC', 'DMC', 'DC'};
subLabels = {'(A)', '(B)', '(C)', '(D)', '(E)', '(F)', '(G)', '(H)', '(I)', '(J)', '(K)', '(L)'};
xLimits   = [ 
    0,   8;     % VPDmax (kPa)
    0,   40;    % Tmean
   -5,   30;    % Tmin
    0,   50;    % Tmax
    0,   3000;  % Precip
    0,   10;    % VS
    0.2, 1;     % NDVI
    0,   35;    % Slope
   -500, 3000;  % Elevation
   70,   101;   % FFMC
    0,   900;   % DMC
    0,   1800]; % DC

nVars     = numel(varNames);
lc        = {'Forest', 'ShrubGrass'};
fireType  = {'Wildland', 'Urban-edge'};
decade    = {'1990s', '2000s', '2010s', '2020s'};
cWild     = [66, 157, 143]/255;   
cUrban    = [231, 111, 81]/255;   

Data.Forest.Urban = table(); Data.Forest.Wild  = table();
Data.Shrub.Urban  = table(); Data.Shrub.Wild   = table();
fprintf('Loading data...\n');
for i = 1:numel(lc)
    for j = 1:numel(fireType)
        fwiName = sprintf('dataFig/variable/%s_%s_withFWI.csv', lc{i}, fireType{j});
        rawName = sprintf('dataFig/variable/%s_%s.csv', lc{i}, fireType{j});
        if isfile(fwiName)
            tmp = readtable(fwiName);
        elseif isfile(rawName)
            tmp = readtable(rawName);
        else
            tmp = table();
            for k = 1:numel(decade)
                fName = sprintf('dataFig/variable/%s_%s_%s.csv', lc{i}, fireType{j}, decade{k});
                if isfile(fName)
                    tmp = [tmp; readtable(fName)];
                end
            end
        end
        if strcmp(lc{i}, 'Forest')
            if strcmp(fireType{j}, 'Urban-edge')
                Data.Forest.Urban = tmp;
            else
                Data.Forest.Wild = tmp;
            end
        elseif strcmp(lc{i}, 'ShrubGrass')
            if strcmp(fireType{j}, 'Urban-edge')
                Data.Shrub.Urban = tmp;
            else
                Data.Shrub.Wild = tmp;
            end
        end
    end
end

%% Plotting Variable Histograms: Forest Only
figure('Color', 'w', 'Position', [100, 100, 1400, 900]);
t = tiledlayout(3, 4, 'TileSpacing', 'loose', 'Padding', 'compact');
for i = 1:nVars
    vName = varNames{i};
    ax = nexttile;
    hold(ax, 'on');
    
    val_U = Data.Forest.Urban.(vName); val_U = val_U(~isnan(val_U));
    val_W = Data.Forest.Wild.(vName);  val_W = val_W(~isnan(val_W));
    
    if strcmp(vName, 'vpdmax')
        val_U = val_U / 10;
        val_W = val_W / 10;
    end
    
    if isempty(val_U) || isempty(val_W), continue; end

    mu_U = mean(val_U); sig_U = std(val_U);
    mu_W = mean(val_W); sig_W = std(val_W);
    
    diff_val = mu_U - mu_W;    
    [~, p_val] = ttest2(val_U, val_W, 'Vartype', 'unequal');    
    xRange = linspace(xLimits(i,1), xLimits(i,2), 100);    

    % Wildland (Teal)
    hW = histogram(ax, val_W, 'Normalization', 'pdf', ...
        'FaceColor', cWild, 'EdgeColor', [0.9 0.9 0.9], 'FaceAlpha', 0.4);
    plot(ax, xRange, normpdf(xRange, mu_W, sig_W), 'Color', cWild, 'LineWidth', 2);
    xline(ax, mu_W, ':', 'Color', cWild, 'LineWidth', 2);
    
    % Urban-edge (Orange)
    hU = histogram(ax, val_U, 'Normalization', 'pdf', ...
        'FaceColor', cUrban, 'EdgeColor', [0.9 0.9 0.9], 'FaceAlpha', 0.4);
    plot(ax, xRange, normpdf(xRange, mu_U, sig_U), 'Color', cUrban, 'LineWidth', 2);
    xline(ax, mu_U, ':', 'Color', cUrban, 'LineWidth', 2);
    
    % Formatting
    xlabel(ax, yLabels{i});
    xlim(ax, xLimits(i, :));
    
    % Generate Title String with Difference
    if p_val < 0.001
        pStr = '{\itp}<0.001';
    elseif p_val < 0.01
        pStr = sprintf('{\\itp}=%.3f', p_val);
    else
        pStr = sprintf('{\\itp}=%.2f', p_val);
    end
    
    % Format: (A) Diff:+1.5 (p<0.001)
    titleStr = sprintf('%s  Diff:%+.1f  (%s)', subLabels{i}, diff_val, pStr);
    
    title(ax, titleStr, 'FontWeight', 'normal', ...
        'Units', 'normalized', 'Position', [0.02, 1.02, 0], 'HorizontalAlignment', 'left');
    
    set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10, 'YColor', 'none');
    
    if i == 1
        lgd = legend([hW, hU], {'Wildland', 'Urban-edge'}, ...
            'Location', 'northeast', 'Orientation', 'vertical');
        lgd.Box = 'off';
    end
end
ylabel(t, 'Probability Density', 'FontSize', 12);
exportgraphics(gcf, 'Fig/FigSI_histogram_forest.png');
exportgraphics(gcf, 'Fig/FigS5.png');
%% Plotting Variable Histograms: ShrubGrass Only
xLimits   = [ 
    0,   8;     % VPDmax (kPa)
    0,   40;    % Tmean
   -5,   30;    % Tmin
    0,   50;    % Tmax
    0,   1500;  % Precip
    0,   11;    % VS
    0,   1;     % NDVI
    0,   35;    % Slope
    -500, 3000];  % Elevation

figure('Color', 'w', 'Position', [100, 100, 1200, 900]);
t = tiledlayout(3, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
for i = 7 %1:nVars
    vName = varNames{i};
    ax = nexttile;
    hold(ax, 'on');
    
    % --- CHANGE: Switch from Forest to Shrub data ---
    val_U = Data.Shrub.Urban.(vName); val_U = val_U(~isnan(val_U));
    val_W = Data.Shrub.Wild.(vName);  val_W = val_W(~isnan(val_W));
    
    % Unit Conversion (VPD hPa -> kPa)
    if strcmp(vName, 'vpdmax')
        val_U = val_U / 10;
        val_W = val_W / 10;
    end
    
    if isempty(val_U) || isempty(val_W), continue; end
    
    % Statistics
    mu_U = mean(val_U); sig_U = std(val_U);
    mu_W = mean(val_W); sig_W = std(val_W);
    
    diff_val = mu_U - mu_W;    
    [~, p_val] = ttest2(val_U, val_W, 'Vartype', 'unequal');    
    
    % Plotting Ranges
    xRange = linspace(xLimits(i,1), xLimits(i,2), 100);    
    
    % Wildland (Teal)
    hW = histogram(ax, val_W, 'Normalization', 'pdf', ...
        'FaceColor', cWild, 'EdgeColor', [0.9 0.9 0.9], 'FaceAlpha', 0.4);
    plot(ax, xRange, normpdf(xRange, mu_W, sig_W), 'Color', cWild, 'LineWidth', 2);
    xline(ax, mu_W, ':', 'Color', cWild, 'LineWidth', 2);
    
    % Urban-edge (Orange)
    hU = histogram(ax, val_U, 'Normalization', 'pdf', ...
        'FaceColor', cUrban, 'EdgeColor', [0.9 0.9 0.9], 'FaceAlpha', 0.4);
    plot(ax, xRange, normpdf(xRange, mu_U, sig_U), 'Color', cUrban, 'LineWidth', 2);
    xline(ax, mu_U, ':', 'Color', cUrban, 'LineWidth', 2);
    
    % Formatting
    xlabel(ax, yLabels{i});
    xlim(ax, xLimits(i, :));
    
    % Title with Difference and P-value
    if p_val < 0.001
        pStr = '{\itp}<0.001';
    elseif p_val < 0.01
        pStr = sprintf('{\\itp}=%.3f', p_val);
    else
        pStr = sprintf('{\\itp}=%.2f', p_val);
    end
    
    titleStr = sprintf('%s  Diff:%+.4f  (%s)', subLabels{i}, diff_val, pStr);
    
    title(ax, titleStr, 'FontWeight', 'normal', ...
        'Units', 'normalized', 'Position', [0.02, 1.02, 0], 'HorizontalAlignment', 'left');
    
    set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10, 'YColor', 'none');
    
    % Legend on first tile
    if i == 1
        lgd = legend([hW, hU], {'Wildland', 'Urban-edge'}, ...
            'Location', 'northeast', 'Orientation', 'vertical');
        lgd.Box = 'off';
    end
end

ylabel(t, 'Probability Density', 'FontSize', 12);
exportgraphics(gcf, 'Fig/FigSI_histogram_shrub.png');