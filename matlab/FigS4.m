clear;clc;

elevations = readmatrix('dataFig/elevation/elevations.csv');
eleCenters = readmatrix('dataFig/elevation/eleCenters.csv');

figure; clf; set(gcf, 'Position', [100, 50, 600, 600]); 
t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = {[216, 118, 89]/255; [41, 157, 143]/255}; % Orange, Teal

% elevation histogram
ax1 = nexttile; hold(ax1, 'on');
histogram(ax1, elevations, 40, 'Normalization', 'pdf', ...
    'FaceColor', '#7785ac', 'EdgeColor', [0.9 0.9 0.9], 'FaceAlpha', 0.4);
for k = 1:length(eleCenters)
    xline(ax1, eleCenters(k),'LineStyle', '-.','Color', '#7209b7','LineWidth', 1.5);
end
xlim([0 3000]);
ylabel(ax1, 'Probability density', 'FontSize', 11);
set(ax1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);

% beta and climate variables
Vars = {'beta', 'VPDmax', 'Tmean', 'Tmin', 'Tmax'}; nVar = numel(Vars);
yLabels   = {'\beta value', 'VPD_{max} (kPa)', 'T_{mean} (°C)', 'T_{min} (°C)', 'T_{max} (°C)'};
ylims     = {[-1.6 -0.8], [2 5], [15 25], [8 16], [20 35]};
yTickVals = {-1.6:0.2:-0.8, 2:1:5, 15:5:25, 8:4:16, 20:5:35};
for k = 1:nVar
    ax = nexttile; hold(ax, 'on');
    data = readtable(sprintf('dataFig/elevation/%s.csv', Vars{k}));
    % Urban
    errorbar(ax, eleCenters, data.Mean_Urban_edge, data.Err_Urban_edge, '-o', ...
        'Color', colors{1},'MarkerSize', 5, 'MarkerFaceColor', colors{1});
    
    % Wildland
    errorbar(ax, eleCenters, data.Mean_Wildland, data.Err_Wildland, '-^', ...
        'Color', colors{2},'MarkerSize', 5, 'MarkerFaceColor', colors{2});
    
    if k >3, xlabel(ax, 'Elevation (m)', 'FontSize', 11);end
    ylabel(ax, yLabels{k}, 'FontSize', 11);
    set(ax, 'Box', 'off', 'TickDir', 'out');
    xlim(ax, [0 3000]);
    ylim(ax, ylims{k});
    yticks(ax, yTickVals{k});
end
exportgraphics(gcf, 'Fig/Fig4.pdf');