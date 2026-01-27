%% Fig1A, fire type
clear; clc;

fireType = {'Urban-edge', 'Wildland'};     
decades  = {'1990s', '2000s', '2010s', '2020s'}; 
edgeCut  = {'<100 km^2', '<1000 km^2', 'Full size'}; 

nFire = numel(fireType);
nDec  = numel(decades);
nCut  = numel(edgeCut);

Y_data = nan(nCut, nDec, nFire);
E_data = nan(nCut, nDec, nFire);
for f = 1:nFire
    for d = 1:nDec
        data = readtable(sprintf('dataFig/beta/3cutSizes-%s-%s.csv', fireType{f}, decades{d})); 
        Y_data(:, d, f) = data.beta;
        E_data(:, d, f) = data.betaErr;
    end
end

colors = {[236, 190, 178;  205, 224, 199] / 255;  % cut 1 (1-100)
          [227, 164, 151;  155, 184, 156] / 255;  % cut 2 (1-1000)
          [216 118 89;  41 157 143] / 255};       % cut 3 (full size)

x_offsets = [0.2, 0.1, 0];
markers   = {'o', '^'};

figure('Position', [100, 100, 400, 350]); 
ax = gca; hold on; 
hStore = gobjects(nCut, nFire);
for f = 1:nFire
    for c = 1:nCut
        x = (1:nDec) + x_offsets(c);
        y = Y_data(c, :, f);
        e = E_data(c, :, f);
        color = colors{c}(f,:);
        errorbar(x, y, e, 'Color', color, 'LineStyle', 'none', ...
            'LineWidth', 1.5, 'CapSize', 4, 'HandleVisibility', 'off');

        h = plot(x, y, 'Marker', markers{f}, 'Color',  color, 'MarkerFaceColor',  color, ...
            'MarkerEdgeColor',  color, 'LineWidth', 1, 'MarkerSize', 5);
        hStore(c, f) = h;
    end
end

box off;
xlim([0.5, 4.5]); 
xticks(1:4); 
xticklabels(decades); 
xtickangle(20);
ylim([-2.0, -1.0]); 
ylabel('\beta value', 'Interpreter', 'tex', 'FontSize', 12);
yticks(-2:0.2:-1);


hUrban = flipud(hStore(:, 1));
hWild  = flipud(hStore(:, 2));
legendLabels = flipud(edgeCut'); 
lgd2 = legend(ax, hWild, legendLabels, 'Location', 'southeast','Box', 'off', 'FontSize', 12);
lgd2.ItemTokenSize = [15, 18];
ax2 = axes('Position', ax.Position, 'Visible', 'off');
lgd1 = legend(ax2, hUrban, {'','',''}, 'Location', 'southeast', 'Box', 'off', 'Color', 'none'); 
lgd1.ItemTokenSize = [15, 18]; 
drawnow;
pos2 = lgd2.Position; 
pos1 = lgd1.Position;
shift_amount = 0.10; 
lgd1.Position = [pos2(1) - shift_amount, pos2(2), pos1(3), pos2(4)];

exportgraphics(gcf, 'Fig/Fig2A.pdf');
%% Figure 2b
clear; clc;
data = readtable('dataFig/vpd/betaVPD.xlsx');
vpd_fire   = data.vpd_fire_kPa;
vpd_region = data.vpd_region_kPa;
beta_CA    = data.beta_CA;
beta_err   = data.beta_err;
periods    = data.Period; 

vpd_wui = data.vpd_wui_kPa;
beta_urban = data.beta_urban;
betaErr_urban = data.beta_err_1;

vpd_wild = data.vpd_wild_kPa;
beta_wild = data.beta_wild;
betaErr_wild = data.beta_err_2;
colorRegion = 'k';                     % Black
colorUrban  = [216, 118, 89]/255;      % Red-Orange
colorWild   = [41, 157, 143]/255;      % Green

figure('Position', [100, 100, 400, 350]); hold on;
% --- 1. Plot Urban (Red-Orange) ---
h1 = errorbar(vpd_wui, beta_urban, betaErr_urban, 'o', ...
    'Color', colorUrban, 'MarkerFaceColor', colorUrban, ...
    'MarkerSize', 6, 'LineWidth', 1, 'CapSize', 0, 'LineStyle', 'none', ...
    'DisplayName', 'Urban-Edge Fires');

idx3 = ~isnan(vpd_wui) & ~isnan(beta_urban);
if any(idx3)
    p = polyfit(vpd_wui(idx3), beta_urban(idx3), 1);
    x_fit = linspace(min(vpd_wui(idx3)), max(vpd_wui(idx3)), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, '-', 'Color', colorUrban, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Calc stats
    [R, P] = corrcoef(vpd_wui(idx3), beta_urban(idx3));
    % Display stats (Middle Left)
    text(0.05, 0.90, sprintf('Urban-edge: r=%.2f, p=%.3f', R(1,2), P(1,2)), ...
        'Units', 'normalized', 'Color', colorUrban, 'FontSize', 12);
    for i = 1:length(vpd_wui)
        if idx3(i)
            txt = string(periods{i}); 
            text(vpd_wui(i)-0.02, beta_urban(i) - 0.03, txt, ...
                'Color', [0.6 0.6 0.6], 'FontSize', 8, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
end

% --- 2. Plot Wildland (Green) ---
h2 = errorbar(vpd_wild, beta_wild, betaErr_wild, '^', ...
    'Color', colorWild, 'MarkerFaceColor', colorWild, ...
    'MarkerSize', 6, 'LineWidth', 1, 'CapSize', 0, 'LineStyle', 'none', ...
    'DisplayName', 'Wildland Fires');

idx4 = ~isnan(vpd_wild) & ~isnan(beta_wild);
if any(idx4)
    p = polyfit(vpd_wild(idx4), beta_wild(idx4), 1);
    x_fit = linspace(min(vpd_wild(idx4)), max(vpd_wild(idx4)), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, '-', 'Color', colorWild, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Calc stats
    [R, P] = corrcoef(vpd_wild(idx4), beta_wild(idx4));
    % Display stats (Bottom Left)
    text(0.05, 0.80, sprintf('Wildland: r=%.2f, p=%.3f', R(1,2), P(1,2)), ...
        'Units', 'normalized', 'Color', colorWild, 'FontSize', 12);
    
    % --- ADD LABELS FOR WILDLAND (Placed Below) ---
    for i = 1:length(vpd_wild)
        if idx4(i)
            txt = string(periods{i});
            text(vpd_wild(i)+0.02, beta_wild(i) - 0.01, txt, ...
                'Color', [0.6 0.6 0.6], 'FontSize', 8, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        end
    end
end

xlabel('VPD_{max} (kPa)', 'FontSize', 12, 'Interpreter', 'tex');
ylabel('\beta value', 'FontSize', 12, 'Interpreter', 'tex');
xlim([1.7 2.1]); xticks(1.7:0.1:2.1); 
ylim([-1.7 -1]); yticks(-1.7:0.2:-1);
box off;
hold off;
exportgraphics(gcf,'Fig/Fig2B.pdf');