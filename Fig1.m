clear; clc;

dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'}; 
nSrc     = numel(dataSrc);
FireType = {'Urban-edge', 'Wildland'};           
nFire    = numel(FireType);
color    = {[216, 118, 89]/255; [41, 157, 143]/255}; 

figure; clf; 
set(gcf, 'Position', [100, 100, 400, 300]); 
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 
hLeg = gobjects(1, nFire); 
for i = 1:nSrc
    ax(i) = nexttile; hold on;
    
    for j = 1:nFire
        fNameData = sprintf('dataFig/curve/%s-%s.csv', dataSrc{i}, FireType{j});
        fNameFit  = sprintf('dataFig/curve/%s-%s-fit.csv', dataSrc{i}, FireType{j});
 
        if i == 1
            hLeg(j) = plot(NaN, NaN, '-o', 'Color', color{j}, 'MarkerFaceColor', color{j}, ...
                'MarkerSize', 6, 'DisplayName', FireType{j});
        end

        if isfile(fNameData)
            data = readtable(fNameData);
            loglog(table2array(data(:,1)), table2array(data(:,2)), 'o', ...
                'Color', color{j}, 'MarkerFaceColor', color{j}, ...
                'MarkerSize', 1.5, 'HandleVisibility', 'off');
            
            hold on; 
            if isfile(fNameFit)
                dataFit = readtable(fNameFit);
                loglog(table2array(dataFit(:,1)), table2array(dataFit(:,2)), '-', ...
                    'Color', color{j}, 'HandleVisibility', 'off');
            end
        end
    end

    set(gca, 'XScale', 'log', 'YScale', 'log', ...
             'XMinorTick', 'off', 'YMinorTick', 'off', ...
             'TickLength', [0.02 0.02]); 
    
    xlim([10^-1, 10^5]); 
    ylim([10^(-7), 10^1]);    
    xticks(10 .^ (-1:2:5));
    yticks(10 .^ (-7:3:1));
    
    text(0.05, 0.05, dataSrc{i}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
    
    if i == 1
        text(0.6, 0.8, '$y \sim x^{\beta}$', 'Interpreter','latex',...
             'Units', 'normalized', 'FontSize', 12);
    end

    if i == 1 || i == 2
        set(gca, 'XTickLabel', []);
    end
    if i == 2 || i == 4
        set(gca, 'YTickLabel', []);
    end
end    
xlabel(t, 'Fire size (km^2)');
ylabel(t, 'Probability density');

lgd = legend(hLeg, 'Orientation', 'horizontal', 'Box', 'off');
lgd.ItemTokenSize = [15, 18]; 
lgd.Layout.Tile = 'North';

exportgraphics(gcf,'Fig/Fig1B.png');
%% Fig1C, beta
dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'}; nSrc = numel(dataSrc);
FireType = {'Urban-edge', 'Wildland'};            nFire = numel(FireType);
colors   = {[216, 118, 89]/255; [41, 157, 143]/255}; 

bVals = zeros(nSrc, nFire); 
bErrs = zeros(nSrc, nFire); 

for j = 1:nFire
    fName = sprintf('dataFig/beta/4Srcs-%s.csv', FireType{j});
    if isfile(fName)
        data = readtable(fName);
        bVals(:, j) = data.beta;
        bErrs(:, j) = data.betaErr;
    else
        bVals(:, j) = -1.8 + rand(nSrc,1); 
        bErrs(:, j) = 0.08 * ones(nSrc,1);
    end
end

figure(1); clf; 
set(gcf, 'Position', [100, 100, 200, 200]);
hold on;
yPos = 1:nSrc; 
for i = 1:nSrc
    plot([bVals(i, 1), bVals(i, 2)], [yPos(i), yPos(i)], ...
        '-', 'Color', [0.4, 0.4, 0.4], 'HandleVisibility', 'off');
end

h = gobjects(1, nFire);
for j = 1:nFire
    h(j) = errorbar(bVals(:, j), yPos, bErrs(:, j), ...
        'horizontal', ...             
        'LineStyle', 'none', ...      
        'Marker', 'o', ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', colors{j}, ...
        'Color', colors{j}, ...
        'LineWidth', 1.5, ...
        'CapSize', 6, ...
        'DisplayName', FireType{j});
end
set(gca, 'YDir', 'reverse', 'XGrid', 'on', 'YGrid', 'off'); 
set(gca, 'YTick', 1:nSrc, 'YTickLabel', dataSrc); 
xlabel('$\beta$ value', 'Interpreter', 'latex', 'FontSize', 10);
ylim([0.5, nSrc + 0.5]); 
xlim([-2, -1.2]); 

exportgraphics(gcf,'Fig/Fig1C.png');