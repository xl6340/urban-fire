%% Fig. 5A, trend
clear;clc;

FireType = {'Urban-edge','Wildland'};   nFire = numel(FireType);
FireTypeLabel = {'WUI fire','Wildland fire'}; 
igType   = {'Human', 'Natural'};        nIg      = numel(igType);
varList = {'Fire number', 'Burned area'};
var      = {'count', 'area'};           nVar = numel(var);  
colors = {'#fdd85d', '#99d6ea'};

yLable = {'Fire number (#)', 'Burned area (10^3 km^2)'};
figure('Position', [100, 100, 500, 350], 'Color', 'w');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 
for i = 1:nVar
    for j = 1:nFire
        ax = nexttile(t); hold on;  
        data = readtable(sprintf('dataFig/ba/%s-%s.csv', var{i}, FireType{j}));
        x = data.Row; 
        xfit = linspace(min(x), max(x), 100)';
        x_design = [ones(size(x)), x];
        for k = 1:nIg
            y = data.(igType{k});
            if i == 2, y = y./1000; end
            thisColor = colors{k};
            [b,~,~,~,stats] = regress(y, x_design);
            yfit = b(1) + b(2) * xfit;
            plot(x, y, '-o', 'Color', thisColor, ...
                    'MarkerFaceColor', thisColor, 'MarkerSize', 4, ...
                    'LineWidth', 0.01, 'DisplayName', igType{k});
            hold on;
            if stats(3) < 0.01
                plot(xfit, yfit, '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');           
            end   

            xlim([1990 2025]); xticks(1990:5:2025);xtickangle(30);  

            if j == 1, txtY = 0.4; ylabel(sprintf(yLable{i}));
            else, txtY = 0.6; yticklabels([]);end

            if i == 1
                ylim([0 400]);
                title(FireTypeLabel{j}, 'FontWeight','normal');xticklabels([]); 
                if stats(3) < 0.01
                text(0.2, txtY, sprintf('s = %.2f', b(2)), ...
                    'Interpreter','none','Units','normalized','Color','k','FontSize',10);
                end
            else 
                ylim([0 12]); yticks(0:4:12); 
                if j == 2, legend({'Human ignited', 'Natural ignited'}, 'Box','off', 'Location','northeast'); end 
                if stats(3) < 0.01
                    text(0.2, txtY, sprintf('$s = %.1f$', b(2)*1000),'Interpreter', 'latex', ...
                        'Units','normalized','Color','k','FontSize',10); 
                end
            end
             
        end
    end
end
fontsize(gcf, 9, "points");
exportgraphics(gcf,'Fig/Fig5A.pdf');
%% Fig.5B, beta, ignition, firetype
figure; clf; 

FireTypeLabel = {'Wildland','WUI'}; 
set(gcf, 'Position', [100, 100, 250, 150]); 
hold on; 
yPos = nFire:-1:1;
for j = 1:nFire
    data = readtable(sprintf('dataFig/beta/2Igns-%s.csv', FireType{j}));
    bVals = data.beta;     
    bErrs = data.betaErr;
    plot([bVals(1, 1), bVals(2, 1)], [yPos(j), yPos(j)], ...
        '-', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1, 'HandleVisibility', 'off');    
    for i = 1:nIg
        thisColor = colors{i};
        if j == 1, lbl = igType{i}; handVis = 'on';
        else, lbl = '';handVis = 'off';end
        errorbar(bVals(i, 1), yPos(j), bErrs(i, 1), 'horizontal', ...
            'LineStyle', 'none', ...      
            'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', thisColor, ...
            'Color', thisColor, 'LineWidth', 1.5, 'CapSize', 6, ...
            'DisplayName', lbl, 'HandleVisibility', handVis);
    end
end
set(gca, 'YTick', 1:nFire, 'YTickLabel', FireTypeLabel, 'FontSize', 9); 
ylim([0.5, nFire + 0.5]); 
xlabel('$\beta$ value', 'Interpreter', 'latex', 'FontSize', 11);
xlim([-1.8, -1]); 
xticks(-1.8:0.2:-1);
% legend('Location', 'best', 'Box', 'off');
hold off;
exportgraphics(gcf,'Fig/Fig5B.pdf');
%% Fig.5C,facility density
clear;clc;
facility = readtable('dataFig/facility.xlsx');
landscape = facility.landscape;
density = facility.density;
figure; clf; 
set(gcf, 'Position', [100, 100, 250, 150]); 
barh(landscape, density, 'FaceColor', '#f2e9e4');
set(gca, 'XTick', 0:4, 'FontSize', 9, 'YDir', 'reverse'); 
xlabel('Facility density (#/100-km^2)');
box off; hold off;
exportgraphics(gcf,'Fig/Fig5C.pdf');