%% NDVI calculation, ShrubGrass only
clear;clc;

lcs      = {'Shrub & grassland'};
nlc      = numel(lcs);
decades  = {'1990s', '2000s', '2010s', '2020s'};  nDec  = numel(decades);
fireType = {'Urban-edge','Wildland'};             nFire = numel(fireType);

data     = shaperead('dataPrc/firePrmt/CalFire.shp');
ndviAll  = [data.ndviM]';
lc       = {data.lndcvr}';
Dec      = {data.decade}';
Fire     = {data.FireType}';

colors   = {[216, 118, 89]/255; [41, 157, 143]/255};

figure('Position', [100, 100, 400, 500]); 
t = tiledlayout(nDec-1, nlc, 'TileSpacing', 'compact', 'Padding', 'compact');
binW = 0.02;
edge = 0:binW:1;
xFit = 0:0.005:1; 
for d = 2:nDec
    for l = 1:nlc
        ax = nexttile;
        hold(ax, 'on');      
        maskLC = strcmp(lc, 'shrub') | strcmp(lc, 'grass'); 
        
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
        
        xlim(ax, [0.1 0.9]); 
        ylim(ax, [0 0.15]); 
        
        if d == 2, title(ax, lcs{l}, 'FontWeight','normal'); end        
        
        text(ax, 0.95, 0.80, decades{d}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right'); 
       
        if d ~= nDec, xticks(ax, []); else, xlabel('NDVI'); end
        set(ax, 'FontSize', 10, 'TickDir', 'out', 'Box', 'off');
        
        if d == 2
            lgd = legend(ax, 'Location', 'east', 'Box', 'off');
            lgd.ItemTokenSize = [15, 15];
        end
    end
end
ylabel(t, 'Probability');

exportgraphics(gcf, 'Fig/Fig6.pdf');