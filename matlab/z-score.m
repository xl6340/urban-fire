%% z-score
clear; clc;
varNames  = {'ndviM', 'elevation', 'slope', 'vs', 'vd', ... 
            'ppt', 'tmean', 'vpdmax','tmin', 'tmax'};  

xLabels   = {'NDVI', 'Elevation', 'Slope', 'Wind Spd', 'Wind Dir', ...
             'P', 'T_{mean}', 'VPD_{max}', 'T_{min}', 'T_{max}'};
nVars     = numel(varNames);

cForest = [50, 100, 50]/255;    % Dark Green (Forest)
cShrub  = [100, 180, 100]/255;  % Light Green (Shrub)

lc       = {'Forest', 'ShrubGrass'};
fireType = {'Urban-edge','Wildland'};
decade   = {'1990s', '2000s', '2010s', '2020s'};

% --- 2. DATA LOADING ---
Data.Forest.Urban = table(); Data.Forest.Wild  = table();
Data.Shrub.Urban  = table(); Data.Shrub.Wild   = table();

fprintf('Loading data...\n');
for i = 1:numel(lc)
    for j = 1:numel(fireType)
        for k = 1:numel(decade)
            fName = sprintf('dataFig/variable/%s_%s_%s.csv', lc{i}, fireType{j}, decade{k});
            if isfile(fName)
                tmp = readtable(fName);
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
    end
end

% ... (Keep Setup and Data Loading sections from previous code) ...

% --- 3. CALCULATE Z-SCORES & SIGNIFICANCE ---
stats = zeros(nVars, 4); % [Mean_F, SE_F, Mean_S, SE_S]
pVals = zeros(nVars, 2); % [p_Forest, p_Shrub]

for i = 1:nVars
    vName = varNames{i};
    
    % --- FOREST ---
    uVec = Data.Forest.Urban.(vName); uVec = uVec(~isnan(uVec));
    wVec = Data.Forest.Wild.(vName);  wVec = wVec(~isnan(wVec));
    
    % Baseline
    mu = mean(wVec);
    sigma = std(wVec);
    
    % Z-Scores
    zScores = (uVec - mu) / sigma;
    
    stats(i, 1) = mean(zScores);
    stats(i, 2) = std(zScores) / sqrt(length(zScores)); % SE
    
    % T-TEST (Urban vs Wildland)
    [~, p] = ttest2(uVec, wVec, 'Vartype','unequal');
    pVals(i, 1) = p;
    
    % --- SHRUB ---
    uVec = Data.Shrub.Urban.(vName); uVec = uVec(~isnan(uVec));
    wVec = Data.Shrub.Wild.(vName);  wVec = wVec(~isnan(wVec));
    
    mu = mean(wVec);
    sigma = std(wVec);
    
    zScores = (uVec - mu) / sigma;
    
    stats(i, 3) = mean(zScores);
    stats(i, 4) = std(zScores) / sqrt(length(zScores));
    
    [~, p] = ttest2(uVec, wVec, 'Vartype','unequal');
    pVals(i, 2) = p;
end

% --- 4. PLOTTING WITH TIERED STARS ---
figure('Position', [100, 100, 800, 400], 'Color', 'w');
hold on;

b = bar(stats(:, [1, 3]), 'grouped');

% Styling
b(1).FaceColor = cForest; b(1).EdgeColor = 'none';
b(2).FaceColor = cShrub;  b(2).EdgeColor = 'none';

% Add Error Bars & Stars
ngroups = nVars;
nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    % Calculate center X positions
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    
    meanVal = stats(:, (i-1)*2 + 1);
    seVal   = stats(:, (i-1)*2 + 2);
    pVal    = pVals(:, i);
    
    % Draw Error Bars
    errorbar(x, meanVal, seVal, 'k', 'linestyle', 'none', 'LineWidth', 1, 'CapSize', 8);
    
    % Draw Significance Stars
    for j = 1:ngroups
        starStr = '';
        if pVal(j) < 0.001
            starStr = '***';
        elseif pVal(j) < 0.01
            starStr = '**';
        elseif pVal(j) < 0.05
            starStr = '*';
        end
        
        if ~isempty(starStr)
            % Determine Y position: 
            % If positive bar: Above Mean + SE
            % If negative bar: Below Mean - SE (or keep above axis if preferred)
            
            offset = 0.05; % Space between error bar and star
            
            if meanVal(j) >= 0
                y_pos = meanVal(j) + seVal(j) + offset;
                va = 'bottom';
            else
                % Option A: Put star below negative bars
                y_pos = meanVal(j) - seVal(j) - offset;
                va = 'top'; 
                
                % Option B: If you prefer all stars on top regardless of bar direction
                % y_pos = max(0, meanVal(j) + seVal(j)) + offset;
                % va = 'bottom';
            end
            
            text(x(j), y_pos, starStr, 'FontSize', 12, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', va);
        end
    end
end

% --- FORMATTING ---
yline(0, 'k-', 'LineWidth', 1);
yline(0.5, ':', 'Color', [0.6 0.6 0.6]); 
yline(-0.5, ':', 'Color', [0.6 0.6 0.6]); 

ylabel('Deviation from Wildland (\sigma)', 'FontSize', 12);
title('Environmental Anomalies (*p<0.05, **p<0.01, ***p<0.001)', 'FontSize', 14);
xticks(1:nVars);
xticklabels(xLabels);
legend({'Forest', 'Shrub'}, 'Location', 'best');

box on;
hold off;
exportgraphics(gcf, 'Fig/ZScore_Stars.png');