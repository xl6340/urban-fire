%% beta for 4-datasets, 4-decades
clear; clc;

dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};   nSrc  = numel(dataSrc);
decades  = {'1990s', '2000s', '2010s', '2020s'};    nDec  = numel(decades);

count    = NaN(nSrc, nDec);
area     = NaN(nSrc, nDec);
alfa     = NaN(nSrc, nDec);
beta     = NaN(nSrc, nDec);
betaErr  = NaN(nSrc, nDec);
pBeta    = NaN(nSrc, nDec);
R2       = NaN(nSrc, nDec);

edges = 10 .^ (-1:0.05:5); 
for i = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    Size = [data.size]';
    Dec  = {data.decade}';
    for d = 1:nDec
        sz = Size(strcmp(Dec, decades{d}));
        if isempty(sz), continue; end
        [counts, binEdges] = histcounts(sz, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
    
        x = log10(binCenters(probDensity > 0)); 
        y = log10(probDensity(probDensity > 0));
        writetable(table(binCenters', probDensity', 'VariableNames', {'binCenter', 'probDensity'}), ...
            sprintf('dataFig/curve/%s-%s.csv', dataSrc{i}, decades{d}));

        mdl = fitlm(x, y);
        x_fit = linspace(min(x), max(x), 100)';
        y_fit = predict(mdl, x_fit);
        writetable(table(10.^(x_fit), 10.^(y_fit), 'VariableNames', {'x_fit', 'y_fit'}), ...
            sprintf('dataFig/curve/%s-%s-fit.csv', dataSrc{i}, decades{d}));

        count(i, d)   = numel(sz);
        area(i, d)    = sum(sz);
        alfa(i, d)    = mdl.Coefficients.Estimate(1);
        beta(i, d)    = mdl.Coefficients.Estimate(2);
        betaErr(i, d) = mdl.Coefficients.SE(2);
        pBeta(i, d)   = mdl.Coefficients.pValue(2);
        R2(i, d)      = mdl.Rsquared.Ordinary;
    end
end

for i = 1:nSrc
        countVec = squeeze(count(i, :))';
        areaVec = squeeze(area(i, :))';
        alfaVec = squeeze(alfa(i, :))';
        betaVec = squeeze(beta(i, :))'; 
        betaErrVec = squeeze(betaErr(i, :))';
        pBetaVec = squeeze(pBeta(i, :))';
        R2Vec = squeeze(R2(i, :))';

        writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
            'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
            'RowNames', decades), sprintf('dataFig/beta/4Decs-%s.csv', dataSrc{i}), 'WriteRowNames', true);
end
%% beta for 4-datasets, 2-fireTypes (OLS vs MLE + Bootstrap + Risk Prob)
clear; clc;
dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};   nSrc  = numel(dataSrc);
fireType = {'Urban-edge','Wildland'};               nFire = numel(fireType);
edges    = 10 .^ (-1:0.05:5);                       % Bins for visualization
nBoot    = 1000;                                    % Bootstrap iterations
dirs = {'dataFig/curve', 'dataFig/beta', ...
        'dataFig/curve_MLE', 'dataFig/beta_MLE'};
for k = 1:numel(dirs)
    if ~exist(dirs{k}, 'dir'), mkdir(dirs{k}); end
end

% OLS Storage
count_OLS    = NaN(nSrc, nFire);
area_OLS     = NaN(nSrc, nFire);
alfa_OLS     = NaN(nSrc, nFire); 
beta_OLS     = NaN(nSrc, nFire); 
betaErr_OLS  = NaN(nSrc, nFire);
pBeta_OLS    = NaN(nSrc, nFire);
R2_OLS       = NaN(nSrc, nFire);
prob_500_OLS = NaN(nSrc, nFire); 
prob_5000_OLS = NaN(nSrc, nFire);

% MLE Storage
count_MLE    = NaN(nSrc, nFire);
area_MLE     = NaN(nSrc, nFire);
xmin_MLE     = NaN(nSrc, nFire);
beta_MLE     = NaN(nSrc, nFire); 
betaErr_MLE  = NaN(nSrc, nFire);
prob_500_MLE  = NaN(nSrc, nFire); 
prob_5000_MLE = NaN(nSrc, nFire);

% Bootstrap Result Tables
res_Bootstrap_MLE = table();
res_Bootstrap_OLS = table(); % <--- NEW TABLE FOR OLS

for i = 1:nSrc
    fprintf('Processing %s...\n', dataSrc{i});
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    SizeVec = [data.size]';
    FireVec = {data.FireType}';
    
    % Store observed values for bootstrap comparison later
    obs_mle_alphas = NaN(1, nFire);     
    obs_ols_betas  = NaN(1, nFire);
    
    for j = 1:nFire
        sz = SizeVec(strcmp(FireVec, fireType{j}));
        if isempty(sz), continue; end
        
        % --- PREPARE DATA (BINNING) ---
        [counts, binEdges] = histcounts(sz, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
        
        valid = probDensity > 0;
        x_vis = log10(binCenters(valid))'; 
        y_vis = log10(probDensity(valid))';
        
        % Save raw data
        T_curve = table(binCenters(valid)', probDensity(valid)', 'VariableNames', {'binCenter', 'probDensity'});
        writetable(T_curve, sprintf('dataFig/curve/%s-%s.csv', dataSrc{i}, fireType{j}));
        writetable(T_curve, sprintf('dataFig/curve_MLE/%s-%s.csv', dataSrc{i}, fireType{j}));
        
        % =================================================================
        % METHOD 1: OLS
        % =================================================================
        mdl = fitlm(x_vis, y_vis);
        
        intercept = mdl.Coefficients.Estimate(1); 
        slope     = mdl.Coefficients.Estimate(2); 
        obs_ols_betas(j) = slope; 
        
        x_fit_ols = linspace(min(x_vis), max(x_vis), 100)';
        y_fit_ols = predict(mdl, x_fit_ols);
        writetable(table(10.^(x_fit_ols), 10.^(y_fit_ols), 'VariableNames', {'x_fit', 'y_fit'}), ...
            sprintf('dataFig/curve/%s-%s-fit.csv', dataSrc{i}, fireType{j}));
        
        count_OLS(i, j)   = numel(sz);
        area_OLS(i, j)    = sum(sz);
        alfa_OLS(i, j)    = 10^intercept; 
        beta_OLS(i, j)    = slope;
        betaErr_OLS(i, j) = mdl.Coefficients.SE(2);
        pBeta_OLS(i, j)   = mdl.Coefficients.pValue(2);
        R2_OLS(i, j)      = mdl.Rsquared.Ordinary;            
        
        % OLS PDF Calculation
        prob_500_OLS(i, j)  = 10^(intercept + slope * log10(500));
        prob_5000_OLS(i, j) = 10^(intercept + slope * log10(5000));
        
        % =================================================================
        % METHOD 2: MLE
        % =================================================================
        xmin = min(sz); 
        n = numel(sz);
        
        alpha_val = 1 + n / sum(log(sz ./ xmin));
        sigma_val = (alpha_val - 1) / sqrt(n); 
        obs_mle_alphas(j) = alpha_val; 
        
        % MLE PDF Calculation
        if 500 >= xmin
            prob_500_MLE(i, j) = ((alpha_val - 1) / xmin) * (500 / xmin) ^ (-alpha_val);
        else, prob_500_MLE(i, j) = 0; end
        
        if 5000 >= xmin
            prob_5000_MLE(i, j) = ((alpha_val - 1) / xmin) * (5000 / xmin) ^ (-alpha_val);
        else, prob_5000_MLE(i, j) = 0; end
        
        % MLE Fit Line
        x_fit_mle = linspace(min(sz), max(sz), 100)';
        term1 = (alpha_val - 1) / xmin;
        term2 = (x_fit_mle ./ xmin) .^ (-alpha_val);
        y_fit_mle = term1 .* term2;
        writetable(table(x_fit_mle, y_fit_mle, 'VariableNames', {'x_fit', 'y_fit'}), ...
            sprintf('dataFig/curve_MLE/%s-%s-fit.csv', dataSrc{i}, fireType{j}));
        
        count_MLE(i, j)   = n;
        area_MLE(i, j)    = sum(sz);
        xmin_MLE(i, j)    = xmin;
        beta_MLE(i, j)    = -alpha_val; 
        betaErr_MLE(i, j) = sigma_val;
    end
    
    % =================================================================
    % SIGNIFICANCE TEST: BOOTSTRAP (MLE & OLS)
    % =================================================================
    idx_urb = find(strcmp(fireType, 'Urban-edge'));
    idx_wld = find(strcmp(fireType, 'Wildland'));
    sz_urb = SizeVec(strcmp(FireVec, 'Urban-edge'));
    sz_wld = SizeVec(strcmp(FireVec, 'Wildland'));
    
    if ~isempty(sz_urb) && ~isempty(sz_wld)
        boot_diffs_MLE = zeros(nBoot, 1);
        boot_diffs_OLS = zeros(nBoot, 1); 
        
        edges_local = edges; 
        
        parfor b = 1:nBoot
            % 1. Resample Data
            samp_u = datasample(sz_urb, numel(sz_urb));
            samp_w = datasample(sz_wld, numel(sz_wld));
            
            % --- MLE Bootstrap ---
            xm_u = min(samp_u); 
            a_u = 1 + numel(samp_u) / sum(log(samp_u ./ xm_u));
            xm_w = min(samp_w); 
            a_w = 1 + numel(samp_w) / sum(log(samp_w ./ xm_w));
            
            % MLE Difference (keeping your original alpha logic for MLE)
            boot_diffs_MLE(b) = a_w - a_u; 
            
            % --- OLS Bootstrap ---
            % Urban OLS
            [c_u, ~] = histcounts(samp_u, edges_local);
            pd_u = c_u ./ (sum(c_u) * diff(edges_local));
            bc_u = sqrt(edges_local(1:end-1) .* edges_local(2:end));
            valid_u = pd_u > 0;
            if sum(valid_u) > 1
                p_u = polyfit(log10(bc_u(valid_u)), log10(pd_u(valid_u)), 1);
                slope_u = p_u(1);
            else
                slope_u = NaN;
            end
            
            % Wildland OLS
            [c_w, ~] = histcounts(samp_w, edges_local);
            pd_w = c_w ./ (sum(c_w) * diff(edges_local));
            bc_w = sqrt(edges_local(1:end-1) .* edges_local(2:end));
            valid_w = pd_w > 0;
            if sum(valid_w) > 1
                p_w = polyfit(log10(bc_w(valid_w)), log10(pd_w(valid_w)), 1);
                slope_w = p_w(1);
            else
                slope_w = NaN;
            end
            
            % --- CORRECTION HERE: Urban - Wildland ---
            boot_diffs_OLS(b) = slope_u - slope_w; 
        end
        
        % --- Process MLE Results ---
        obs_diff_mle = obs_mle_alphas(idx_wld) - obs_mle_alphas(idx_urb);
        ci_low_mle  = prctile(boot_diffs_MLE, 2.5);
        ci_high_mle = prctile(boot_diffs_MLE, 97.5);
        if obs_diff_mle > 0, p_val_mle = mean(boot_diffs_MLE <= 0); else, p_val_mle = mean(boot_diffs_MLE >= 0); end
        
        row_mle = table(dataSrc(i), -obs_mle_alphas(idx_urb), -obs_mle_alphas(idx_wld), ...
            obs_diff_mle, ci_low_mle, ci_high_mle, p_val_mle, ...
            'VariableNames', {'Dataset', 'Beta_MLE_Urban', 'Beta_MLE_Wild', ...
            'Alpha_Diff', 'CI_2_5', 'CI_97_5', 'P_Value'});
        res_Bootstrap_MLE = [res_Bootstrap_MLE; row_mle];
        
        % --- Process OLS Results ---
        boot_diffs_OLS = boot_diffs_OLS(~isnan(boot_diffs_OLS)); 
        
        % --- CORRECTION HERE: Calculate observed diff as Urban - Wild ---
        obs_diff_ols = obs_ols_betas(idx_urb) - obs_ols_betas(idx_wld);
        
        ci_low_ols  = prctile(boot_diffs_OLS, 2.5);
        ci_high_ols = prctile(boot_diffs_OLS, 97.5);
        
        % P-value logic: If Diff > 0, we test how often bootstrap is <= 0
        if obs_diff_ols > 0
            p_val_ols = mean(boot_diffs_OLS <= 0);
        else
            p_val_ols = mean(boot_diffs_OLS >= 0);
        end
        
        row_ols = table(dataSrc(i), obs_ols_betas(idx_urb), obs_ols_betas(idx_wld), ...
            obs_diff_ols, ci_low_ols, ci_high_ols, p_val_ols, ...
            'VariableNames', {'Dataset', 'Beta_OLS_Urban', 'Beta_OLS_Wild', ...
            'Beta_Diff_Urb_Min_Wild', 'CI_2_5', 'CI_97_5', 'P_Value'}); % Renamed column for clarity
        res_Bootstrap_OLS = [res_Bootstrap_OLS; row_ols];
    end
end

% Save results OLS
for i = 1:nSrc
    writetable(table(count_OLS(i,:)', area_OLS(i,:)', alfa_OLS(i,:)', beta_OLS(i,:)', betaErr_OLS(i,:)', pBeta_OLS(i,:)', R2_OLS(i,:)', ...
        prob_500_OLS(i,:)', prob_5000_OLS(i,:)', ...
        'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2', 'PDF_at_500', 'PDF_at_5000'}, ...
        'RowNames', fireType), sprintf('dataFig/beta/2Fires-%s.csv', dataSrc{i}), 'WriteRowNames', true);
end

% Save results MLE
for i = 1:nSrc
    writetable(table(count_MLE(i,:)', area_MLE(i,:)', xmin_MLE(i,:)', beta_MLE(i,:)', betaErr_MLE(i,:)', ...
        prob_500_MLE(i,:)', prob_5000_MLE(i,:)', ...
        'VariableNames', {'count', 'area', 'xmin', 'beta', 'betaErr', 'PDF_at_500', 'PDF_at_5000'}, ...
        'RowNames', fireType), sprintf('dataFig/beta_MLE/2Fires-%s.csv', dataSrc{i}), 'WriteRowNames', true);
end

% Save Bootstrap Summaries
writetable(res_Bootstrap_MLE, 'dataFig/beta_MLE/Significance_Summary.csv');
writetable(res_Bootstrap_OLS, 'dataFig/beta/Significance_Summary.csv');
disp('Processing Complete. OLS and MLE Bootstrap Significance Calculated.');
%% beta for 4-datasets, 4-decades, 2-fireTypes
clear; clc;

dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};   nSrc  = numel(dataSrc);
decades  = {'1990s', '2000s', '2010s', '2020s'};    nDec  = numel(decades);
fireType = {'Urban-edge','Wildland'};               nFire = numel(fireType);
count    = NaN(nSrc, nDec, nFire);
area     = NaN(nSrc, nDec, nFire);
alfa     = NaN(nSrc, nDec, nFire);
beta     = NaN(nSrc, nDec, nFire);
betaErr  = NaN(nSrc, nDec, nFire);
pBeta    = NaN(nSrc, nDec, nFire);
R2       = NaN(nSrc, nDec, nFire);

edges = 10 .^ (-1:0.05:5); 
for i = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    Size = [data.size]';
    Fire = {data.FireType}';
    Dec  = {data.decade}';
    for d = 1:nDec
        for j = 1:nFire
            sz = Size(strcmp(Fire, fireType{j}) & strcmp(Dec, decades{d}));
            if isempty(sz), continue; end
            [counts, binEdges] = histcounts(sz, edges);
            binWidth = diff(binEdges);
            probDensity = counts ./ (sum(counts) * binWidth); 
            binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
        
            x = log10(binCenters(probDensity > 0)); 
            y = log10(probDensity(probDensity > 0));
            writetable(table(x', y', 'VariableNames', {'binCenter', 'probDensity'}), ...
                sprintf('dataFig/curve/%s-%s-%s.csv', dataSrc{i}, fireType{j}, decades{d}));

            mdl = fitlm(x, y);
            x_fit = linspace(min(x), max(x), 100)';
            y_fit = predict(mdl, x_fit);
            writetable(table(x_fit, y_fit, 'VariableNames', {'x_fit', 'y_fit'}), ...
                sprintf('dataFig/curve/%s-%s-%s-fit.csv', dataSrc{i}, fireType{j}, decades{d}));

            count(i, d, j)   = numel(sz);
            area(i, d, j)  = sum(sz);
            alfa(i, d, j)   = mdl.Coefficients.Estimate(1);
            beta(i, d, j)   = mdl.Coefficients.Estimate(2);
            betaErr(i, d, j) = mdl.Coefficients.SE(2);
            pBeta(i, d, j)   = mdl.Coefficients.pValue(2);
            R2(i, d, j)      = mdl.Rsquared.Ordinary;            
        end
    end
end

for i = 1:nSrc
    for j = 1:nFire
        countVec = squeeze(count(i, :, j))';
        areaVec = squeeze(area(i, :, j))';
        alfaVec = squeeze(alfa(i, :, j))';
        betaVec = squeeze(beta(i, :, j))'; 
        betaErrVec = squeeze(betaErr(i, :, j))';
        pBetaVec = squeeze(pBeta(i, :, j))';
        R2Vec = squeeze(R2(i, :, j))';

        writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
            'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
            'RowNames', decades), ...
            sprintf('dataFig/beta/4Decs-%s-%s.csv', dataSrc{i}, fireType{j}), "WriteRowNames", true);
    end
end
%% beta for CalFire, 2-fireTypes, 2-ignition
clear; clc;

fireType = {'Urban-edge','Wildland'}; nFire = numel(fireType);
igType   = {'Human', 'Natural'};      nIg   = numel(igType);

data = shaperead('dataPrc/firePrmt/CalFire.shp');
Size = [data.size]';
Fire = {data.FireType}';
Ign  = {data.Ignition}';

count    = NaN(nFire ,nIg);
area     = NaN(nFire ,nIg);
alfa     = NaN(nFire ,nIg);
beta     = NaN(nFire ,nIg);
betaErr  = NaN(nFire ,nIg);
pBeta    = NaN(nFire ,nIg);
R2       = NaN(nFire ,nIg);
edges = 10 .^ (0:0.05:5); 
for i = 1:nFire
    for j = 1:nIg
        sz = Size(strcmp(Fire, fireType{i}) & strcmp(Ign, igType{j}));
        if isempty(sz), continue; end
        [counts, binEdges] = histcounts(sz, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
    
        x = log10(binCenters(probDensity > 0)); 
        y = log10(probDensity(probDensity > 0));
        mdl = fitlm(x, y);

        count(i,j)   = numel(sz);
        area(i,j)    = sum(sz);
        alfa(i,j)    = mdl.Coefficients.Estimate(1);
        beta(i,j)    = mdl.Coefficients.Estimate(2);
        betaErr(i,j) = mdl.Coefficients.SE(2);
        pBeta(i,j)   = mdl.Coefficients.pValue(2);
        R2(i,j)      = mdl.Rsquared.Ordinary;            
    end
end

for i = 1:nFire
    countVec = count(i,:)';  
    areaVec = area(i,:)';    
    alfaVec = alfa(i,:)';    
    betaVec = beta(i,:)';    
    betaErrVec = betaErr(i,:)'; 
    pBetaVec = pBeta(i,:)';   
    R2Vec = R2(i,:)';         

    writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
        'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
        'RowNames', igType), sprintf('dataFig/beta/2Igns-%s.csv', fireType{i}), "WriteRowNames", true);
end
%% beta for CalFire, 4-decade, 2-ignition
clear; clc;
fireType = {'Urban-edge','Wildland'}; 
decades  = {'1990s', '2000s', '2010s', '2020s'};  nDec  = numel(decades);
igType   = {'Human', 'Natural'};                  nIg   = numel(igType);

data = shaperead('dataPrc/firePrmt/CalFire.shp');
Size = [data.size]';
Ign  = {data.Ignition}';
Dec  = {data.decade}';
Fire = {data.FireType}';

count    = NaN(nDec ,nIg);
area     = NaN(nDec ,nIg);
alfa     = NaN(nDec ,nIg);
beta     = NaN(nDec ,nIg);
betaErr  = NaN(nDec ,nIg);
pBeta    = NaN(nDec ,nIg);
R2       = NaN(nDec ,nIg);
edges = 10 .^ (0:0.05:5); 
for i = 1:nDec
    for j = 1:nIg
        sz = Size(strcmp(Dec, decades{i}) & strcmp(Ign, igType{j}) & strcmp(Fire, fireType{1}));
        if isempty(sz), continue; end
        [counts, binEdges] = histcounts(sz, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
    
        x = log10(binCenters(probDensity > 0)); 
        y = log10(probDensity(probDensity > 0));
        mdl = fitlm(x, y);

        count(i,j)   = numel(sz);
        area(i,j)    = sum(sz);
        alfa(i,j)    = mdl.Coefficients.Estimate(1);
        beta(i,j)    = mdl.Coefficients.Estimate(2);
        betaErr(i,j) = mdl.Coefficients.SE(2);
        pBeta(i,j)   = mdl.Coefficients.pValue(2);
        R2(i,j)      = mdl.Rsquared.Ordinary;            
    end
end

for i = 1:nDec
    countVec = count(i,:)';  
    areaVec = area(i,:)';    
    alfaVec = alfa(i,:)';    
    betaVec = beta(i,:)';    
    betaErrVec = betaErr(i,:)'; 
    pBetaVec = pBeta(i,:)';   
    R2Vec = R2(i,:)';         

    writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
        'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
        'RowNames', igType), sprintf('dataFig/beta/2Igns-%s.csv', decades{i}), "WriteRowNames", true);
end
%% beta for CalFire, 4 decades, 3 cut-off sizes + Probability + Bootstrap Comparison
clear; clc;

edgeCut  = {'<100 km^2', '<1000 km^2', 'Full size'}; 
nEdge    = numel(edgeCut);

edges = cell(nEdge,1); 
edges{1} = 10.^(-1:0.05:2); 
edges{2} = 10.^(-1:0.05:3); 
edges{3} = 10.^(-1:0.05:5); 

fireType = {'Urban-edge','Wildland'};               nFire = numel(fireType);
decades  = {'1990s', '2000s', '2010s', '2020s'};    nDec  = numel(decades);
nBoot    = 1000; % Bootstrap iterations

data = shaperead('dataPrc/firePrmt/CalFire.shp');
Size = [data.size]';
Fire = {data.FireType}';
Dec  = {data.decade}';

count    = NaN(nEdge, nFire, nDec);
area     = NaN(nEdge, nFire, nDec);
alfa     = NaN(nEdge, nFire, nDec);
beta     = NaN(nEdge, nFire, nDec);
betaErr  = NaN(nEdge, nFire, nDec);
pBeta    = NaN(nEdge, nFire, nDec);
R2       = NaN(nEdge, nFire, nDec);

prob_500  = NaN(nEdge, nFire, nDec); 
prob_5000 = NaN(nEdge, nFire, nDec);
res_Bootstrap = table();

for i = 1:nEdge
    current_edges = edges{i};
    max_val = max(current_edges); 
    
    for j = 1:nFire
        for d = 1:nDec
            sz = Size(strcmp(Fire, fireType{j}) & strcmp(Dec, decades{d}));
            sz = sz(sz <= max_val);             
            if numel(sz) < 10, continue; end            
           
            [counts, binEdges] = histcounts(sz, current_edges); 
            binWidth = diff(binEdges);
            probDensity = counts ./ (sum(counts) * binWidth); 
            binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
              
            valid = probDensity > 0;
            x = log10(binCenters(valid))'; 
            y = log10(probDensity(valid))';

            mdl = fitlm(x, y); 
            intercept = mdl.Coefficients.Estimate(1);
            slope     = mdl.Coefficients.Estimate(2);

            count(i,j,d)   = numel(sz);
            area(i,j,d)    = sum(sz);
            alfa(i,j,d)    = intercept;
            beta(i,j,d)    = slope;
            betaErr(i,j,d) = mdl.Coefficients.SE(2);
            pBeta(i,j,d)   = mdl.Coefficients.pValue(2);
            R2(i,j,d)      = mdl.Rsquared.Ordinary;       
            
            % Formula: P(x) = 10^(intercept + slope * log10(x))
            if max_val >= 500
                prob_500(i,j,d) = 10^(intercept + slope * log10(500));
            else
                prob_500(i,j,d) = 0; % Cannot predict 500km2 if max size is 100km2
            end
            
            if max_val >= 5000
                prob_5000(i,j,d) = 10^(intercept + slope * log10(5000));
            else
                prob_5000(i,j,d) = 0;
            end
        end    
    end
end

for j = 1:nFire
    for d = 1:nDec
        countVec = squeeze(count(:,j,d));  
        areaVec = squeeze(area(:,j,d));      
        alfaVec = squeeze(alfa(:,j,d));     
        betaVec = squeeze(beta(:,j,d));      
        betaErrVec = squeeze(betaErr(:,j,d));  
        pBetaVec = squeeze(pBeta(:,j,d));  
        R2Vec = squeeze(R2(:,j,d));
        
        prob500Vec = squeeze(prob_500(:,j,d));
        prob5000Vec = squeeze(prob_5000(:,j,d));
    
        writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
            prob500Vec, prob5000Vec, ...
            'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2', 'Prob_500', 'Prob_5000'}, ...
            'RowNames', edgeCut), sprintf('dataFig/beta/3cutSizes-%s-%s.csv',fireType{j}, decades{d}), 'WriteRowNames', true);
    end
end

for i = 1:nEdge
    current_edges = edges{i};
    max_val = max(current_edges);
    
    for d = 1:nDec
        sz_u = Size(strcmp(Fire, 'Urban-edge') & strcmp(Dec, decades{d}));
        sz_w = Size(strcmp(Fire, 'Wildland') & strcmp(Dec, decades{d}));
        
        sz_u = sz_u(sz_u <= max_val);
        sz_w = sz_w(sz_w <= max_val);
        
        if numel(sz_u) < 10 || numel(sz_w) < 10, continue; end
        
        beta_u_obs = beta(i, 1, d); 
        beta_w_obs = beta(i, 2, d);

        if isnan(beta_u_obs) || isnan(beta_w_obs), continue; end
        
        obs_diff = beta_u_obs - beta_w_obs;
        
        % Bootstrap Prep
        boot_diffs = zeros(nBoot, 1);
        binDiffs = diff(current_edges);
        % Calculate log centers once
        binCenters = sqrt(current_edges(1:end-1) .* current_edges(2:end));
        logBinCenters = log10(binCenters);
        
        parfor b = 1:nBoot
            samp_u = datasample(sz_u, numel(sz_u));
            samp_w = datasample(sz_w, numel(sz_w));
            
            % -- Fit Urban --
            [c_u, ~] = histcounts(samp_u, current_edges);
            pd_u = c_u ./ (sum(c_u) * binDiffs);
            valid_u = pd_u > 0;
            if sum(valid_u) > 1
                p_u = polyfit(logBinCenters(valid_u), log10(pd_u(valid_u)), 1);
                b_u = p_u(1);
            else
                b_u = NaN;
            end
            
            % -- Fit Wildland --
            [c_w, ~] = histcounts(samp_w, current_edges);
            pd_w = c_w ./ (sum(c_w) * binDiffs);
            valid_w = pd_w > 0;
            if sum(valid_w) > 1
                p_w = polyfit(logBinCenters(valid_w), log10(pd_w(valid_w)), 1);
                b_w = p_w(1);
            else, b_w = NaN; 
            end
            
            boot_diffs(b) = b_u - b_w;
        end
        
        % Process Bootstrap Results
        boot_diffs = boot_diffs(~isnan(boot_diffs));
        if isempty(boot_diffs), continue; end
        
        ci_low  = prctile(boot_diffs, 2.5);
        ci_high = prctile(boot_diffs, 97.5);
        
        if obs_diff > 0, p_val = mean(boot_diffs <= 0);
        else, p_val = mean(boot_diffs >= 0); end

        row = table(edgeCut(i), decades(d), beta_u_obs, beta_w_obs, obs_diff, ci_low, ci_high, p_val, ...
            'VariableNames', {'Cutoff', 'Decade', 'Beta_Urban', 'Beta_Wild', 'Beta_Diff', 'CI_2_5', 'CI_97_5', 'P_Value'});
        res_Bootstrap = [res_Bootstrap; row];
    end
end
writetable(res_Bootstrap, 'dataFig/beta/Significance_Summary_Cutoffs.csv');
%% climate, CalFire
varNames = {'tmean', 'vpdmax','tmin', 'tmax'}; nvar = numel(varNames);
fireType = {'Urban-edge','Wildland'};          nFire = numel(fireType);

data = shaperead('dataPrc/firePrmt/CalFire.shp');
Fire = {data.FireType}';

data_u = data(strcmp(Fire, fireType{1}));
data_w = data(strcmp(Fire, fireType{2}));

Var_u = NaN(numel(data_u), nvar);
Var_w = NaN(numel(data_w), nvar);
for i =1:nvar
    Var_u(:,i) = [data_u.(varNames{i})]';
    Var_w(:,i) = [data_w.(varNames{i})]';
end
writetable(array2table(Var_u, 'VariableNames', string(varNames)),'dataFig/climate/Urban-edge.csv');
writetable(array2table(Var_w, 'VariableNames', string(varNames)),'dataFig/climate/Wildland.csv');
%% NDVI, CalFire
ndviAll = [data.ndviM]';
for i = 1:nFire
    for d = 2:nDec
        ndvi = ndviAll(strcmp(Fire, fireType{i}) & strcmp(Dec, decades{d}) & ~isnan(ndviAll));
        writetable(table(ndvi), sprintf('dataFig/ndvi/%s-%s.csv', fireType{i}, decades{d}));
    end
end
%% year: fire number & burned area, CalFire
igType   = {'Human', 'Natural'};      nIg   = numel(igType);
years = 1990:2024;                    nyr = numel(years);
fireType = {'Urban-edge','Wildland'}; nFire = numel(fireType);

data = shaperead('dataPrc/firePrmt/CalFire_fullSize.shp');
Size = [data.size]';
Ign  = {data.Ignition}';
year = [data.year]';
Fire = {data.FireType}';  

count    = NaN(nIg, nFire, nyr);
area     = NaN(nIg, nFire, nyr);
for i = 1:nIg
    for j = 1:nFire
        for yr = 1:nyr
            sz = Size(strcmp(Ign, igType{i}) & strcmp(Fire, fireType{j}) & year == years(yr));        
            count(i,j,yr)   = numel(sz);
            area(i,j,yr)    = sum(sz);
        end
    end
end
for j = 1:nFire
    countVec = squeeze(count(:,j,:))'; 
    writetable(array2table([countVec, fracCount], 'VariableNames',{'Human', 'Natural', 'F_Human', 'F_Natur'}, ...
        'RowNames',string(years)), sprintf('dataFig/ba/count-%s.csv', fireType{j}), "WriteRowNames",true);

    areaVec = squeeze(area(:,j,:))'; 
    writetable(array2table([areaVec, fracArea], 'VariableNames', {'Human', 'Natural', 'F_Human', 'F_Natur'}, ...
        'RowNames',string(years)), sprintf('dataFig/ba/area-%s.csv', fireType{j}), "WriteRowNames",true);
end

frac_Count_Human = zeros(nFire, 1);
frac_Count_Natur = zeros(nFire, 1);
frac_Area_Human  = zeros(nFire, 1);
frac_Area_Natur  = zeros(nFire, 1);

for j = 1:nFire
    % Count Totals
    countVec = squeeze(count(:,j,:))'; 
    total_Count_Human = sum(countVec(:,1), 'omitnan');
    total_Count_Natur = sum(countVec(:,2), 'omitnan');
    grand_Total_Count = total_Count_Human + total_Count_Natur;
    
    frac_Count_Human(j) = total_Count_Human / grand_Total_Count;
    frac_Count_Natur(j) = total_Count_Natur / grand_Total_Count;
    
    % Area Totals
    areaVec = squeeze(area(:,j,:))';
    total_Area_Human = sum(areaVec(:,1), 'omitnan');
    total_Area_Natur = sum(areaVec(:,2), 'omitnan');
    grand_Total_Area = total_Area_Human + total_Area_Natur;
    
    frac_Area_Human(j) = total_Area_Human / grand_Total_Area;
    frac_Area_Natur(j) = total_Area_Natur / grand_Total_Area;
end

results = table(frac_Count_Human, frac_Count_Natur, frac_Area_Human, frac_Area_Natur, ...
    'RowNames', fireType', ...
    'VariableNames', {'Count_Frac_Human', 'Count_Frac_Natur', 'Area_Frac_Human', 'Area_Frac_Natur'});

writetable(results, 'dataFig/ba/fractions.csv', "WriteRowNames", true);
%% urban vs. wildland, forest vs shrubGrass
% NDVI, climate, terrain
clear; clc;

lc       = {'Forest', 'ShrubGrass'};                    nlc   = numel(lc);
fireType = {'Urban-edge','Wildland'};                   nFire = numel(fireType);
varName  = {'ndviM', 'elevation', 'slope', 'vs', 'vd', ... 
            'ppt', 'tmean', 'vpdmax','tmin', 'tmax', 'ignition'};   nvar  = numel(varName);

data     = shaperead('dataPrc/firePrmt/CalFire.shp');
lcs       = {data.lc}';
Fires     = {data.FireType}';

allData = zeros(numel(data), nvar);
for v = 1:nvar-1
    allData(:, v) = [data.(varName{v})]';
end

Ignitions = {data.Ignition}'; 
ignitionVals(strcmpi(Ignitions, 'Human')) = 1;
ignitionVals(strcmpi(Ignitions, 'Natural')) = 0;
ignitionVals(strcmpi(Ignitions, 'Unknown')) = nan;
allData(:,nvar) = ignitionVals;

summaryCount = {};
for i = 1:nlc
    for j = 1:nFire
            idx = strcmp(lcs, lc{i}) & strcmp(Fires, fireType{j});            
            if ~any(idx), continue; end
            subsetData = allData(idx, :);

            T = array2table(subsetData, 'VariableNames', string(varName));
            writetable(T, sprintf('dataFig/variable/%s_%s.csv', lc{i}, fireType{j}));
    end
end
%% forest: elevatiion stratification, climate, beta
clear; clc;

data = shaperead('dataPrc/firePrmt/CalFire.shp');

lcs        = {data.lc}';
Fires      = {data.FireType}';
elevations = [data.elevation]';
FireSizes  = [data.size]'; 
Tmeans     = [data.tmean]'; 
Tmaxs      = [data.tmax]';
Tmins      = [data.tmin]';
VPDmaxs    = [data.vpdmax]' / 10; % Convert VPD from hPa to kPa

lc        = {'Forest', 'Shrub', 'Grass'};                                         
fireType  = {'Urban-edge','Wildland'};              
nFire     = numel(fireType);
nBins     = 8; 

% Define Elevation Bins
eleEdges = unique(quantile(elevations(strcmp(lcs, lc{1})), linspace(0, 1, nBins+1)));

beta     = NaN(nFire, nBins);
betaErr  = NaN(nFire, nBins);
envMeans = NaN(nFire, nBins, 4); % 1=Tmean, 2=Tmax, 3=Tmin, 4=VPDmax
envErrs  = NaN(nFire, nBins, 4);

sizeEdges = 10 .^ (0:0.05:5); 
for i = 1:nFire
    for j = 1:nBins
        minEle = eleEdges(j);
        maxEle = eleEdges(j+1);
        
        idx = strcmp(Fires, fireType{i}) & strcmp(lcs, lc{1}) & ...
              elevations >= minEle & elevations < maxEle;          
        
        sz = FireSizes(idx);
        if numel(sz) < 10, continue; end
        
        % Calculate Beta
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

% export data
eleCenters = (eleEdges(1:end-1) + eleEdges(2:end)) / 2;
writematrix(eleCenters', 'dataFig/elevation/eleCenters.csv');
writematrix(elevations(strcmp(lcs, lc{1})), 'dataFig/elevation/elevations.csv');

valMean = cat(3, beta, envMeans);
valErr  = cat(3, betaErr, envErrs);
variable = {'beta', 'VPDmax', 'Tmean', 'Tmin', 'Tmax'};
for v = 1:5 
    dataMean = squeeze(valMean(:,:,v))'; 
    dataErr  = squeeze(valErr(:,:,v))';     
    T = array2table([dataMean, dataErr], ...
        'VariableNames', {'Mean_Urban-edge', 'Mean_Wildland', 'Err_Urban-edge', 'Err_Wildland'});
    writetable(T, sprintf('dataFig/elevation/%s.csv', variable{v}));
end
