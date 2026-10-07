%% beta for 4-datasets, 2-fireTypes (OLS vs MLE + Bootstrap + Risk Prob)
clear; clc;

dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};   nSrc  = numel(dataSrc);
fireType = {'Urban-edge','Wildland'};               nFire = numel(fireType);
edges    = 10 .^ (-1:0.05:5);                       % Bins for visualization
nBoot    = 1000;                                    % Bootstrap iterations

dirs = {'dataFig/curve_OLS', 'dataFig/beta_OLS', ...
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

% MLE Storage
count_MLE    = NaN(nSrc, nFire);
area_MLE     = NaN(nSrc, nFire);
xmin_MLE     = NaN(nSrc, nFire);
beta_MLE     = NaN(nSrc, nFire); % Slope (Negative Alpha)
betaErr_MLE  = NaN(nSrc, nFire);

prob_500_MLE  = NaN(nSrc, nFire); % P(Size >= 500)
prob_5000_MLE = NaN(nSrc, nFire); % P(Size >= 5000)

res_Bootstrap = table();
for i = 1:nSrc
    fprintf('Processing %s...\n', dataSrc{i});
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    SizeVec = [data.size]';
    FireVec = {data.FireType}';
    
    obs_mle_alphas = NaN(1, nFire);     
    for j = 1:nFire
        sz = SizeVec(strcmp(FireVec, fireType{j}));
        if isempty(sz), continue; end
        
        [counts, binEdges] = histcounts(sz, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
        
        valid = probDensity > 0;
        x_vis = log10(binCenters(valid))'; 
        y_vis = log10(probDensity(valid))';
        
        % Save raw data
        T_curve = table(binCenters(valid)', probDensity(valid)', 'VariableNames', {'binCenter', 'probDensity'});
        writetable(T_curve, sprintf('dataFig/curve_OLS/%s-%s.csv', dataSrc{i}, fireType{j}));
        writetable(T_curve, sprintf('dataFig/curve_MLE/%s-%s.csv', dataSrc{i}, fireType{j}));
        
        % =================================================================
        % METHOD 1: OLS
        % =================================================================
        mdl = fitlm(x_vis, y_vis);
        
        x_fit_ols = linspace(min(x_vis), max(x_vis), 100)';
        y_fit_ols = predict(mdl, x_fit_ols);
        writetable(table(10.^(x_fit_ols), 10.^(y_fit_ols), 'VariableNames', {'x_fit', 'y_fit'}), ...
            sprintf('dataFig/curve_OLS/%s-%s-fit.csv', dataSrc{i}, fireType{j}));
        
        count_OLS(i, j)   = numel(sz);
        area_OLS(i, j)    = sum(sz);
        alfa_OLS(i, j)    = 10^mdl.Coefficients.Estimate(1);
        beta_OLS(i, j)    = mdl.Coefficients.Estimate(2);
        betaErr_OLS(i, j) = mdl.Coefficients.SE(2);
        pBeta_OLS(i, j)   = mdl.Coefficients.pValue(2);
        R2_OLS(i, j)      = mdl.Rsquared.Ordinary;            
        
        % =================================================================
        % METHOD 2: MLE + PROBABILITY CALCULATION
        % =================================================================
        xmin = min(sz); 
        n = numel(sz);
        
        % MLE Calculation: Alpha = 1 + n / sum(ln(x/xmin))
        alpha_val = 1 + n / sum(log(sz ./ xmin));
        sigma_val = (alpha_val - 1) / sqrt(n); 
        
        obs_mle_alphas(j) = alpha_val; 
        
        % --- NEW: CALCULATE PROBABILITIES (RISK) ---
        % Formula: P(X >= x) = (x / xmin) ^ -(alpha - 1)
        % This gives the probability that a fire will be AT LEAST size X
        if 500 >= xmin
            prob_500_MLE(i, j) = (500 / xmin) ^ -(alpha_val - 1);
        else
            prob_500_MLE(i, j) = 1; % If cutoff > 500, prob is 1 (technically undefined below xmin)
        end
        
        if 5000 >= xmin
            prob_5000_MLE(i, j) = (5000 / xmin) ^ -(alpha_val - 1);
        else
            prob_5000_MLE(i, j) = 1; 
        end
        
        % MLE Fit Line
        x_fit_mle = linspace(min(sz), max(sz), 100)';
        term1 = (alpha_val - 1) / xmin;
        term2 = (x_fit_mle ./ xmin) .^ (-alpha_val);
        y_fit_mle = term1 .* term2;
        
        writetable(table(x_fit_mle, y_fit_mle, 'VariableNames', {'x_fit', 'y_fit'}), ...
            sprintf('dataFig/curve_MLE/%s-%s-fit.csv', dataSrc{i}, fireType{j}));
        
        % Store Stats
        count_MLE(i, j)   = n;
        area_MLE(i, j)    = sum(sz);
        xmin_MLE(i, j)    = xmin;
        beta_MLE(i, j)    = -alpha_val; 
        betaErr_MLE(i, j) = sigma_val;
    end
    
    % =================================================================
    % SIGNIFICANCE TEST: BOOTSTRAP
    % =================================================================
    idx_urb = find(strcmp(fireType, 'Urban-edge'));
    idx_wld = find(strcmp(fireType, 'Wildland'));
    sz_urb = SizeVec(strcmp(FireVec, 'Urban-edge'));
    sz_wld = SizeVec(strcmp(FireVec, 'Wildland'));
    
    if ~isempty(sz_urb) && ~isempty(sz_wld)
        boot_diffs = zeros(nBoot, 1);
        
        parfor b = 1:nBoot
            samp_u = datasample(sz_urb, numel(sz_urb));
            samp_w = datasample(sz_wld, numel(sz_wld));
            
            xm_u = min(samp_u); 
            a_u = 1 + numel(samp_u) / sum(log(samp_u ./ xm_u));
            
            xm_w = min(samp_w); 
            a_w = 1 + numel(samp_w) / sum(log(samp_w ./ xm_w));
            
            boot_diffs(b) = a_w - a_u; 
        end
        
        obs_diff = obs_mle_alphas(idx_wld) - obs_mle_alphas(idx_urb);
        ci_low  = prctile(boot_diffs, 2.5);
        ci_high = prctile(boot_diffs, 97.5);
        
        if obs_diff > 0
            p_val = mean(boot_diffs <= 0);
        else
            p_val = mean(boot_diffs >= 0);
        end
        
        row = table(dataSrc(i), -obs_mle_alphas(idx_urb), -obs_mle_alphas(idx_wld), ...
            obs_diff, ci_low, ci_high, p_val, ...
            'VariableNames', {'Dataset', 'Beta_MLE_Urban', 'Beta_MLE_Wild', ...
            'Alpha_Diff', 'CI_2_5', 'CI_97_5', 'P_Value'});
        res_Bootstrap = [res_Bootstrap; row];
    end
end

% save results from OLS Summary
for i = 1:nSrc
    writetable(table(count_OLS(i,:)', area_OLS(i,:)', alfa_OLS(i,:)', beta_OLS(i,:)', betaErr_OLS(i,:)', pBeta_OLS(i,:)', R2_OLS(i,:)', ...
        'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
        'RowNames', fireType), sprintf('dataFig/beta_OLS/2Fires-%s.csv', dataSrc{i}), 'WriteRowNames', true);
end

% save results from MLE Summary
for i = 1:nSrc
    writetable(table(count_MLE(i,:)', area_MLE(i,:)', xmin_MLE(i,:)', beta_MLE(i,:)', betaErr_MLE(i,:)', ...
        prob_500_MLE(i,:)', prob_5000_MLE(i,:)', ...
        'VariableNames', {'count', 'area', 'xmin', 'beta', 'betaErr', 'Prob_Ge_500', 'Prob_Ge_5000'}, ...
        'RowNames', fireType), sprintf('dataFig/beta_MLE/2Fires-%s.csv', dataSrc{i}), 'WriteRowNames', true);
end

% save results from Bootstrap Significance Summary
writetable(res_Bootstrap, 'dataFig/beta_MLE/Significance_Summary.csv');
disp('Processing Complete. Risk probabilities added to MLE summary.');