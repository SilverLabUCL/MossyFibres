function behaviour_results = analyze_behaviours_MFnum(X_raw, L_state, valid_t, behaviours, behaviour_names, selection_mode, file)
    if nargin < 7
        selection_mode = 'sorted';
    end
    lambda = 0.01;
    n_fold = 5;
    repeats = 3;
    neuron_range = [1:4,5,10,15,20:10:120];
    trace_neurons = [1,5,10,30,60,120];
    block_size = 300;
    
    behaviour_results = struct();
    fprintf('Will analyze %d behaviours: %s\n', length(behaviours), strjoin(behaviour_names, ', '));
   
    [~, zc_wL] = bootstrap_cc(X_raw, L_state', 300, 200);
    for b = 1:length(behaviours)
        Y_all = behaviours{b};
        try
            Y = Y_all(valid_t)';
            X = X_raw(valid_t,:);
        catch
            Y = Y_all';
            X = X_raw(:,:);
        end
        [chopt1, chopt2] = chop_by_file(file, valid_t);
        Y = Y(chopt1:chopt2);
        X = X(chopt1:chopt2, :);
        Y = zscore(Y);
        num_trials = length(Y);
        num_neurons = size(X, 2);
        
        CVEV_curve = zeros(length(neuron_range), repeats);
        Ypred_trace = nan(length(trace_neurons), num_trials, repeats);
        Ytest_trace = nan(length(trace_neurons), num_trials, repeats);
        
        %% ======= Step 1: Univariate Ridge Sorting =======
        fprintf('  [Behaviour: %s] Sorting neurons by Univariate Ridge...\n', behaviour_names{b});
        CVEV_univar = zeros(num_neurons, 1);
        
        for i = 1:num_neurons
            Xi = X(:, i);
            CVEV_univar(i) = compute_cvev(Xi, Y, lambda, n_fold);
        end
        
        [sorted_CVEV, sorted_idx] = sort(CVEV_univar, 'descend');
        behaviour_results.(behaviour_names{b}).BestMFidx = sorted_idx;
        behaviour_results.(behaviour_names{b}).EachMFCVEV = sorted_CVEV;

        %% ======= Step 2: Record full-population info =======
        best_MF_num = num_neurons;
        behaviour_results.(behaviour_names{b}).OptMFnum = best_MF_num;
        selected_best_idx = 1:num_neurons;

        %% ======= Step 3: Full-population CVEV and block-shuffle control =======
        X_best = X;
        behaviour_results.(behaviour_names{b}).OptMFCVEV = compute_cvev(X_best, Y, lambda, n_fold);

        % Block-shuffle control: shuffle Y in blocks to destroy fine temporal structure
        n_shuf = 20;
        CVEV_shuf_all = nan(n_shuf, 1);
        for s = 1:n_shuf
            Y_sh = block_shuffle(Y, block_size);
            CVEV_shuf_all(s) = compute_cvev(X_best, Y_sh, lambda, n_fold);
        end
        
        behaviour_results.(behaviour_names{b}).OptMF_ShuffleAcross      = mean(CVEV_shuf_all);
        behaviour_results.(behaviour_names{b}).OptMF_ShuffleAcross_std  = std(CVEV_shuf_all);
        behaviour_results.(behaviour_names{b}).OptMF_ShuffleAcross_all  = CVEV_shuf_all;
        
        %% ======= Step 4: Decoding curve over increasing neuron counts =======
        for r = 1:repeats
            fprintf('  Repeat %d/%d...\n', r, repeats);
            for n = 1:length(neuron_range)
                k = neuron_range(n);
        
                if strcmp(selection_mode, 'random')
                    selected_idx = randperm(num_neurons, k);
                elseif strcmp(selection_mode, 'sorted')
                    selected_idx = sorted_idx(1:k);
                else
                    error('Unknown selection_mode');
                end
        
                X_sub = X(:, selected_idx);
        
                if ismember(k, trace_neurons)
                    [cvev_k, y_pred_all, y_test_all] = compute_cvev(X_sub, Y, lambda, n_fold);
                else
                    cvev_k = compute_cvev(X_sub, Y, lambda, n_fold);
                end
        
                CVEV_curve(n, r) = cvev_k;
        
                if ismember(k, trace_neurons)
                    idx_tr = find(trace_neurons == k);
                    Ypred_trace(idx_tr, :, r) = y_pred_all;
                    Ytest_trace(idx_tr, :, r) = y_test_all;
                end
            end
        end

        mean_CVEV = mean(CVEV_curve, 2);
        sem_CVEV = std(CVEV_curve, 0, 2) / sqrt(repeats);
        behaviour_results.(behaviour_names{b}).mean_CVEV = mean_CVEV;
        behaviour_results.(behaviour_names{b}).sem_CVEV  = sem_CVEV;
        behaviour_results.(behaviour_names{b}).neuron_range = neuron_range;
        behaviour_results.(behaviour_names{b}).trace_neurons = trace_neurons;
        behaviour_results.(behaviour_names{b}).Ypred_trace = mean(Ypred_trace, 3);
        behaviour_results.(behaviour_names{b}).Ytest_trace = mean(Ytest_trace, 3);
        
        %% ======= Step 5: Decoding by neuron property (PM / NM / NSM) =======
        fprintf('  Analyzing pm/nm/nsm subsets...\n');
        neuron_counts = [1:4, 5, 10, 15, 20, 25, 30];
        neuron_counts = neuron_counts(neuron_counts <= size(X, 2));
        behaviour_results.(behaviour_names{b}).property.neuron_counts = neuron_counts;
        
        beta_sort = CVEV_univar; 
        
        % Classify neurons by bootstrap z-score of correlation with L_state
        sig_lev = 2;
        pm_idx = find(zc_wL > sig_lev);
        nm_idx = find(zc_wL < -sig_lev);
        nsm_idx = find(abs(zc_wL) < sig_lev);
        
        prop_types = {'pm', 'nm', 'nsm'};
        prop_groups = {pm_idx, nm_idx, nsm_idx};
        
        for p = 1:length(prop_types)
            group = prop_groups{p};
            type_name = prop_types{p};
            
            % Sort neurons within this group by univariate CVEV
            [~, idx_sort] = sort(beta_sort(group), 'descend');
            group_sorted = group(idx_sort);
            
            CVEV_prop = nan(length(neuron_counts), 1);
            CVEV_shuf = nan(length(neuron_counts), 1);
            
            for idx_k = 1:length(neuron_counts)
                k = neuron_counts(idx_k);
                if k > length(group_sorted)
                    continue;
                end
                sub_idx = group_sorted(1:k);
                X_sub = X(:, sub_idx);
                CVEV_prop(idx_k) = compute_cvev(X_sub, Y, lambda, n_fold);
                Y_shuf = block_shuffle(Y, block_size);
                CVEV_shuf(idx_k) = compute_cvev(X_sub, Y_shuf, lambda, n_fold);
            end
            behaviour_results.(behaviour_names{b}).property.(type_name).CVEV = CVEV_prop;
            behaviour_results.(behaviour_names{b}).property.([type_name '_shuf']).CVEV = CVEV_shuf;
        end
        
        % Random neuron selection as baseline control
        rand_repeat = 5;
        rand_idx_all = find(CVEV_univar > -inf);
        
        rand_CVEV = nan(length(neuron_counts), 1);
        rand_CVEV_shuf = nan(length(neuron_counts), 1);
        
        for idx_k = 1:length(neuron_counts)
            k = neuron_counts(idx_k);
            if k > length(rand_idx_all)
                continue;
            end
        
            vals = nan(rand_repeat, 1);
            vals_shuf = nan(rand_repeat, 1);
        
            for r = 1:rand_repeat
                sel_idx = randsample(rand_idx_all, k);
                X_sub = X(:, sel_idx);
        
                vals(r) = compute_cvev(X_sub, Y, lambda, n_fold);
        
                Y_shuf = block_shuffle(Y, block_size);
                vals_shuf(r) = compute_cvev(X_sub, Y_shuf, lambda, n_fold);
            end
        
            rand_CVEV(idx_k) = mean(vals);
            rand_CVEV_shuf(idx_k) = mean(vals_shuf);
        end
        
        behaviour_results.(behaviour_names{b}).property.random.CVEV = rand_CVEV;
        behaviour_results.(behaviour_names{b}).property.random_shuf.CVEV = rand_CVEV_shuf;
    end
end

%% ======================== Helper functions ========================

function [cvev, y_pred_all, y_test_all] = compute_cvev(X, Y, lambda, n_fold)
% COMPUTE_CVEV  Cross-validated explained variance via ridge regression.
%   Accumulates SS_res and SS_tot across folds and returns 1 - SS_res/SS_tot.
%   Optionally returns full predicted and observed vectors when called with
%   two or more output arguments.

    num_trials = length(Y);
    gap_blocks = 0;
    [cv_folds, ~] = blockCV(num_trials, 100, n_fold, gap_blocks);
    
    want_traj  = nargout >= 2;
    if want_traj
        y_pred_all = nan(num_trials, 1);
        y_test_all = nan(num_trials, 1);
    else
        y_pred_all = []; y_test_all = [];
    end
    
    SS_res_total = 0;
    SS_tot_total = 0;
    
    for i = 1:n_fold
        te = cv_folds{i};
        tr = setdiff(1:num_trials, te);
        X_train = X(tr, :);
        X_test  = X(te, :);
        y_train = Y(tr);
        y_test  = Y(te);
        
        % ridge with scaled=0 returns [intercept; betas]
        beta   = ridge(y_train, X_train, lambda, 0);
        y_pred = [ones(numel(te), 1), X_test] * beta;

        SS_res = sum((y_test - y_pred).^2);
        mu_tr  = mean(y_train);
        SS_tot = sum((y_test - mu_tr).^2);
        
        SS_res_total = SS_res_total + SS_res;
        SS_tot_total = SS_tot_total + SS_tot;
        
        if want_traj
            y_pred_all(te) = y_pred;
            y_test_all(te) = y_test;
        end
    end
    
    if SS_tot_total < eps
        cvev = 0;
    else
        cvev = 1 - SS_res_total / SS_tot_total;
    end
end

function y_shuf = block_shuffle(y, block_size)
% BLOCK_SHUFFLE  Permute contiguous blocks of y to destroy fine temporal
%   structure while preserving within-block autocorrelation.
%   Residual samples beyond the last full block are kept in place.

    num_blocks = floor(length(y) / block_size);
    if num_blocks == 0
        y_shuf = y(randperm(length(y)));
        return;
    end
    idx = reshape(1:num_blocks*block_size, block_size, []);
    idx = idx(:, randperm(num_blocks));
    y_shuf = y(idx(:));
    if num_blocks*block_size < length(y)
        y_shuf = [y_shuf; y(num_blocks*block_size+1:end)];
    end
end