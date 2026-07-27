function CVEV_mean = pls_CVEV(X, Y, ncomp, K)
    tic
    [N, p] = size(X);
    M = size(Y, 2);
    
    if size(Y, 1) ~= N
        error('X and Y have mismatched number of samples');
    end
    
    ncomp = min([ncomp, p, N-2]);
    
    [cv_folds, embargo_trains] = blockCV(N, 100, K, 0, 42);
    folds = struct('train_idx', cell(K,1), 'test_idx', cell(K,1));
    for k = 1:K
        folds(k).test_idx  = cv_folds{k};
        folds(k).train_idx = embargo_trains{k};
    end
    
    Y_pred_all = nan(N, M, ncomp);  % N × M × ncomp

    for k = 1:K
        tr = folds(k).train_idx(:);
        te = folds(k).test_idx(:);
        
        Xtr = X(tr, :);
        Xte = X(te, :);
        Ytr = Y(tr, :);
        
        xmean = mean(Xtr, 1);
        ymean = mean(Ytr, 1);
        ones_te = ones(size(Xte, 1), 1);
        
        [XL, YL, ~, ~, ~, ~, ~, stats] = plsregress(Xtr, Ytr, ncomp);
        P = XL;
        Q = YL;
        W = stats.W;

        for a = 1:ncomp
            Wm = W(:, 1:a);
            Pm = P(:, 1:a);
            Qm = Q(:, 1:a);
            
            A = (Pm' * Wm);
            if rcond(A) < 1e-10
                coef = pinv(A) * Qm';
            else
                coef = A \ Qm';
            end
            
            Bm = Wm * coef;
            beta0 = ymean - xmean * Bm;
            
            Yhat = [ones_te, Xte] * [beta0; Bm];
            
            Y_pred_all(te, :, a) = Yhat;
        end
    end
    
    %% ===== compute overall CVEV =====
    CVEV_mean = nan(ncomp, M);
    
    for a = 1:ncomp
        for j = 1:M
            y_true = Y(:, j);
            y_pred = Y_pred_all(:, j, a);
            
            valid_idx = ~isnan(y_pred);
            y_true_valid = y_true(valid_idx);
            y_pred_valid = y_pred(valid_idx);
            
            if isempty(y_true_valid) || var(y_true_valid, 1) == 0
                CVEV_mean(a, j) = NaN;
            else
                SS_res = sum((y_true_valid - y_pred_valid).^2);
                SS_tot = sum((y_true_valid - mean(y_true_valid)).^2);
                CVEV_mean(a, j) = 1 - SS_res / SS_tot;
            end
        end
    end
    
    toc

    fprintf('\n=== CVEV ===\n');
    fprintf('ncomp.| ');
    for j = 1:M
        fprintf('Out%d     ', j);
    end
    fprintf('| avg. \n');
    fprintf('-------|');
    for j = 1:M
        fprintf('---------|');
    end
    fprintf('----------\n');
    
    show_comps = unique([1, round(ncomp/2), ncomp]);
    for a = show_comps
        fprintf('  %2d   | ', a);
        for j = 1:M
            fprintf('%8.4f | ', CVEV_mean(a, j));
        end
        fprintf('%8.4f\n', mean(CVEV_mean(a, :)));
    end
end