results_dir = 'X:\MFB\MFB_AH_2023\PLSresults';
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

folder_names = {
    '171212_16_19_37';
    '191018_13_39_41';
    %'191018_13_56_55'; % remove loco, state (only small portion at the begin)
    '191018_14_30_00';
    %'191018_14_11_33'; % remove loco, state (only small portion at the end)
    '191209_13_44_12';
    %'191209_14_04_14'; % remove loco, state (only small portion at the end)
    %'191209_14_32_39'; % remove loco
    %'191209_15_01_22'; % remove loco
    '191209_14_18_13';
    '191209_14_46_58';
    '200130_13_21_13 FunctAcq';
    '200130_13_36_14 FunctAcq';
    '200130_13_49_09 FunctAcq';
    '200130_14_15_24 FunctAcq';
    '200130_14_29_30 FunctAcq';
};

n_components_all = 1:30;
K                = 5;      % K-fold CV
block_size       = 300;    % block shuffle
n_shuffles       = 5;
N_sub            = 100;
n_repeats        = 5;

for file_i = 1:length(folder_names)
    file = char(folder_names(file_i));
    fprintf('\n=== [%d/%d] Processing %s ===\n', file_i, length(folder_names), file);

    quickAnalysis

    main_path = fullfile('X:\MFB\MFB_AH_2023\Data', file, [file '_wsk_analysis.mat']);
    if exist(main_path, 'file')
        load(main_path);
    else
        alt_path = fullfile('X:\MFB\MFB_AH_2023\Data', file, 'Functional_data', [file '_wsk_analysis.mat']);
        if exist(alt_path, 'file')
            load(alt_path);
        else
            warning('File not found for %s', file);
            continue;
        end
    end

    X_full = dff_rz';

    [t1, t2] = chop_by_file(file, valid_t);
    X_full   = X_full(t1:t2, :);

    behaviours = {zscore(MI_wheel_r(t1:t2)), zscore(dis_L2(t1:t2)), zscore(dis_R2(t1:t2)), ...
                  zscore(MI_whisker_r(t1:t2)), zscore(wa_r(t1:t2)), ...
                  zscore(wam_r(t1:t2)), zscore(wsp_r(t1:t2))};
    behaviour_names = {'Loco', 'LFL', 'RFL', 'WMI', 'WAngle', 'WAmp', 'WSP'};

    Y = [];
    for i = 1:length(behaviours)
        Y = [Y, behaviours{i}(:)];
    end

    ncomp                = numel(n_components_all);
    n_behaviours         = size(Y, 2);
    CVEV_all_reps        = zeros(ncomp, n_behaviours, n_repeats);
    CVEV_single_shuf_all = zeros(ncomp, n_behaviours, n_shuffles, n_behaviours, n_repeats);  % [k x beh x shuffle x target_b x rep]

    selected_neurons_all = {};

    for rep = 1:n_repeats
        fprintf('  > Subsample %d/%d ', rep, n_repeats);
        rng(2025 + file_i*1000 + rep);

        N_all  = size(X_full, 2);
        n_take = min(N_sub, N_all);
        if N_all < N_sub
            warning('%s: only %d neurons (<%d), using all of them', file, N_all, N_sub);
        end
        idx_keep                  = randperm(N_all, n_take);
        selected_neurons_all{rep} = idx_keep;
        X                         = X_full(:, idx_keep);

        % ----- 1) Real data -----
        CVEV_all_reps(:,:,rep) = pls_CVEV(X, Y, ncomp, K);

        % ----- 2) Shuffle control -----
        for target_b = 1:n_behaviours
            fprintf('.');

            for sh = 1:n_shuffles
                % Shuffle target behavior
                Y_shuf = Y;
                Y_shuf(:, target_b) = blockShuffle(Y(:, target_b), block_size);

                % Compute CVEV
                CVEV_shuf = pls_CVEV(X, Y_shuf, ncomp, K);

                % behaviorCVEV
                CVEV_single_shuf_all(:, target_b, sh, target_b, rep) = CVEV_shuf(:, target_b);
            end
        end
        fprintf(' done\n');
    end

    CVEV_all = mean(CVEV_all_reps, 3, 'omitnan');       % [k x beh]

    temp = mean(CVEV_single_shuf_all, 5, 'omitnan');    % [k x beh x sh x target_b]

    CVEV_shuf_diag_mean = zeros(ncomp, n_behaviours);
    CVEV_shuf_diag_std  = zeros(ncomp, n_behaviours);
    for bb = 1:n_behaviours
        % temp(:, bb, :, bb) -> [k x 1 x sh x 1]
        data_bb = squeeze(temp(:, bb, :, bb));          % [k x sh]
        CVEV_shuf_diag_mean(:, bb) = mean(data_bb, 2, 'omitnan');
        CVEV_shuf_diag_std(:, bb)  = std(data_bb, 0, 2, 'omitnan');
    end

    save_path = fullfile(results_dir, ['PLS4_' file '.mat']);
    save(save_path, 'file', 'n_components_all', 'behaviour_names', ...
        'CVEV_all', 'CVEV_single_shuf_all', ...
        'CVEV_shuf_diag_mean', 'CVEV_shuf_diag_std', ...
        'selected_neurons_all', 'N_sub', 'n_repeats', 'K', 'block_size', 'n_shuffles', ...
        '-v7.3');
    fprintf('Saved results to %s\n', save_path);

end

%% plot pls vs beh

results_dir = 'X:\MFB\MFB_AH_2023\PLSresults';
files       = dir(fullfile(results_dir, 'PLS4_*.mat'));

exclude_idx = [3, 4, 7, 9, 11];
exclude_idx = [];

all_CVEV_real = [];
all_CVEV_shuf = [];

valid_idx = setdiff(1:numel(files), exclude_idx);

for i = 1:numel(valid_idx)
    f    = valid_idx(i);
    data = load(fullfile(results_dir, files(f).name));

    if i == 1
        n_components_all = data.n_components_all;
        behaviour_names  = {'Loco', 'LFL', 'RFL', 'WMI', 'WAngle', 'WAmp', 'WSP'};
        n_behaviours     = length(behaviour_names);
    end

    all_CVEV_real(i,:,:) = data.CVEV_all;
    all_CVEV_shuf(i,:,:) = data.CVEV_shuf_diag_mean;
end

mean_real = squeeze(mean(all_CVEV_real, 1));            % [n_components x n_behaviours]
sem_real  = squeeze(std(all_CVEV_real, 0, 1) ./ sqrt(size(all_CVEV_real, 1)));

mean_shuf = squeeze(mean(all_CVEV_shuf, 1));
sem_shuf  = squeeze(std(all_CVEV_shuf, 0, 1) ./ sqrt(size(all_CVEV_shuf, 1)));

n_sessions   = size(all_CVEV_real, 1);
n_components = size(all_CVEV_real, 2);

beh_real_eachSession = squeeze(mean(all_CVEV_real, 3));                 % [n_sessions x n_components]
beh_shuf_eachSession = squeeze(mean(all_CVEV_shuf, 3));                 % [n_sessions x n_components]

beh_mean_real = mean(beh_real_eachSession, 1)';                         % [n_components x 1]
beh_sem_real  = (std(beh_real_eachSession, 0, 1) / sqrt(n_sessions))';  % [n_components x 1]

beh_mean_shuf = mean(beh_shuf_eachSession, 1)';                         % [n_components x 1]
beh_sem_shuf  = (std(beh_shuf_eachSession, 0, 1) / sqrt(n_sessions))';  % [n_components x 1]

mean_real = [mean_real, beh_mean_real];                                 % [n_components x (n_behaviours+1)]
sem_real  = [sem_real,  beh_sem_real];
mean_shuf = [mean_shuf, beh_mean_shuf];
sem_shuf  = [sem_shuf,  beh_sem_shuf];

behaviour_names{end+1} = 'All';
n_behaviours           = numel(behaviour_names);

colors = lines(n_behaviours);
figure('Position', [10 10 1550 750]);
for b = 1:n_behaviours
    subplot(2, ceil(n_behaviours/2), b);
    hold on;

    fill([n_components_all, fliplr(n_components_all)], ...
         [ (mean_real(:,b)'-sem_real(:,b)') , fliplr(mean_real(:,b)'+sem_real(:,b)') ], ...
         colors(b,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h1 = plot(n_components_all, mean_real(:,b), 'Color', colors(b,:), 'LineWidth', 2);

    fill([n_components_all, fliplr(n_components_all)], ...
         [ (mean_shuf(:,b)'-sem_shuf(:,b)') , fliplr(mean_shuf(:,b)'+sem_shuf(:,b)') ], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h2 = plot(n_components_all, mean_shuf(:,b), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

    title(behaviour_names{b}, 'Interpreter', 'none');
    xlabel('Number of PLS components');
    ylabel('CVEV');
    legend([h1, h2], {'PLS', 'Across block shuffle'}, 'Location', 'best', 'box', 'off');
    ylim([0, 1]);
    set(gca, 'FontSize', 15, ...
         'TickDir', 'out', ...
         'TickLength', [0.025 0.025], ...
         'LineWidth', 1.2);
    box off;
end

%%

% new_order = [1 2 3 5 6 7 8 4];
% behaviour_names = behaviour_names(new_order);
% mean_real = mean_real(:,new_order);
% sem_real  = sem_real(:,new_order);
% mean_shuf = mean_shuf(:,new_order);
% sem_shuf  = sem_shuf(:,new_order);
% all_CVEV_real = all_CVEV_real(:,:,new_order);
% all_CVEV_shuf = all_CVEV_shuf(:,:,new_order);

%%
trace_colors = [
    [0 153 76]/255
    [0 204 102]/255
    [0 204 204]/255
    [255 163 26]/255
    0.5 0.2 0.7
    0.1 0.2 0.7
    0.1 0.2 0.4
    0.3 0.3 0.3
];

n_sessions   = size(all_CVEV_real, 1);
n_components = size(all_CVEV_real, 2);

beh_real_eachSession = squeeze(mean(all_CVEV_real, 3, 'omitnan'));  % [S x C]
beh_shuf_eachSession = squeeze(mean(all_CVEV_shuf, 3, 'omitnan'));  % [S x C]

all_CVEV_real_aug = cat(3, all_CVEV_real, reshape(beh_real_eachSession, [n_sessions, n_components, 1]));
all_CVEV_shuf_aug = cat(3, all_CVEV_shuf, reshape(beh_shuf_eachSession, [n_sessions, n_components, 1]));

behaviour_names_aug = behaviour_names;
n_behaviours_aug    = numel(behaviour_names_aug);

mean_real_aug = squeeze(mean(all_CVEV_real_aug, 1, 'omitnan'));      % [C x B+1]
sem_real_aug  = squeeze(std(all_CVEV_real_aug, 0, 1, 'omitnan') ./ sqrt(n_sessions));

mean_shuf_aug = squeeze(mean(all_CVEV_shuf_aug, 1, 'omitnan'));
sem_shuf_aug  = squeeze(std(all_CVEV_shuf_aug, 0, 1, 'omitnan') ./ sqrt(n_sessions));

base_colors = trace_colors;
if size(base_colors, 1) < n_behaviours_aug
    extra        = lines(n_behaviours_aug - size(base_colors, 1));
    trace_colors = [base_colors; extra];
else
    trace_colors = base_colors;
end

nrows = 2;
ncols = ceil(n_behaviours_aug/nrows);
fig = figure('Color', 'w', 'Position', [10 10 850 400]);
tl  = tiledlayout(fig, nrows, ncols, 'TileSpacing', 'loose', 'Padding', 'compact');

for b = 1:n_behaviours_aug
    ax = nexttile(tl); hold(ax, 'on');

    base_col  = trace_colors(b,:);
    light_col = 0.65 + 0.35*base_col;

    for f = 1:n_sessions
        this_trace = squeeze(all_CVEV_real_aug(f,:,b));  % [1 x C]
        plot(n_components_all, this_trace, ...
             'Color', light_col, 'LineWidth', 0.75, ...
             'HandleVisibility', 'off');
    end

    fill([n_components_all, fliplr(n_components_all)], ...
         [mean_real_aug(:,b)'-sem_real_aug(:,b)', ...
          fliplr(mean_real_aug(:,b)'+sem_real_aug(:,b)')], ...
         base_col, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
    h1 = plot(n_components_all, mean_real_aug(:,b), 'Color', base_col, 'LineWidth', 2);

    % shuffle +/- SEM
    fill([n_components_all, fliplr(n_components_all)], ...
         [mean_shuf_aug(:,b)'-sem_shuf_aug(:,b)', ...
          fliplr(mean_shuf_aug(:,b)'+sem_shuf_aug(:,b)')], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h2 = plot(n_components_all, mean_shuf_aug(:,b), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

    title(ax, behaviour_names_aug{b}, 'Interpreter', 'none', 'FontWeight', 'normal');
    ylim(ax, [-0.05, 0.8]);
    yticks(ax, 0:0.4:0.8);
    xlim(ax, [1 8]);
    xticks(ax, [1 8]);
    box(ax, 'off');
    set(ax, 'FontSize', 18, ...
         'TickDir', 'out', ...
         'TickLength', [0.015 0.015], ...
         'LineWidth', 1.2);

    % --- drop repeated tick labels: keep only on the outer edge ---
    row = ceil(b/ncols);
    col = mod(b-1, ncols) + 1;
    if col ~= 1,     ax.YTickLabel = []; end
    if row ~= nrows, ax.XTickLabel = []; end

    % legend([h1, h2], {'PLS (mean \pm SEM)','Across block shuffle'}, 'Location','best');
end

% --- one shared label per axis for the whole grid ---
xlabel(tl, 'Number of PLS components', 'FontSize', 19.8);
ylabel(tl, 'CVEV',                     'FontSize', 19.8);


exportgraphics(fig, fullfile('C:\Users\Yizhou Xie\Downloads\beahviour  new analysis', 'PLS all.pdf'), 'ContentType', 'vector');

%% Optimal mode
threshold            = 0.01;
consecutive_required = 5;

n_files       = size(all_CVEV_real_aug, 1);
n_behaviours  = numel(behaviour_names);
best_mode_all = nan(n_behaviours, n_files);     % [behaviour x file]

for f = 1:n_files
    for b = 1:n_behaviours
        curve   = squeeze(all_CVEV_real_aug(f,:,b));
        diffs   = diff(curve);
        dim_opt = n_components_all(end);
        for i = 1:(numel(diffs) - consecutive_required + 1)
            if all(abs(diffs(i:i+consecutive_required-1)) < threshold)
                dim_opt = n_components_all(i);
                break;
            end
        end
        best_mode_all(b,f) = dim_opt;
    end
end

mean_best_mode = mean(best_mode_all, 2, 'omitnan');

% --- guard against silently clipped points ---
y_top = max([15, ceil(max(best_mode_all(:))/5)*5]);
if max(best_mode_all(:)) > 15
    warning('%d point(s) exceed 15 (max = %g); y-axis extended to %g', ...
        sum(best_mode_all(:) > 15), max(best_mode_all(:)), y_top);
end

% --- per-behaviour colours, matching the other panels ---
if exist('trace_colors', 'var') && size(trace_colors, 1) >= n_behaviours
    bar_cols = trace_colors(1:n_behaviours, :);
else
    bar_cols = lines(n_behaviours);
end
bar_cols(end,:) = [0.45 0.45 0.45];             % 'All' = derived summary

% --- x positions: small gap before the last ('All') category ---
xpos        = 1:n_behaviours;

figure('Position', [10 10 560 420]); hold on;

for b = 1:n_behaviours
    bar(xpos(b), mean_best_mode(b), 0.62, ...
        'FaceColor', 0.45*bar_cols(b,:) + 0.55, ...   % lightened fill
        'EdgeColor', 'none');
end

dot_step = 0.1;                                % fixed spacing -> density is readable
for b = 1:n_behaviours
    y = best_mode_all(b,:)';  y = y(~isnan(y));
    if isempty(y), continue; end
    x = xpos(b) * ones(size(y));
    [vals,~,grp] = unique(y, 'stable');
    for g = 1:numel(vals)
        idx    = find(grp == g);
        m      = numel(idx);
        offs   = ((1:m) - (m+1)/2) * dot_step;   % centred, constant spacing
        x(idx) = x(idx) + offs(:);
    end
    scatter(x, y, 26, 'filled', ...
        'MarkerFaceColor', bar_cols(b,:), 'MarkerFaceAlpha', 0.85, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
end

xticks(xpos);
xticklabels(behaviour_names);
xlim([xpos(1)-0.7, xpos(end)+0.7]);
ylabel('Optimal PLS mode');
ylim([0 y_top]);
yticks(0:5:y_top);
set(gca, 'FontSize', 18, ...
     'TickDir', 'out', ...
     'TickLength', [0.015 0.015], ...
     'LineWidth', 1.2);
box off;
exportgraphics(gcf, fullfile('C:\Users\Yizhou Xie\Downloads\beahviour  new analysis', ...
    'PLS opt modes.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');

%% bar +scatter

n_files      = size(all_CVEV_real_aug, 1);
n_behaviours = numel(behaviour_names_aug);

if ~exist('best_mode_all', 'var') || isempty(best_mode_all) || size(best_mode_all, 1) ~= n_behaviours
    threshold = 0.01; consecutive_required = 5;
    best_mode_all = nan(n_behaviours, n_files);
    for b = 1:n_behaviours
        for f = 1:n_files
            curve   = squeeze(all_CVEV_real_aug(f,:,b));
            diffs   = diff(curve);
            dim_opt = n_components_all(end);
            for i = 1:(numel(diffs) - consecutive_required + 1)
                if all(abs(diffs(i:i+consecutive_required-1)) < threshold)
                    dim_opt = n_components_all(i); break;
                end
            end
            best_mode_all(b,f) = dim_opt;
        end
    end
end

col_first   = [0.7 0.7 1.0];
col_opt     = [1.0 0.7 0.4];
col_shuffle = [0.6 0.6 0.6];

first_vals  = nan(n_behaviours, n_files);
opt_vals    = nan(n_behaviours, n_files);
shuf_vals   = nan(n_behaviours, n_files);

for b = 1:n_behaviours
    idx_first = 1;
    for f = 1:n_files
        opt_mode = best_mode_all(b,f);
        idx_opt  = find(n_components_all == opt_mode, 1, 'first');
        if isempty(idx_opt), idx_opt = numel(n_components_all); end

        first_vals(b,f) = all_CVEV_real_aug(f, idx_first, b);
        opt_vals(b,f)   = all_CVEV_real_aug(f, idx_opt,   b);
        shuf_vals(b,f)  = all_CVEV_shuf_aug(f, idx_opt,   b);
    end
end

mean_first = mean(first_vals, 2, 'omitnan');
mean_opt   = mean(opt_vals,   2, 'omitnan');
mean_shuf  = mean(shuf_vals,  2, 'omitnan');

sem_first = std(first_vals, 0, 2, 'omitnan')/sqrt(n_files);
sem_opt   = std(opt_vals,   0, 2, 'omitnan')/sqrt(n_files);
sem_shuf  = std(shuf_vals,  0, 2, 'omitnan')/sqrt(n_files);

groupMeans = [mean_first, mean_opt, mean_shuf];
groupSEMs  = [sem_first,  sem_opt,  sem_shuf];

figure('Position', [80 80 680 400]); hold on;

bh = bar(groupMeans, 'grouped', 'EdgeColor', 'none');
bh(1).FaceColor = col_first;
bh(2).FaceColor = col_opt;
bh(3).FaceColor = col_shuffle;

if isprop(bh(1), 'XEndPoints')
    x_first = bh(1).XEndPoints(:);
    x_opt   = bh(2).XEndPoints(:);
    x_shuf  = bh(3).XEndPoints(:);
else
    x = 1:n_behaviours; w = 0.25;
    x_first = x - w; x_opt = x; x_shuf = x + w;
end

rng(0); jitter = 0.0;
for b = 1:n_behaviours
    for f = 1:n_files
        y1  = first_vals(b,f);
        y2  = opt_vals(b,f);
        y3  = shuf_vals(b,f);
        dx1 = (rand-0.5)*2*jitter; dx2 = (rand-0.5)*2*jitter; dx3 = (rand-0.5)*2*jitter;

        plot([x_first(b)+dx1, x_opt(b)+dx2, x_shuf(b)+dx3], [y1, y2, y3], ...
             '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.7, 'HandleVisibility', 'off');

        scatter(x_first(b)+dx1, y1, 13, 'filled', ...
                'MarkerFaceColor', col_first*0.9,   'MarkerEdgeColor', [0.3,0.3,0.3], 'LineWidth', 0.1);
        scatter(x_opt(b)+dx2,   y2, 13, 'filled', ...
                'MarkerFaceColor', col_opt*0.9,     'MarkerEdgeColor', [0.3,0.3,0.3], 'LineWidth', 0.1);
        scatter(x_shuf(b)+dx3,  y3, 13, 'filled', ...
                'MarkerFaceColor', col_shuffle*0.9, 'MarkerEdgeColor', [0.3,0.3,0.3], 'LineWidth', 0.1);
    end
end

errorbar(x_first, groupMeans(:,1), groupSEMs(:,1), 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 5);
errorbar(x_opt,   groupMeans(:,2), groupSEMs(:,2), 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 5);
errorbar(x_shuf,  groupMeans(:,3), groupSEMs(:,3), 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 5);

xticks(1:n_behaviours);
xticklabels(behaviour_names_aug);
ylabel('CVEV');
ylim([-0.2, 1]);
set(gca, 'FontSize', 18, ...
     'TickDir', 'out', ...
     'TickLength', [0.015 0.015], ...
     'LineWidth', 1.2);
box off;
legend({'First mode', 'Optimal mode', 'Shuffle'}, 'Location', 'best', 'Box', 'off');
yticks([0 0.5 1])
exportgraphics(gcf, fullfile('C:\Users\Yizhou Xie\Downloads\beahviour  new analysis', 'PLS summary.pdf'), 'ContentType', 'vector');
