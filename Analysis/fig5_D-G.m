
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_13_56_55'; % remove loco, state (only small portion at the begin)
   '191018_14_30_00';
   '191018_14_11_33'; % remove loco, state (only small portion at the end)
   '191209_13_44_12'; 
   '191209_14_04_14'; % remove loco, state (only small portion at the end)
   '191209_14_32_39'; % remove loco
   '191209_15_01_22'; % remove loco
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };


num_files = length(folder_names);
for file_i = 1% 2:num_files
    file = char(folder_names(file_i));
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
    valid_t= 1:length(L_state);


    dis = cumsum(MI_wheel_r) * (1/100);

    acc_W = [0 diff(MI_wheel_r)] / (1/100);
    acc_R = [0 diff(dis_R2)] / (1/100);
    acc_L = [0 diff(dis_L2)] / (1/100);
    acc_wsk = [0 diff(MI_whisker_r)] / (1/100);

    win = 10; % smoothing window, adjust as needed

acc_W = smoothdata([0, diff(MI_wheel_r)] * 100, 'gaussian', win);
acc_R = smoothdata([0, diff(dis_R2)] * 100, 'gaussian', win);
acc_L = smoothdata([0, diff(dis_L2)] * 100, 'gaussian', win);


    behaviours = {MI_wheel_r,dis_L2,dis_R2,...
        L_state,MI_whisker_r,wa_r,wam_r,wsp_r,acc_W, acc_R2, acc_L2,dis,acc_wsk};
    behaviour_names = {'Loco','LFL','RFL' ,'STATE','WMI','WAngle','WAmp','WSP','accW','accR','accL','dis','accwsk'};
    results2 = analyze_behaviours_MFnum(dff_rz', L_state, valid_t, behaviours, behaviour_names, 'sorted',file);
    %save(['Ridgebeh_' file '.mat'], 'results');
end

%% plot fig 6a

save_dir = 'X:\MFB\MFB_AH_2023\update fig 5_6';
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

trace_colors = [
    [0 153 76]/255
    [0 204 102]/255
    [0 204 204]/255
    0.8 0.1 0.6
    0.5 0.2 0.7
    0.1 0.2 0.7
    0.1 0.2 0.4
    1.0 0.5 0.1
];

t1 = 10000;
t2 = 40400;
behaviour_names = {'Loco','LFL','RFL','WMI','WAngle','WAmp','WSP'};

T_min = inf;
for b = 1:numel(behaviour_names)
    T_min = min(T_min, size(results.(behaviour_names{b}).Ytest_trace, 2));
end
t2_glob = min(t2, T_min);
t1_glob = max(1, min(t1, t2_glob));
win_len = t2_glob - t1_glob + 1;
pad_left     = max(5,  round(0.02 * win_len));
pad_right    = max(20, round(0.10 * win_len));
x_text_right = pad_left + win_len + round(0.02 * win_len);

n_beh = numel(behaviour_names);
nrows = 2; ncols = ceil(n_beh/2);
fig = figure('Name','All Behaviour Traces','Color','w','Position',[100 0 1250 850]);
tl = tiledlayout(fig, nrows, ncols, 'TileSpacing','compact', 'Padding','compact');

spacing = 1.2;

ax_all  = gobjects(1, n_beh);
txt_all = gobjects(0);

for b = 1:n_beh
    name = behaviour_names{b};
    res  = results.(name);

    Ypred = res.Ypred_trace;      % [num_traces x T]
    Ytest = res.Ytest_trace;      % [num_traces x T]
    trace_neurons = res.trace_neurons;
    neuron_range  = res.neuron_range;
    mean_CVEV     = res.mean_CVEV;

    num_traces = numel(trace_neurons);
    T_beh = size(Ytest, 2);

    idx_end = min(T_beh, t1_glob + win_len - 1);
    idx_t = t1_glob:idx_end;
    win_len_b = numel(idx_t);
    X = pad_left + (1:win_len_b);

    ax = nexttile(tl, b); hold(ax, 'on');

    for i = 1:num_traces
        offset = (num_traces - i) * spacing;

        if win_len_b > 0
            yt = Ytest(i, idx_t);
            yp = Ypred(i, idx_t);
           

            if b~=8
            %yt= movmean(yt,5);
            end
    
           % yp = movmean(yp,5);

            m = min(yt); M = max(yt);
            if isfinite(M) && (M > m)
        m = min(yt); M = max(yt);
        yt_s = (yt - m) / (M - m);
        yp_s = (yp - m) / (M - m);
            else
                yt_s = yt*0 + 0.5;
                yp_s = yp*0 + 0.5;
            end

            plot(ax, X, yt_s + offset, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            plot(ax, X, yp_s + offset, 'Color', trace_colors(b,:), 'LineWidth', 1);
        end

        k = trace_neurons(i);
        cvev_idx = find(neuron_range == k, 1, 'first');
        if ~isempty(cvev_idx)
            cvev_val = mean_CVEV(cvev_idx);
            ht = text(ax, x_text_right, offset, sprintf('%.2f', cvev_val), ...
                 'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
                 'FontSize', 18, 'Clipping','on');
            txt_all(end+1) = ht;
        end
    end

    ylim(ax, [-0.2*spacing, (num_traces-1)*spacing + 1.2*spacing]);

    axis(ax, 'off'); box(ax, 'off');
    title(ax, name, 'Interpreter','none', 'FontSize', 18, 'FontWeight','normal');
    set(ax, 'PositionConstraint','innerposition');
    ax_all(b) = ax;
end

x_lo = 1;
x_hi = pad_left + win_len + pad_right;
set(ax_all, 'XLim', [x_lo, x_hi]);
drawnow;

w_frac = 0;
for h = txt_all
    e = get(h, 'Extent');
    w_frac = max(w_frac, e(3));
end
w_frac = w_frac / (x_hi - x_lo);

gap_frac = 0.02;
x_hi_new = x_lo + (x_text_right - x_lo) / max(0.05, 1 - w_frac - gap_frac);

set(ax_all, 'XLim', [x_lo, x_hi_new]);

 exportgraphics(fig, fullfile(save_dir, 'Behaviour_traces_all_righttext.pdf'), 'ContentType','vector');
% exportgraphics(fig, fullfile(save_dir, 'Behaviour_traces_all_righttext.png'), 'Resolution', 300);


%% fig 6b and c
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_14_11_33'; % remove loco, state (only small portion at the end)
   '191018_14_30_00';
   '191018_13_56_55'; % remove loco, state (only small portion at the begin)
   '191209_13_44_12'; 
   '191209_14_04_14'; % remove loco, state (only small portion at the end)
   '191209_14_18_13';
   '191209_14_32_39'; % remove loco
   '191209_14_46_58';
   '191209_15_01_22'; % remove loco
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };

savepath = 'X:\MFB';

behavior_fields = {'Loco','LFL','RFL' ,'STATE','WMI','WAngle','WAmp','WSP'};
%behavior_fields = {'Loco','LFL','RFL' ,'STATE','WMI','WAngle','WAmp','WSP','accW','accR','accL','dis'};
num_files = length(folder_names);
results_all = struct();

for file_i = 1:num_files
    file = char(folder_names(file_i));
    main_path = fullfile('X:\MFB\MFB_AH_2023\Data', file, ['Ridgebeh_' file '.mat']);
    
    if exist(main_path, 'file')
        load(main_path);
    else
        alt_path = fullfile('X:\MFB\MFB_AH_2023\Data', file, 'Functional_data', ['Ridgebeh_' file '.mat']);
        if exist(alt_path, 'file')
            load(alt_path);
        else
            warning('File not found for %s', file);
            continue;
        end
    end

    for b = 1:length(behavior_fields)
        field = behavior_fields{b};
        if strcmp(field,'Loco') && strcmp(file, '191209_15_01_22')
            continue
        end
        if strcmp(field,'Spd_r') && strcmp(file, '191209_15_01_22')
            continue

        end
        
%         if ~isfield(results_all.(field), 'CVEV')
%             results_all.(field).CVEV = nan(num_files, length(results.(field).mean_CVEV));
%         end
            
    results_all.(field).CVEV(file_i, :) = results.(field).mean_CVEV(:)';
    results_all.(field).neurons = results.(field).neuron_range(:)';
    
    results_all.(field).BestMFidx(file_i, :) = results.(field).BestMFidx(1:100);
    results_all.(field).BestMFnum(file_i,1) = results.(field).OptMFnum;
    results_all.(field).BestMFCVEV(file_i,1) = results.(field).OptMFCVEV;
    results_all.(field).EachMFCVEV(file_i, :) = results.(field).EachMFCVEV(1:100);
    
    results_all.(field).OptMF_ShuffleAcross(file_i,1) = results.(field).OptMF_ShuffleAcross;
%    results_all.(field).OptMF_ShuffleWithin(file_i,1) = results.(field).OptMF_ShuffleWithin;
    
    results_all.(field).pm(file_i, :)         = results.(field).property.pm.CVEV(:)';
    results_all.(field).pm_shuf(file_i, :)    = results.(field).property.pm_shuf.CVEV(:)';
    results_all.(field).nm(file_i, :)         = results.(field).property.nm.CVEV(:)';
    results_all.(field).nm_shuf(file_i, :)    = results.(field).property.nm_shuf.CVEV(:)';
    results_all.(field).nsm(file_i, :)        = results.(field).property.nsm.CVEV(:)';
    results_all.(field).nsm_shuf(file_i, :)   = results.(field).property.nsm_shuf.CVEV(:)';
    results_all.(field).random(file_i, :)     = results.(field).property.random.CVEV(:)';
    results_all.(field).random_shuf(file_i,:) = results.(field).property.random_shuf.CVEV(:)';
 
    end
end


n_beh = numel(behavior_fields)-1;
plot_order = [1 2 3 5 6 7 8];
plot_fields = behavior_fields(plot_order);

neuron_range = results_all.(plot_fields{1}).neurons;
n_neurons = length(neuron_range);

num_files = size(results_all.(plot_fields{1}).CVEV, 1);

all_behaviors_CVEV_curve = nan(num_files, n_neurons, n_beh);

for b = 1:n_beh
    field = plot_fields{b};
    all_behaviors_CVEV_curve(:, :, b) = results_all.(field).CVEV;
end

all_mean_per_session = mean(all_behaviors_CVEV_curve, 3, 'omitnan');  % [files x neurons]

all_behaviors_CVEV_curve_aug = cat(3, all_behaviors_CVEV_curve, all_mean_per_session);

plot_fields_aug = [plot_fields, {'All'}];
n_beh_aug = n_beh+1;

nrows = 2;
ncols = ceil(n_beh_aug/2);
figure('Position',[80 80 1150 450]);

for b = 1:n_beh_aug
    field = plot_fields_aug{b};
    subplot(nrows, ncols, b); hold on;
    
    CVEV_mat = all_behaviors_CVEV_curve_aug(:, :, b);  % [files x neurons]
    
    if strcmp(field, 'All')
        valid_idx = setdiff(1:num_files, [3 5 7 9 11]);
    elseif strcmp(field,'Loco') || strcmp(field,'LFL') || strcmp(field,'RFL') || strcmp(field,'STATE')
        valid_idx = setdiff(1:num_files, [3 5 7 9 11]);
    else
        valid_idx = 1:num_files;
    end
    
    for file_i = valid_idx
        plot(neuron_range, CVEV_mat(file_i, :), ...
             'Color', [0 0.6 0.6], 'LineWidth', 1);
    end
    
    CVEV_mean = mean(CVEV_mat(valid_idx, :), 1, 'omitnan');
    CVEV_sem  = std(CVEV_mat(valid_idx, :), 0, 1, 'omitnan') / sqrt(numel(valid_idx));
    
    fill([neuron_range, fliplr(neuron_range)], ...
         [CVEV_mean + CVEV_sem, fliplr(CVEV_mean - CVEV_sem)], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    plot(neuron_range, CVEV_mean, 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
    
    row = ceil(b/ncols);
    col = mod(b-1, ncols) + 1;
    
    if col == 1
        ylabel('CVEV', 'FontSize', 18);
    end
    
    if row == nrows
        xlabel('MFA #', 'FontSize', 18);
    end
    
    ylim([0 1]); 
    yticks([0 0.5 1]);
    xlim([1 120]); 
    xticks([1 60 120]);
    box off;
    title(field, 'FontWeight', 'normal');
    set(gca, 'FontSize', 18, 'LineWidth', 1);
end




 exportgraphics(gcf, fullfile(save_dir, 'CVEV_all_behaviours_subplot.pdf'), 'ContentType','vector');

%% 6c
% ===== 3) plot (including 'All') =====

nrows = 2;

ncols = ceil(n_beh_aug/2);

figure('Position', [80 80 700 400]);

t = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');



for b = 1:n_beh_aug

    ax = nexttile; hold on;

    field      = behavior_fields_aug{b};

    top10_CVEV = all_behaviors_CVEV_aug(:,:,b);   % [files x top_n]



    if strcmp(field, 'All')

        valid_idx = 1:size(top10_CVEV, 1);

    elseif ismember(field, {'Loco','LFL','RFL','STATE'})

        valid_idx = setdiff(1:size(top10_CVEV,1), [3 5 7 9 11]);

    else

        valid_idx = 1:size(top10_CVEV, 1);

    end



    valid_data = top10_CVEV(valid_idx, :);

    mean_CVEV  = mean(valid_data, 1, 'omitnan');

    sem_CVEV   = std(valid_data, 0, 1, 'omitnan') / sqrt(numel(valid_idx));



    bar(1:top_n, mean_CVEV, 'FaceColor', [0 0.6 0.6], 'EdgeColor', 'none');

    errorbar(1:top_n, mean_CVEV, sem_CVEV, 'k.', 'LineWidth', 1.5);



    ylim([0 0.5]);   yticks([0 0.25 0.5]);

    xlim([0.5 top_n+0.5]);   xticks([1 5 10]);

    title(field, 'FontWeight', 'normal');

    box off;

    set(ax, 'FontSize', 18, 'LineWidth', 1, 'TickDir', 'out');



    % --- drop repeated tick labels: keep only on the outer edge ---

    row = ceil(b/ncols);

    col = mod(b-1, ncols) + 1;

    if col ~= 1,     ax.YTickLabel = []; end

    if row ~= nrows, ax.XTickLabel = []; end

end



% --- one shared label for the whole grid ---

xlabel(t, 'Number of most informative MFAs', 'FontSize', 18);

ylabel(t, 'CVEV',      'FontSize', 18);



exportgraphics(gcf, fullfile(save_dir, 'CVEV_indiv_behaviours_subplot.pdf'), ...
    'ContentType', 'vector');

%% select behaviour for summary

selected_behaviours = {'Loco','LFL','RFL','WMI','WAngle','WAmp','WSP'};

%selected_behaviours = {'Loco','LFL','RFL' ,'STATE','WMI','WAngle','WAmp','WSP','accW','accR','accL','dis'};



plot_order = [1 2 3 5 6 7 8];
plot_fields = behavior_fields(plot_order);
num_behaviours = length(plot_fields);

CVEV_first_all = cell(num_behaviours, 1);
CVEV_best_all = cell(num_behaviours, 1);
CVEV_within_all = cell(num_behaviours, 1);
CVEV_across_all = cell(num_behaviours, 1);
optnum_all = cell(num_behaviours, 1);

valid_idx = [1,2,4,6,8,10,12,13,14,15,16];
n_valid = length(valid_idx);

all_first_mat = nan(n_valid, num_behaviours);
all_best_mat = nan(n_valid, num_behaviours);
all_within_mat = nan(n_valid, num_behaviours);
all_across_mat = nan(n_valid, num_behaviours);

for b = 1:num_behaviours
    field = plot_fields{b};

    if ~isfield(results_all, field)
        warning(['Field [', field, '] not found, skipping']);
        continue;
    end

    CVEV_mat = results_all.(field).CVEV;
    CVEV_best = results_all.(field).BestMFCVEV;
    CVEV_across = results_all.(field).OptMF_ShuffleAcross;
%    CVEV_within = results_all.(field).OptMF_ShuffleWithin;
    CVEV_first = CVEV_mat(:, 1);

    CVEV_first = CVEV_first(valid_idx);
    CVEV_best = CVEV_best(valid_idx);
%    CVEV_within = CVEV_within(valid_idx);
    CVEV_across = CVEV_across(valid_idx);

    CVEV_first_all{b} = CVEV_first;
    CVEV_best_all{b} = CVEV_best;
 %   CVEV_within_all{b} = CVEV_within;
    CVEV_across_all{b} = CVEV_across;

    all_first_mat(:, b) = CVEV_first;
    all_best_mat(:, b) = CVEV_best;
%    all_within_mat(:, b) = CVEV_within;
    all_across_mat(:, b) = CVEV_across;

    optnum = results_all.(field).BestMFnum;
    optnum_all{b} = optnum(valid_idx,:);

    if length(CVEV_first) > 5
        p1 = signrank(CVEV_first, CVEV_best);
%        p2 = signrank(CVEV_best, CVEV_within);
        p3 = signrank(CVEV_best, CVEV_across);
        disp(['[', field, '] first vs best p = ', num2str(p1), ...
              ', best vs across p = ', num2str(p3)]);
    else
        disp(['[', field, '] too few samples, skipping stats']);
    end
end

all_first_mean = mean(all_first_mat, 2, 'omitnan');    % [n_valid x 1]
all_best_mean = mean(all_best_mat, 2, 'omitnan');
all_within_mean = mean(all_within_mat, 2, 'omitnan');
all_across_mean = mean(all_across_mat, 2, 'omitnan');

CVEV_first_all{end+1} = all_first_mean;
CVEV_best_all{end+1} = all_best_mean;
CVEV_within_all{end+1} = all_within_mean;
CVEV_across_all{end+1} = all_across_mean;

plot_fields_aug = [plot_fields, {'All'}];
num_behaviours_aug = num_behaviours + 1;

xtick_labels = plot_fields_aug;

figure('Position', [100 100 700 400]); hold on;

col_first   = [0.7 0.7 1.0];
col_opt     = [1.0 0.7 0.4];
col_within  = [0.80 0.80 0.80];
col_across  = [0.5 0.5 0.5];
col_link    = [0.70 0.70 0.70];

num_b = num_behaviours_aug;
mean_first  = nan(num_b,1); sem_first  = nan(num_b,1);
mean_best   = nan(num_b,1); sem_best   = nan(num_b,1);
mean_within = nan(num_b,1); sem_within = nan(num_b,1);
mean_across = nan(num_b,1); sem_across = nan(num_b,1);

for b = 1:num_b
    y1 = CVEV_first_all{b};
    y2 = CVEV_best_all{b};
    y3 = CVEV_within_all{b};
    y3 = CVEV_across_all{b};

    y1n = y1(~isnan(y1)); n1 = numel(y1n);
    y2n = y2(~isnan(y2)); n2 = numel(y2n);
    y3n = y3(~isnan(y3)); n3 = numel(y3n);
    y4n = y3(~isnan(y3)); n4 = numel(y4n);

    mean_first(b)  = mean(y1n, 'omitnan');
    mean_best(b)   = mean(y2n, 'omitnan');
    mean_within(b) = mean(y3n, 'omitnan');
    mean_across(b) = mean(y4n, 'omitnan');

    sem_first(b)  = std(y1n, 0, 'omitnan')/max(1,sqrt(n1));
    sem_best(b)   = std(y2n, 0, 'omitnan')/max(1,sqrt(n2));
    sem_within(b) = std(y3n, 0, 'omitnan')/max(1,sqrt(n3));
    sem_across(b) = std(y4n, 0, 'omitnan')/max(1,sqrt(n4));
end

bar_vals = [mean_first, mean_best, mean_across];

bar_handle = bar(bar_vals, 'grouped', 'BarWidth', 0.7, 'EdgeColor','none');
bar_handle(1).FaceColor = col_first;
bar_handle(2).FaceColor = col_opt;
bar_handle(3).FaceColor = col_across;

if isprop(bar_handle(1), 'XEndPoints')
    x_first  = bar_handle(1).XEndPoints(:);
    x_best   = bar_handle(2).XEndPoints(:);
    x_across = bar_handle(3).XEndPoints(:);

else
    x = (1:num_b)'; w = 0.22;
    x_first  = x + (-1.5)*w;
    x_best   = x + (-0.5)*w;
    x_across = x + ( 0.5)*w;
end

for b = 1:num_b
    y1 = CVEV_first_all{b};
    y2 = CVEV_best_all{b};
    y3 = CVEV_across_all{b};
    n  = max([numel(y1), numel(y2), numel(y3), numel(y3)]);
    
    for i = 1:n
        v1 = i <= numel(y1) && ~isnan(y1(i));
        v2 = i <= numel(y2) && ~isnan(y2(i));
        v3 = i <= numel(y3) && ~isnan(y3(i));
        
        if v1 && v2 && v3
            plot([x_first(b), x_best(b), x_across(b)], ...
                 [y1(i),      y2(i),     y3(i)], ...
                 '--', 'Color', col_link, 'LineWidth', 0.7, 'HandleVisibility','off');
        end
        
        if v1, scatter(x_first(b),  y1(i), 16, 'filled', ...
                       'MarkerFaceColor', col_first*0.9,  'MarkerEdgeColor',[0.3 0.3 0.3], 'LineWidth',0.1, 'HandleVisibility','off'); end
        if v2, scatter(x_best(b),   y2(i), 16, 'filled', ...
                       'MarkerFaceColor', col_opt*0.9,    'MarkerEdgeColor',[0.3 0.3 0.3], 'LineWidth',0.1, 'HandleVisibility','off'); end
        if v3, scatter(x_across(b), y3(i), 16, 'filled', ...
                       'MarkerFaceColor', col_across*0.9, 'MarkerEdgeColor',[0.7 0.7 0.7], 'LineWidth',0.1, 'HandleVisibility','off'); end
    end
end

eh(1) = errorbar(x_first,  mean_first,  sem_first,  'k', 'LineStyle','none', 'LineWidth',1.2, 'CapSize',5, 'HandleVisibility','off');
eh(2) = errorbar(x_best,   mean_best,   sem_best,   'k', 'LineStyle','none', 'LineWidth',1.2, 'CapSize',5, 'HandleVisibility','off');
eh(3) = errorbar(x_across, mean_across, sem_across, 'k', 'LineStyle','none', 'LineWidth',1.2, 'CapSize',5, 'HandleVisibility','off');
for i = 1:numel(eh), uistack(eh(i), 'top'); end

xlim([0.5, num_b + 0.5]);
ylim([-0.3 1]);
yticks([0 0.5 1]);
xticks(1:num_b);
xticklabels(xtick_labels);
ylabel('CVEV');
legend({'First MF', 'All MFs', 'Shuffle'}, 'Location', 'best', 'Box', 'off');
set(gca, 'FontSize', 18, 'LineWidth', 1);
box off;
filename = ['SelectedBehaviours_FirstVsBestMF_Shuffle_', '.pdf'];
exportgraphics(gcf, fullfile(save_dir, filename), 'ContentType', 'vector');


%% opt num

num_groups = length(optnum_all);

opt_mean = zeros(num_groups, 1);
opt_sem  = zeros(num_groups, 1);

for i = 1:num_groups
    data_i = optnum_all{i};
    
    opt_mean(i) = mean(data_i, 'omitnan');
    
    n_valid = sum(~isnan(data_i));
    opt_sem(i)  = std(data_i, 0, 'omitnan') / sqrt(n_valid);
end

figure('Color', 'w', 'Position', [100, 100, 600, 400]); 
hold on;

bh = bar(opt_mean, 'FaceColor', [0.77 0.77 0.77], ...
         'EdgeColor', 'none', 'BarWidth', 0.6);

er = errorbar(1:num_groups, opt_mean, opt_sem, ...
              'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

jitter_width = 0.1; 
point_color = [0.3 0.3 0.3];
point_size = 25;

for i = 1:num_groups
    y_data = optnum_all{i};
    num_points = length(y_data);
    
    x_jitter = i + (rand(num_points, 1) - 0.5) * jitter_width;
    
    scatter(x_jitter, y_data, point_size, ...
            'MarkerFaceColor', point_color, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.5);
end
xticks(1:num_behaviours)

xticklabels(xtick_labels)

xlim([0.5, num_groups + 0.5])
xticks(1:num_groups)
xticklabels(xtick_labels) 
xtickangle(45)

ylabel('Optimal MF #')
set(gca, 'FontSize', 18, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off')

hold off;

filename_mfnum = ['OptimalMFnum_byBehaviour', '.pdf'];
exportgraphics(gcf, fullfile(save_dir, filename_mfnum), 'ContentType', 'vector');


%%
% ===== 3) update behaviour list and plot =====
behavior_fields2_aug = [behavior_fields2, {'All'}];
n_beh_aug = n_beh + 1;

nrows = 2;
ncols = ceil(n_beh_aug/nrows);
fig = figure('Color','w','Position',[80 80 700 420]);
tl = tiledlayout(fig, nrows, ncols, 'TileSpacing','compact','Padding','compact');

for b = 1:n_beh_aug
    field = behavior_fields2_aug{b};
    ax = nexttile(tl); hold(ax,'on');

    for group = {'pm','nm','nsm','random'}
        g = group{1};
        mat      = results_all.(field).(g);
        mat_shuf = results_all.(field).([g '_shuf']);

        mat      = mat(valid_idx,:);
        mat_shuf = mat_shuf(valid_idx,:);

        m = mean(mat, 1, 'omitnan');
        s = std(mat, 0, 1, 'omitnan') / sqrt(size(mat, 1));

        fill(ax, [n_range, fliplr(n_range)], [m + s, fliplr(m - s)], ...
            colors.([g '_shuf']), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(ax, n_range, m, '-', 'Color', colors.(g), 'LineWidth', 2);

        m_shuf = mean(mat_shuf, 1, 'omitnan');
        plot(ax, n_range, m_shuf, '--', 'Color', colors.(g), 'LineWidth', 1.2);
    end

    ylim(ax, [-0.1 0.6]);   yticks(ax, [0 0.3 0.6]);
    xlim(ax, [1 30]);       xticks(ax, [1 15 30]);
    box(ax, 'off');
    title(ax, field, 'FontWeight', 'normal');
    set(ax, 'FontSize', 18, 'LineWidth', 1, 'TickDir', 'out');

    % --- drop repeated tick labels: keep only on the outer edge ---
    row = ceil(b/ncols);
    col = mod(b-1, ncols) + 1;
    if col ~= 1,     ax.YTickLabel = []; end
    if row ~= nrows, ax.XTickLabel = []; end
end

% --- one shared label per axis for the whole grid ---
xlabel(tl, 'Number of MFAs', 'FontSize', 18);
ylabel(tl, 'CVEV',           'FontSize', 18);

% --- legend: colour encodes group, line style encodes real vs shuffle ---
ax_leg = nexttile(tl, 1);                       % any axes works as the handle parent
h_leg  = gobjects(1, numel(groups) + 1);
lbl    = cell(1, numel(groups) + 1);
for g_idx = 1:numel(groups)
    g = groups{g_idx};
    h_leg(g_idx) = plot(ax_leg, nan, nan, '-', 'Color', colors.(g), 'LineWidth', 2);
    lbl{g_idx}   = legend_names.(g);
end
h_leg(end) = plot(ax_leg, nan, nan, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);
lbl{end}   = 'Shuffle';

lg = legend(h_leg, lbl, 'Box','off', 'FontSize', 16, ...
    'NumColumns', numel(lbl), 'Interpreter', 'tex');
lg.Layout.Tile = 'south';

exportgraphics(fig, fullfile(save_dir, 'all_behaviour_property_cvev.pdf'), 'ContentType','vector');

%% ===== NM vs shuffle: per-behaviour x per-MF-count paired sign-rank =====
p_results = struct();

for b = 1:numel(behavior_fields2)
    field = behavior_fields2{b};

    nm      = results_all.(field).nm;        % [n_sessions × numel(n_range)]
    nm_shuf = results_all.(field).nm_shuf;   % same size

    K = numel(n_range);
    pvec  = nan(1, K);
    zvec  = nan(1, K);
    npair = zeros(1, K);
    meddf = nan(1, K);

    for k = 1:K
        x = nm(:,k);
        y = nm_shuf(:,k);

        valid = ~isnan(x) & ~isnan(y);
        x = x(valid);  y = y(valid);
        npair(k) = numel(x);

        if npair(k) > 0
            [p,~,stats] = signrank(x, y);
            pvec(k) = p;
            if isfield(stats,'zval'), zvec(k) = stats.zval; end
            meddf(k) = median(x - y);
        end
    end

    qvec = nan(size(pvec));
    [ps, idx] = sort(pvec, 'ascend', 'MissingPlacement','last');
    m = sum(~isnan(ps));
    if m > 0
        rk = 1:m;
        q  = zeros(1,m);
        q(end) = ps(end) * m / rk(end);
        for i = m-1:-1:1
            q(i) = min(ps(i) * m / rk(i), q(i+1));
        end
        qvec(idx(1:m)) = q;
    end

    p_results.(field).p            = pvec;
    p_results.(field).q            = qvec;
    p_results.(field).z            = zvec;
    p_results.(field).median_diff  = meddf;       % median(x - y)
    p_results.(field).n_pairs      = npair;
    p_results.(field).n_range      = n_range;
end

%%
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_13_56_55'; % remove loco, state (only small portion at the begin)
   '191018_14_30_00';
   '191018_14_11_33'; % remove loco, state (only small portion at the end)
   '191209_13_44_12'; 
   '191209_14_04_14'; % remove loco, state (only small portion at the end)
   '191209_14_32_39'; % remove loco
   '191209_15_01_22'; % remove loco
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };


num_files = length(folder_names);
all_wsp = cell(num_files, 1);

all_Loco= cell(num_files, 1);
all_Wsk= cell(num_files, 1);
all_state = cell(num_files, 1);

for file_i = 1:num_files
    file = char(folder_names(file_i));
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
    
        [chopt1, chopt2] = chop_by_file(file, dff_rz);
        Y = wsp_r(chopt1:chopt2);

    all_wsp{file_i} = rescale(wsp_r);  %
    all_Loco{file_i} = rescale(MI_wheel_r);
    all_Wsk{file_i} = rescale(MI_whisker_r);
    all_state{file_i} = rescale(dis_L2);
 
end



%%
N = numel(all_wsp);
for i = 1:N
    wsp   = all_wsp{i}(:);
    loco  = all_Loco{i}(:);
    wsk   = all_Wsk{i}(:);
    state = all_state{i}(:);

    L = min([length(wsp), length(loco), length(wsk), length(state)]);
    wsp   = wsp(1:L);
    loco  = loco(1:L);
    wsk   = wsk(1:L);
    state = state(1:L);
    t     = 1:L;

    figure('Color','w', 'Name', sprintf('File %d', i), 'Position', [100 100 1000 600]);
    subplot(4,1,1);
    plot(t, wsp, 'r', 'LineWidth', 1.5); ylabel('WSP'); title(sprintf('File %d', i), 'FontWeight', 'normal');
    ylim([0 1]); grid on;

    subplot(4,1,2);
    plot(t, loco, 'b', 'LineWidth', 1.5); ylabel('Loco');
    ylim([0 1]); grid on;

    subplot(4,1,3);
    plot(t, wsk, 'g', 'LineWidth', 1.5); ylabel('WSK');
    ylim([0 1]); grid on;

    subplot(4,1,4);
    plot(t, state, 'k', 'LineWidth', 1.5); ylabel('State');
    ylim([0 1]); xlabel('Time (frames)'); grid on;
end

