clear all;close all;clc


concat_folder = 'X:\MFB\MFB_AH_2023\Correlation_data\4mice';
confiles = dir(fullfile(concat_folder, 'concat_*.mat'));

folder_names = {
    'Animal1';
    'Animal2';
    'Animal3';
    'Animal4';
};

%% fig. 6a  
file = '191209_13_44_12';
quickAnalysis;

Xall = xyz(:,1);
Yall = xyz(:,2);
Zall = xyz(:,3);

block_size = 100;
dff1 = dff0;
Nmf1 = length(xyz);
dff_r = reshape(permute(dff1, [1,3,2]), Nmf1, []);
dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';

[cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', block_size, 300);

sig_level = 2;

pids = find(zc_wL > sig_level);
nids = find(zc_wL < -sig_lev);
nsids = find(abs(zc_wL) < sig_level);

figure('Position',[360, 278, 900, 900]*.5)
hold on; rotate3d on;
grid on;

pids = logical(zc_wL>sig_lev);
nids = logical(zc_wL<-sig_lev);
nsids = logical(abs(zc_wL)<sig_lev);

h1 = scatter3(Xall(pids), Yall(pids), Zall(pids), 20, 'r', 'LineWidth', 2);
h2 = scatter3(Xall(nids), Yall(nids), Zall(nids), 20, 'b', 'LineWidth', 2);
h3 = scatter3(Xall(nsids), Yall(nsids), Zall(nsids), 20, 'MarkerEdgeColor', [0.3 0.3 0.3],'LineWidth', 2);

legend([h1,h2,h3], {'MFB+', 'MFB-', 'MFB*'}, 'Location','northwest')
legend boxoff

xyz_mx = [250,250, max(Zall)];
xyz_mn = [0,0, min(Zall)];
xlim([0,xyz_mx(1)]);
ylim([0,xyz_mx(2)]);
zlim([xyz_mn(3)-10,xyz_mx(3)+10]);
zticks([-100, -80, -60, -40, -20, 0]);
xlabel('X (μm)'); ylabel('Y (μm)'); zlabel('Z (μm)');
view(-28,28);

set(gca, 'LineWidth', 1, 'FontSize', 18, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'   , ...
    'ZMinorTick'  , 'off'  )

currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure6', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end
fileName = ['Spatial_example_locomotion'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig.6b PM vs NM spatial distribution (permutation test + ECDF, per animal)

file_names = {
    '171212_16_19_37';
    '191018_13_39_41';
    '191209_13_44_12';
    '200130_14_15_24 FunctAcq';
};
sig_level  = 2;
n_shuffles = 10000;
dim_names  = {'X', 'Y', 'Z'};
n_animals  = length(file_names);

% Colormap for MFB+ / MFB-
redToWhite  = [linspace(1, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)'];
whiteToBlue = [linspace(1, 0, 128)', linspace(1, 0, 128)', linspace(1, 1, 128)'];
mycm   = flipud([redToWhite; whiteToBlue]);
cl_on  = mycm(256 - 40, :);  % MFB+
cl_off = mycm(40, :);         % MFB-

% Results table and storage for coordinates
perm_results = table('Size', [n_animals, 6], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'File', 'nPM', 'nNM', 'P_X', 'P_Y', 'P_Z'});
xyz_pm = cell(n_animals, 1);
xyz_nm = cell(n_animals, 1);

% ===== Pass 1: load data, run permutation tests =====
rng(42);
for iii = 1:n_animals
    file = file_names{iii};
    fprintf('Processing %s...\n', file);
    quickAnalysis;
    close all;  % close figures opened by quickAnalysis
    [~, zc_wL] = bootstrap_cc(dff0_rz', L_state', 100, 300);

    pm_idx = find(zc_wL >  sig_level);
    nm_idx = find(zc_wL < -sig_level);

    xyz_pm{iii} = xyz(pm_idx, :);
    xyz_nm{iii} = xyz(nm_idx, :);

    perm_results.File(iii) = string(file);
    perm_results.nPM(iii)  = length(pm_idx);
    perm_results.nNM(iii)  = length(nm_idx);

    for d = 1:3
        pos_pm = xyz(pm_idx, d);
        pos_nm = xyz(nm_idx, d);

        [~, ~, d_obs] = kstest2(pos_pm, pos_nm);
        pool    = [pos_pm; pos_nm];
        n_pm    = length(pos_pm);
        d_shuff = zeros(n_shuffles, 1);
        for s = 1:n_shuffles
            r = randperm(length(pool));
            [~, ~, d_shuff(s)] = kstest2(pool(r(1:n_pm)), pool(r(n_pm+1:end)));
        end
        p_perm = sum(d_shuff >= d_obs) / n_shuffles;

        switch d
            case 1, perm_results.P_X(iii) = p_perm;
            case 2, perm_results.P_Y(iii) = p_perm;
            case 3, perm_results.P_Z(iii) = p_perm;
        end
    end
    fprintf('  P: X=%.4f, Y=%.4f, Z=%.4f\n', ...
        perm_results.P_X(iii), perm_results.P_Y(iii), perm_results.P_Z(iii));
end

disp('--- Permutation test results (MFB+ vs MFB-, 10000 shuffles) ---');
disp(perm_results);

% ===== Pass 2: plot all results in one figure =====
figure('Position', [50, 50, 1100, 500]);

for iii = 1:n_animals
    p_vals = [perm_results.P_X(iii), perm_results.P_Y(iii), perm_results.P_Z(iii)];

    for d = 1:3
        sp_idx = (d - 1) * n_animals + iii;
        subplot(3, n_animals, sp_idx); hold on;

        [f_pm, x_pm] = ecdf(xyz_pm{iii}(:, d));
        [f_nm, x_nm] = ecdf(xyz_nm{iii}(:, d));
        plot(x_pm, f_pm, '-', 'Color', cl_on,  'LineWidth', 2);
        plot(x_nm, f_nm, '-', 'Color', cl_off, 'LineWidth', 2);

        % Permutation p value
        p_perm = p_vals(d);
        if p_perm == 0
            p_str = sprintf('p<%.2f', 1/n_shuffles);
        elseif p_perm < 0.01
            p_str = 'p<0.01';
        else
            p_str = sprintf('p=%.2f', p_perm);
        end
        text(0.95, 0.05, p_str, 'Units', 'normalized', ...
            'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

        if d == 1
            title(sprintf('Animal %d', iii), 'FontSize', 18);
        end
        if iii == 1
            ylabel(dim_names{d}, 'FontSize', 18);
        end
        if d == 3
            xlabel('Position (μm)', 'FontSize', 18);
        end
        if d == 1 && iii == n_animals
            h1 = plot(nan, nan, '-', 'Color', cl_on,  'LineWidth', 2);
            h2 = plot(nan, nan, '-', 'Color', cl_off, 'LineWidth', 2);
            legend([h1, h2], {'MFB+', 'MFB-'}, 'Location', 'best', 'FontSize', 14);
            legend box off;
        end

        set(gca, 'LineWidth', 1, 'FontSize', 18, ...
            'FontName'   , 'Helvetica', ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.02 .02] );
        hold off;
    end
end

fileName = 'PMNMXYZ_ECDF';
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

% ===== Pass 3: pooled histogram (all animals, X/Y/Z) =====
% Concatenate across animals
all_pm = cat(1, xyz_pm{:});
all_nm = cat(1, xyz_nm{:});

figure('Position', [50, 50, 400, 500]);

for d = 1:3
    subplot(3, 1, d); hold on;

    % Determine shared bin edges
    all_vals = [all_pm(:, d); all_nm(:, d)];
    edges = linspace(min(all_vals), max(all_vals), 25);

    h1 = histogram(all_pm(:, d), edges, 'FaceColor', cl_on,  'EdgeColor', 'none', 'FaceAlpha', 0.6);
    h2 = histogram(all_nm(:, d), edges, 'FaceColor', cl_off, 'EdgeColor', 'none', 'FaceAlpha', 0.6);

    ylabel('Count', 'FontSize', 18);
    if d == 3
        xlabel('Position (μm)', 'FontSize', 18);
    end
    if d == 1
        lg = legend([h1, h2], {'MFB+', 'MFB-'}, 'Location', 'best', 'FontSize', 14);
        lg.Box = 'off';
        lg.ItemTokenSize = [15, 10];
    end
    title(dim_names{d}, 'FontSize', 18);

    set(gca, 'LineWidth', 1, 'FontSize', 18, ...
        'FontName'   , 'Helvetica', ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] );
    hold off;
end

fileName = 'PMNMXYZ_Hist_pooled';
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig 6c.Ripley function analysis

% Parameters
sig_level  = 2.0;
fov_size   = [250, 250];
N_perm     = 1000;
cell_types = {'PM', 'NM', 'NSM'};
type_colors = [1.0000  0.3150  0.3150;
               0.3071  0.3071  1.0000;
               0.7000  0.7000  0.7000];

% Containers: store full permutation matrices for correct CI
all_L_obs     = struct('PM', {{}}, 'NM', {{}}, 'NSM', {{}});
all_L_perm    = struct('PM', {{}}, 'NM', {{}}, 'NSM', {{}});
all_weights   = struct('PM', {{}}, 'NM', {{}}, 'NSM', {{}});
r_ref = [];

rng(42);

% File loop
for file_i = 1:length(confiles)
    file     = folder_names{file_i};
    filePath = fullfile(concat_folder, confiles(file_i).name);

    fprintf('\n=== File %d/%d: %s ===\n', file_i, length(confiles), file);

    load(filePath);

    if ~exist('dff_r0_all','var') || ~exist('xyz','var') || ~exist('L_state','var')
        warning('Missing required variables, skipping.');
        continue;
    end

    % Compute correlation-based classification
    dff_r  = dff_r0_all;
    dff_rz = (dff_r - nanmean(dff_r, 2)) ./ nanstd(dff_r, 0, 2);
    [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', 100, 300);

    pm_idx  = find(zc_wL >  sig_level);
    nm_idx  = find(zc_wL < -sig_level);
    nsm_idx = find(abs(zc_wL) < sig_level);
    idx_sets = {pm_idx, nm_idx, nsm_idx};

    % Analyze per z-layer
    z_vals = xyz(:, 3);
    layers = unique(z_vals);

    for li = 1:length(layers)
        z_now     = layers(li);
        layer_idx = find(z_vals == z_now);

        for ct = 1:3
            idx_in_layer = intersect(layer_idx, idx_sets{ct});
            n_pts = length(idx_in_layer);

            if n_pts < 10 || length(layer_idx) < n_pts + 2
                continue;
            end

            xy_obs = xyz(idx_in_layer, 1:2);

            % Observed L(r)
            try
                [r_vals, ~, L_obs] = ripleyKL2D(xy_obs, 0, fov_size);
                if isempty(r_ref), r_ref = r_vals; end
            catch
                continue;
            end

            % Null distribution: randomly sample n_pts from same layer
            L_perm = zeros(N_perm, length(r_vals));
            for p = 1:N_perm
                rand_idx = layer_idx(randperm(length(layer_idx), n_pts));
                xy_rand  = xyz(rand_idx, 1:2);
                try
                    [~, ~, L_tmp] = ripleyKL2D(xy_rand, 0, fov_size);
                    L_perm(p, :) = L_tmp;
                catch
                    L_perm(p, :) = NaN;
                end
            end

            % Remove failed permutations
            L_perm = L_perm(~any(isnan(L_perm), 2), :);
            if size(L_perm, 1) < 10, continue; end

            % Store observed and full permutation matrix
            all_L_obs.(cell_types{ct}){end+1}  = L_obs;
            all_L_perm.(cell_types{ct}){end+1} = L_perm;   % N_valid x n_r
            all_weights.(cell_types{ct}){end+1} = n_pts;

            fprintf('  [%s] z=%g, n=%d, perm OK\n', cell_types{ct}, z_now, n_pts);
        end
    end

    clearvars -except ...
        confiles folder_names concat_folder ...
        sig_level fov_size N_perm cell_types type_colors ...
        all_L_obs all_L_perm all_weights r_ref file_i;
end

%% Plotting
figure('Position', [100 100 750 320]);
labels_names = {'MFB+', 'MFB-'};
plot_types = [1, 2];  % PM and NM only

for pi = 1:2
    ct = plot_types(pi);
    subplot(1, 2, pi);
    hold on;

    L_obs_cell = all_L_obs.(cell_types{ct});
    if isempty(L_obs_cell), continue; end

    % Get weights
    w = cell2mat(all_weights.(cell_types{ct})');  % n_layers x 1
    w = w / sum(w);
    n_layers = length(w);
    n_r = length(r_ref);

    % Observed: weighted mean
    L_obs_mat  = cell2mat(L_obs_cell');
    L_obs_mean = w' * L_obs_mat;

    % Null: aggregate first, then 95% CI from percentiles
    L_perm_cells = all_L_perm.(cell_types{ct});
    min_perm = min(cellfun(@(x) size(x, 1), L_perm_cells));

    L_agg_null = zeros(min_perm, n_r);
    for k = 1:n_layers
        L_agg_null = L_agg_null + w(k) * L_perm_cells{k}(1:min_perm, :);
    end

    L_null_mean = mean(L_agg_null, 1);
    L_null_lo   = prctile(L_agg_null, 2.5, 1);
    L_null_hi   = prctile(L_agg_null, 97.5, 1);

    c = type_colors(ct, :);

    % Null 95% CI (grey shading)
    fill([r_ref, fliplr(r_ref)], ...
         [L_null_hi, fliplr(L_null_lo)], ...
         [0.85 0.85 0.85], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    % Null mean (grey dashed)
    plot(r_ref, L_null_mean, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

    % Observed mean curve
    plot(r_ref, L_obs_mean, '-', 'Color', c, 'LineWidth', 2.5);

    % Mark regions where observed mean exceeds null 95% CI
    above = L_obs_mean > L_null_hi;
    if any(above)
        plot(r_ref(above), L_obs_mean(above), '.', 'Color', c, 'MarkerSize', 8);
    end

    xlabel('r (µm)', 'FontSize', 18);
    ylabel('L(r) - r', 'FontSize', 18);
    title(sprintf('%s', labels_names{pi}), 'FontSize', 18, 'FontWeight', 'bold');

    legend({'Null 95% CI', 'Null mean', ...
            sprintf('%s mean', labels_names{pi})}, ...
           'Location', 'best', 'Box', 'off', 'FontSize', 14);

    set(gca, 'LineWidth', 1, 'FontSize', 18, 'Box', 'off');
end

% Save
savepath    = 'X:/MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
savepath2   = fullfile(savepath, 'Figures', 'Figure6', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end
fullFilePathPDF = fullfile(savepath2, 'ripleyresults.pdf');
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% supplementary part
% --- ensure required params/flags exist even if this cell is run standalone ---
block_size = 100;
sig_level  = 2;
plot_correlationsSpatial      = true;
plot_correlationsSpatial_PMNM = true;

% initialise accumulators (appended to inside the file loop below)
cc_all_qw = []; dd_all_qw = [];
cc_all_as = []; dd_all_as = [];
cc_all_pm = []; dd_all_pm = [];
cc_all_nm = []; dd_all_nm = [];
cc_all_ns = []; dd_all_ns = [];

for file_i = 1:length(folder_names)
    fprintf('=== Processing file %d/%d: %s ===\n', file_i, length(folder_names), folder_names{file_i});
    
    filePath = fullfile(concat_folder, confiles(file_i).name);
    load(filePath)
    dff_r = dff_r0_all;
    cc1 = corr(dff_r', 'rows', 'complete');
    cc1(eye(length(cc1)) == 1) = nan;
    Nmf1 = length(cc1);
    dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';

    if ~all(L_state == 0)
        fprintf('  Running bootstrap_cc...\n');
        [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', block_size, 300);
    else
        fprintf('  Skipping: L_state is all zeros.\n');
        continue
    end

    pm = zc_wL > sig_level;
    [~, pm_ids] = find(pm == 1);
    nm = zc_wL < -sig_level;
    [~, nm_ids] = find(nm == 1);
    nsm = abs(zc_wL) < sig_level;
    [~, ns_ids] = find(nsm == 1);
    fprintf('  PM: %d, NM: %d, NSM: %d boutons\n', length(pm_ids), length(nm_ids), length(ns_ids));

    Nmf = Nmf1;

    if Nmf > 1
        fprintf('  Computing pairwise distances (Nmf = %d)...\n', Nmf);
        dd = nan * ones(Nmf, Nmf);
        for iii = 1:Nmf
            dd(iii,:) = sqrt(nansum((xyz - xyz(iii,:)).^2, 2));
        end

        % === QW ===
        fprintf('  Running QW bootstrap_cc_pw...\n');
        [cc_qw, zc_qw] = bootstrap_cc_pw(dff_rz(:, L_state == 0), block_size, 200);
        if length(cc_qw) > 1
            cc_qw(eye(size(cc_qw)) == 1) = nan;
            cc_qw(abs(zc_qw) < sig_level) = nan;
            cc_all_qw = [cc_all_qw; cc_qw(~isnan(cc_qw))];
            dd_all_qw = [dd_all_qw; dd(~isnan(cc_qw))];
            fprintf('  QW: %d pairs added\n', sum(~isnan(cc_qw(:))));
        end

        % === AS ===
        fprintf('  Running AS bootstrap_cc_pw...\n');
        [cc_as, zc_as] = bootstrap_cc_pw(dff_rz(:, L_state == 1), block_size, 200);
        if length(cc_as) > 1
            cc_as(eye(size(cc_as)) == 1) = nan;
            cc_as(abs(zc_as) < sig_level) = nan;
            cc_all_as = [cc_all_as; cc_as(~isnan(cc_as))];
            dd_all_as = [dd_all_as; dd(~isnan(cc_as))];
            fprintf('  AS: %d pairs added\n', sum(~isnan(cc_as(:))));
        end

        % === PM ===
        if ~isempty(pm_ids)
            fprintf('  Running PM bootstrap_cc_pw...\n');
            dd_pm = dd(pm_ids, pm_ids);
            [cc_pm, zc_pm] = bootstrap_cc_pw(dff_rz(pm_ids, L_state == 1), block_size, 200);
            if length(cc_pm) > 1
                cc_pm(eye(size(cc_pm)) == 1) = nan;
                cc_pm(abs(zc_pm) < sig_level) = nan;
                cc_all_pm = [cc_all_pm; cc_pm(~isnan(cc_pm))];
                dd_all_pm = [dd_all_pm; dd_pm(~isnan(cc_pm))];
                fprintf('  PM: %d pairs added\n', sum(~isnan(cc_pm(:))));
            end
        end

        % === NM ===
        if ~isempty(nm_ids)
            fprintf('  Running NM bootstrap_cc_pw...\n');
            dd_nm = dd(nm_ids, nm_ids);
            [cc_nm, zc_nm] = bootstrap_cc_pw(dff_rz(nm_ids, L_state == 1), block_size, 200);
            if length(cc_nm) > 1
                cc_nm(eye(size(cc_nm)) == 1) = nan;
                cc_nm(abs(zc_nm) < sig_level) = nan;
                cc_all_nm = [cc_all_nm; cc_nm(~isnan(cc_nm))];
                dd_all_nm = [dd_all_nm; dd_nm(~isnan(cc_nm))];
                fprintf('  NM: %d pairs added\n', sum(~isnan(cc_nm(:))));
            end
        end

        % === NSM ===
        if ~isempty(ns_ids)
            fprintf('  Running NSM bootstrap_cc_pw...\n');
            dd_ns = dd(ns_ids, ns_ids);
            [cc_ns, zc_ns] = bootstrap_cc_pw(dff_rz(ns_ids, L_state == 1), block_size, 200);
            if length(cc_ns) > 1
                cc_ns(eye(size(cc_ns)) == 1) = nan;
                cc_ns(abs(zc_ns) < sig_level) = nan;
                cc_all_ns = [cc_all_ns; cc_ns(~isnan(cc_ns))];
                dd_all_ns = [dd_all_ns; dd_ns(~isnan(cc_ns))];
                fprintf('  NSM: %d pairs added\n', sum(~isnan(cc_ns(:))));
            end
        end

    else
        fprintf('  Skipping: Nmf <= 1\n');
    end

    fprintf('  File %d/%d done.\n\n', file_i, length(folder_names));
end
fprintf('All files processed.\n');

%% === Plot: QW and AS spatial correlations ===
if plot_correlationsSpatial

    for kk = 1:2

        if kk == 1
            x = dd_all_qw;
            y = cc_all_qw;
            which_state = 'Quiet Wakefulness';
            title_color = 'c';
        elseif kk == 2
            x = dd_all_as;
            y = cc_all_as;
            which_state = 'Active State';
            title_color = 'm';
        end

        figure('Position', [160 360 600 400])
        hold on
        title(which_state, 'Color', title_color)

        xx = 0:10:200;
        clids = find((x > 10) .* (y < 1));
        x = x(clids);
        y = y(clids);

        bins = -1:.1:1;
        yhh = nan([length(bins)-1, length(xx)-1]);

        ym = zeros(1, length(xx)-1);
        ye = zeros(1, length(xx)-1);
        ym2 = zeros(1, length(xx)-1);
        ye2 = zeros(1, length(xx)-1);

        y_all_kk = cell(1, length(xx)-1);
        for j = 1:length(xx)-1
            bin_ids = find((x >= xx(j)) .* (x < xx(j+1)));
            y_all_kk{j} = y(bin_ids);

            ym(j) = nanmean(y(bin_ids));
            ye(j) = nanstd(y(bin_ids));
            ym2(j) = nanmean(abs(y(bin_ids)));
            ye2(j) = nanstd(abs(y(bin_ids)));

            [yh, yx] = histcounts(y(bin_ids), bins);
            yhh(:, j) = yh;
        end

        xc = xx(1:end-1) + diff(xx(1:2))/2;
        yxc = yx(1:end-1) + diff(yx(1:2))/2;

        if kk == 1
            ym_sta = ym; ye_sta = ye;
            ym2_sta = ym2; ye2_sta = ye2;
        elseif kk == 2
            ym_run = ym; ye_run = ye;
            ym2_run = ym2; ye2_run = ye2;
        end

        [a, b] = meshgrid(yxc, xc);
        surfc(a, b, yhh');
        colormap('redblue')
        rotate3d on
        view(-15, 60)

        ylabel({'MFB pairwise', 'distance (μm)'}, 'Rotation', -60)
        xlabel({'MFB pairwise correlation'}, 'rotation', 0)
        zlabel('# MFB pairs')

        set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
            'TickDir', 'out', 'TickLength', [.01 .01])

        currentDate = datestr(now, 'yyyy-mm-dd');
        savepath2 = fullfile(savepath, 'Figures', 'Figure6', currentDate);
        if ~exist(savepath2, 'dir')
            mkdir(savepath2);
        end
        fileName = ['CC_spatial_' which_state];
        fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
        exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    end

    %% Average CC plot
    figure('Position', [160 260 700 350])

    subplot(1,2,1);
    title('Q-state', 'color', 'c')
    hold on
    h1 = plot(xc, ym_sta, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_sta, ye_sta, 'k-')
    h2 = plot(xc, ym2_sta, 'r-', 'LineWidth', 2);
    errorbar(xc, ym2_sta, ye2_sta, 'r-')
    legend([h1, h2], {'Average CC', 'Average |CC|'}, 'Location', 'best')
    legend boxoff
    xlim([0, 200]); ylim([-1, 1])
    xlabel({'MFB pairwise distance (μm)'})
    ylabel({'MFB pairwise correlation'})
    set(gca, 'LineWidth', 1, 'FontSize', 18, 'Box', 'off', ...
        'TickDir', 'out', 'TickLength', [.01 .01])

    subplot(1,2,2);
    title('A-State', 'color', 'm')
    hold on
    plot(xc, ym_run, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_run, ye_run, 'k-')
    plot(xc, ym2_run, 'r-', 'LineWidth', 2);
    errorbar(xc, ym2_run, ye2_run, 'r-')
    xlim([0, 200]); ylim([-1, 1])
    set(gca, 'LineWidth', 1, 'FontSize', 18, 'Box', 'off', ...
        'TickDir', 'out', 'TickLength', [.01 .01])

    fileName = 'CC_spatial_avg_AQ';
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
end

%% === Plot: PM / NM / NSM ===
if plot_correlationsSpatial_PMNM

    for kk = 1:3

        if kk == 1
            x = dd_all_pm; y = cc_all_pm; which_state = 'PM';
        elseif kk == 2
            x = dd_all_nm; y = cc_all_nm; which_state = 'NM';
        elseif kk == 3
            x = dd_all_ns; y = cc_all_ns; which_state = 'NSM';
        end

        figure('Position', [160 360 600 400])
        hold on
        title(which_state)

        xx = 0:10:200;
        clids = find((x > 10) .* (y < 1));
        x = x(clids);
        y = y(clids);

        bins = -1:.1:1;
        yhh = nan([length(bins)-1, length(xx)-1]);

        ym = zeros(1, length(xx)-1);
        ye = zeros(1, length(xx)-1);
        ym2 = zeros(1, length(xx)-1);
        ye2 = zeros(1, length(xx)-1);

        for j = 1:length(xx)-1
            bin_ids = find((x >= xx(j)) .* (x < xx(j+1)));
            ym(j) = nanmean(y(bin_ids));
            ye(j) = nanstd(y(bin_ids));
            ym2(j) = nanmean(abs(y(bin_ids)));
            ye2(j) = nanstd(abs(y(bin_ids)));
            [yh, yx] = histcounts(y(bin_ids), bins);
            yhh(:, j) = yh;
        end

        xc = xx(1:end-1) + diff(xx(1:2))/2;
        yxc = yx(1:end-1) + diff(yx(1:2))/2;

        if kk == 1
            ym_pm = ym; ye_pm = ye; ym2_pm = ym2; ye2_pm = ye2;
        elseif kk == 2
            ym_nm = ym; ye_nm = ye; ym2_nm = ym2; ye2_nm = ye2;
        elseif kk == 3
            ym_ns = ym; ye_ns = ye; ym2_ns = ym2; ye2_ns = ye2;
        end

        [a, b] = meshgrid(yxc, xc);
        surfc(a, b, yhh');
        colormap('redblue')
        rotate3d on
        view(-15, 60)

        ylabel({'MFB pairwise', 'distance (μm)'}, 'Rotation', -60)
        xlabel({'MFB pairwise correlation'}, 'rotation', 0)
        zlabel('# MFB pairs')

        set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
            'TickDir', 'out', 'TickLength', [.01 .01])

        fileName = ['CC_spatial_' which_state];
        fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
        exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    end

    figure('Position', [160 60 350 750])

    subplot(311);
    title('PM'); hold on
    h1 = plot(xc, ym_pm, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_pm, ye_pm, 'k-')
    h2 = plot(xc, ym2_pm, 'r-', 'LineWidth', 2);
    errorbar(xc, ym2_pm, ye2_pm, 'r-')
    legend([h1, h2], {'Average CC', 'Average |CC|'}, 'Location', 'best')
    legend boxoff
    xlim([0, 200]); ylim([-1, 1])
    pos = get(gca, 'Position'); pos(1) = 0.2; pos(3) = 0.75; set(gca, 'Position', pos)
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01])

    subplot(312);
    title('NM'); hold on
    h3 = plot(xc, ym_nm, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_nm, ye_nm, 'k-')
    h4 = plot(xc, ym2_nm, 'b-', 'LineWidth', 2);
    errorbar(xc, ym2_nm, ye2_nm, 'b-')
    legend([h3, h4], {'Average CC', 'Average |CC|'}, 'Location', 'best')
    legend boxoff
    xlim([0, 200]); ylim([-1, 1])
    ylabel({'MFB pairwise correlation'})
    pos = get(gca, 'Position'); pos(1) = 0.2; pos(3) = 0.75; set(gca, 'Position', pos)
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01])

    subplot(313);
    title('NSM'); hold on
    h5 = plot(xc, ym_ns, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_ns, ye_ns, 'k-')
    h6 = plot(xc, ym2_ns, '-', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    errorbar(xc, ym2_ns, ye2_ns, '-', 'color', [0.5, 0.5, 0.5])
    legend([h5, h6], {'Average CC', 'Average |CC|'}, 'Location', 'best')
    legend boxoff
    xlim([0, 200]); ylim([-1, 1])
    xlabel({'MFB pairwise distance (μm)'})
    set(get(gca, 'YLabel'), 'Position', [-30 0 0])
    pos = get(gca, 'Position'); pos(1) = 0.2; pos(3) = 0.75; set(gca, 'Position', pos)
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01])

    fileName = 'CC_spatial_avg';
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
end