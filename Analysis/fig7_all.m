savepath = 'X:/MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure7', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end

%% fig. 7 plot shared PC
Nigel_c  = [1,.4,.4];
Jeremy_c = [.1,.5,.8];
Bernie_c = [1,.7,0];
Bill_c   = [.7,.1,.7];
folder = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\cvpca_raw_matlab';
files = dir(fullfile(folder, '*_grc_shared_dim_details.mat'));
nFiles = length(files);

% Nsub 
  if contains(files(1).name, 'grc', 'IgnoreCase', true)
      Nsub_values = 100:50:700;
      N_sub = [200:50:700];
  else
      Nsub_values = 100:25:175;
      N_sub = [100:25:175];
  end
nSub = length(Nsub_values);

figure('Position', [100, 100, 400, 375]);
hold on;

for i = 1:nFiles
    fname = files(i).name;
    parts = strsplit(fname, '_');
    file = strjoin(parts(1:4), '_');
    Mice = getMice(file);
    
    if strcmp(Mice, 'Nigel'),        face_color = Nigel_c;
    elseif strcmp(Mice, 'Jeremy'),    face_color = Jeremy_c;
    elseif strcmp(Mice, 'Bernie'),    face_color = Bernie_c;
    elseif strcmp(Mice, 'Bill'),      face_color = Bill_c;
    else,                             face_color = [0.5, 0.5, 0.5]; Mice = 'Unknown';
    end
    if contains(fname, 'grc', 'IgnoreCase', true)

    data = load(fullfile(folder, fname), 'varexp_Nsub550');
    varexp = data.varexp_Nsub550;
    else
            data = load(fullfile(folder, fname), 'varexp_Nsub125');
    varexp = data.varexp_Nsub125;
    end
    avg_line = nanmean(varexp, 1);
    n_valid  = sum(~isnan(varexp), 1);
    sem_line = nanstd(varexp, 0, 1) ./ sqrt(n_valid);
    x = 1:length(avg_line);
    
    valid = ~isnan(avg_line);
    x_valid = x(valid);
    
    fill([x_valid, fliplr(x_valid)], ...
        [avg_line(valid) + sem_line(valid), fliplr(avg_line(valid) - sem_line(valid))], ...
        face_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h = plot(x_valid, avg_line(valid), '-', 'Color', face_color, 'LineWidth', 3);
    
    [vm, dm] = max(avg_line);
    plot(dm, vm + 0.05, 'v', 'MarkerFaceColor', face_color, ...
        'MarkerEdgeColor', face_color, 'MarkerSize', 12);
end

box off;
xlabel('Num. of components');
ylabel('CVEV');
ylim([0, 0.8]);
yticks([0 0.4 0.8]);
xlim([0 350]);
set(gca, 'FontSize', 18, 'Box', 'off', 'LineWidth', 1);
hold off;

fileName = ['shared PC_curve'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


N_file = nFiles;
slope = zeros(size(N_sub));
sem = zeros(size(N_sub));

file_colors = zeros(N_file, 3);
for i = 1:N_file
    fname = files(i).name;
    parts = strsplit(fname, '_');
    file = strjoin(parts(1:4), '_');
    Mice = getMice(file);
    
    if strcmp(Mice, 'Nigel'),        file_colors(i,:) = Nigel_c;
    elseif strcmp(Mice, 'Jeremy'),    file_colors(i,:) = Jeremy_c;
    elseif strcmp(Mice, 'Bernie'),    file_colors(i,:) = Bernie_c;
    elseif strcmp(Mice, 'Bill'),      file_colors(i,:) = Bill_c;
    else,                             file_colors(i,:) = [0.5, 0.5, 0.5];
    end
end

slope_ind_all = nan(N_file, length(N_sub));

for k = 1:length(N_sub)
    varexp = cell(N_file, 1);
    varname = sprintf('varexp_Nsub%d', N_sub(k));
    
    for dataset_i = 1:N_file
        fname = files(dataset_i).name;
        data = load(fullfile(folder, fname), varname);
        if isfield(data, varname)
            varexp{dataset_i} = data.(varname);
        else
            varexp{dataset_i} = [];
        end
    end
    
    varmax = nan(N_file, 1);
    dimmax = nan(N_file, 1);
    slope_val = nan(N_file, 1);
    
    for dataset_i = 1:N_file
        if ~isempty(varexp{dataset_i})
            [varmax(dataset_i), dimmax(dataset_i)] = max(nanmean(varexp{dataset_i}, 1));
            val = varmax(dataset_i);
            if isnan(val) || isinf(val) || val == 0
                varmax(dataset_i) = nan;
                dimmax(dataset_i) = nan;
            end
        end
        slope_val(dataset_i) = dimmax(dataset_i) ./ varmax(dataset_i);
    end
    
    slope_ind = N_sub(k) ./ slope_val;
    slope_ind_all(:, k) = slope_ind;
    
    ix = find(~isnan(varmax));
    if ~isempty(ix)
        slope(k) = (varmax(ix)' * varmax(ix)) \ (varmax(ix)' * dimmax(ix));
    else
        slope(k) = nan;
    end
end

if contains(fname, 'grc', 'IgnoreCase', true)

    figure('Position', [30 20 580 350]);
else
    figure('Position', [30 20 300 350]);
end

if contains(fname, 'mf', 'IgnoreCase', true)
    mean_color = [62, 130, 65] / 255* 0.5+0.5;
else
    mean_color = [41, 28, 64] / 255 * 0.5+0.5;
end

bar(N_sub', N_sub' ./ slope', 'FaceColor', mean_color, 'EdgeColor', 'w', 'LineWidth', 1);
hold on;

% plot scatter
jitter_width = 15;
for dataset_i = 1:N_file
    jitter = 0;
    scatter(N_sub + jitter, slope_ind_all(dataset_i, :), 40, ...
        file_colors(dataset_i, :), 'filled', 'MarkerFaceAlpha', 0.7, ...
        'MarkerEdgeColor', 'none');
end

set(gca, 'FontSize', 18, 'Box', 'off', 'LineWidth', 1);

if contains(fname, 'grc', 'IgnoreCase', true)
    xlim([min(N_sub) - 50, max(N_sub) + 50]);
    xlabel('Num. of GrCs');
    ylabel('GrCs per dimension');
    hold off;
else
    xlim([min(N_sub) - 10, max(N_sub) + 10]);
    xlabel('Num. of MFAs');
    ylabel('MFAs per dimension');
    hold off;
end


fileName = ['shared PC_axon_dim'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% input combination

topo_folder = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\topology';
topo_files = dir(fullfile(topo_folder, '*.csv*'));
nFiles = length(topo_files);

patterns = {'NNNN', 'NNNS', 'PNNN', 'NNSS', 'PNNS', 'NSSS', 'PNSS', 'PPNN', 'SSSS', 'PSSS', 'PPNS', 'PPSS', 'PPPN', 'PPPS', 'PPPP'};
nPat = length(patterns);


frac_all = nan(nFiles, nPat);
theory_all = nan(nFiles, nPat);

for f = 1:nFiles
    T = readtable(fullfile(topo_folder, topo_files(f).name));
    
    labels = T.label;
    fraction_pct = T.fraction_pct;
    theory_pct = T.theory_pct;
    n_PM = T.n_PM;
    n_NM = T.n_NM;
    n_NSM = T.n_NSM;
    
    for r = 1:height(T)

        pat_str = [repmat('P', 1, n_PM(r)), repmat('N', 1, n_NM(r)), repmat('S', 1, n_NSM(r))];
        pat_sorted = pat_str;
        
        pat_idx = find(strcmp(patterns, pat_sorted));
        if ~isempty(pat_idx)
            frac_all(f, pat_idx) = fraction_pct(r);
            theory_all(f, pat_idx) = theory_pct(r);
        end
    end
end

colors = zeros(nPat, 3);
for i = 1:nPat
    numP = sum(patterns{i} == 'P');
    numN = sum(patterns{i} == 'N');
    numS = sum(patterns{i} == 'S');
    total = numP + numN + numS;
    rR = numP / total;
    bR = numN / total;
    gR = numS / total;
    colors(i, :) = [1, 0, 0] * rR + [0, 0, 1] * bR + [0.6, 0.6, 0.6] * gR;
end
colors_bar = colors * 0.5 + 0.45;  


frac_mean = nanmean(frac_all, 1);
frac_sem  = nanstd(frac_all, 0, 1) ./ sqrt(sum(~isnan(frac_all), 1));
theory_mean = nanmean(theory_all, 1);

figure('Position', [100, 100, 900, 400]);
hold on;

x = 1:nPat;


for i = 1:nPat
    bar(i, frac_mean(i), 0.7, 'FaceColor', colors_bar(i, :), 'EdgeColor', 'none');
end

% Error bar: SEM
errorbar(x, frac_mean, frac_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);

% Scatter
for f = 1:nFiles
    for i = 1:nPat
        if ~isnan(frac_all(f, i))
            scatter(i, frac_all(f, i), 30, colors(i, :), 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
        end
    end
end

plot(x, theory_mean, 'v', 'MarkerFaceColor', [1, 0.6, 0], ...
    'MarkerEdgeColor', [1, 0.6, 0], 'MarkerSize', 8);

set(gca, 'XTick', x, 'XTickLabel', patterns, 'FontSize', 12);
xtickangle(45);
xlim([0.3, nPat + 0.7]);
ylabel('Fraction (%)');
yticks([0 15 30])
set(gca, 'FontSize', 12, 'Box', 'off', 'LineWidth', 1);
hold off;

fileName = ['input_comp'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% Nsub 150 ,550

mf_color  = [62, 130, 65] / 255;
grc_color_150 = [204, 153, 255] / 255;
grc_color_550 = [100, 50, 150] / 255;  % purple

mf_scatter_c  = mf_color;
grc_scatter_150 = grc_color_150;
grc_scatter_550 = grc_color_550;

% load
folder = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\cvpca_raw_matlab';

% GrC Nsub150
grc_files = dir(fullfile(folder, '*_GrC_shared_dim_details.mat'));
grc_varmax_150 = [];
grc_dimmax_150 = [];
for i = 1:length(grc_files)
    data = load(fullfile(folder, grc_files(i).name), 'varexp_Nsub150');
    varexp = data.varexp_Nsub150;
    for k = 1:size(varexp, 1)
        row = varexp(k, :);
        if all(isnan(row)), continue; end
        [varmax_, dimmax_] = max(row);
        grc_varmax_150 = [grc_varmax_150; varmax_];
        grc_dimmax_150 = [grc_dimmax_150; dimmax_];
    end
end

% GrC Nsub550
grc_varmax_550 = [];
grc_dimmax_550 = [];
for i = 1:length(grc_files)
    data = load(fullfile(folder, grc_files(i).name), 'varexp_Nsub550');
    varexp = data.varexp_Nsub550;
    for k = 1:size(varexp, 1)
        row = varexp(k, :);
        if all(isnan(row)), continue; end
        [varmax_, dimmax_] = max(row);
        grc_varmax_550 = [grc_varmax_550; varmax_];
        grc_dimmax_550 = [grc_dimmax_550; dimmax_];
    end
end

% MF Nsub150
mf_files = dir(fullfile(folder, '*_MF_shared_dim_details_Nreps50.mat'));
mf_varmax_all = [];
mf_dimmax_all = [];
for i = 2:length(mf_files)
    data = load(fullfile(folder, mf_files(i).name), 'varexp_Nsub150');
    varexp = data.varexp_Nsub150;
    for k = 1:size(varexp, 1)
        row = varexp(k, :);
        if all(isnan(row)), continue; end
        [varmax_, dimmax_] = max(row);
        mf_varmax_all = [mf_varmax_all; varmax_];
        mf_dimmax_all = [mf_dimmax_all; dimmax_];
    end
end

figure('Position', [100, 100, 620, 400]);
hold on;

% --- GrC Nsub150 ---
ix_grc150 = find(~isnan(grc_varmax_150) & grc_varmax_150 > 0);
plot(grc_varmax_150(ix_grc150), grc_dimmax_150(ix_grc150), 'v', ...
    'MarkerFaceColor', grc_scatter_150, 'MarkerEdgeColor', grc_scatter_150, 'MarkerSize', 3);
slope_grc150 = (grc_varmax_150(ix_grc150)' * grc_varmax_150(ix_grc150)) \ ...
               (grc_varmax_150(ix_grc150)' * grc_dimmax_150(ix_grc150));
plot([0, 1], [0, 1] * slope_grc150, '-', 'Color', grc_color_150, 'LineWidth', 6);
%plot([0, 1], [1, 1] * slope_grc150, ':', 'Color', grc_color_150, 'LineWidth', 1.5);

% --- GrC Nsub550 ---
ix_grc550 = find(~isnan(grc_varmax_550) & grc_varmax_550 > 0);
plot(grc_varmax_550(ix_grc550), grc_dimmax_550(ix_grc550), 'v', ...
    'MarkerFaceColor', grc_scatter_550, 'MarkerEdgeColor', grc_scatter_550, 'MarkerSize', 3);
slope_grc550 = (grc_varmax_550(ix_grc550)' * grc_varmax_550(ix_grc550)) \ ...
               (grc_varmax_550(ix_grc550)' * grc_dimmax_550(ix_grc550));
plot([0, 1], [0, 1] * slope_grc550, '--', 'Color', grc_color_550, 'LineWidth', 2);
%plot([0, 1], [1, 1] * slope_grc550, ':', 'Color', grc_color_550, 'LineWidth', 1.5);

% --- MF Nsub150 ---
ix_mf = find(~isnan(mf_varmax_all) & mf_varmax_all > 0);
plot(mf_varmax_all(ix_mf), mf_dimmax_all(ix_mf), 'v', ...
    'MarkerFaceColor', mf_scatter_c, 'MarkerEdgeColor', mf_scatter_c, 'MarkerSize', 3);
slope_mf = (mf_varmax_all(ix_mf)' * mf_varmax_all(ix_mf)) \ ...
           (mf_varmax_all(ix_mf)' * mf_dimmax_all(ix_mf));
plot([0, 1], [0, 1] * slope_mf, '-', 'Color', mf_color, 'LineWidth', 3);
%plot([0, 1], [1, 1] * slope_mf, ':', 'Color', mf_color, 'LineWidth', 1.5);

xlabel('Max variance explained');
ylabel('Num. components');
set(gca, 'FontSize', 18, 'Box', 'off', 'LineWidth', 1);
yticks([0 100 200 300])

ax = gca;
xr = ax.XLim;
yr = ax.YLim;
tx = xr(2) * 0.65;
ty_start = yr(2) * 0.95;
gap = (yr(2) - yr(1)) * 0.07;

plot(tx - 0.02, ty_start, 'v', 'MarkerFaceColor', grc_color_150, 'MarkerEdgeColor', grc_color_150, 'MarkerSize', 8);
builtin('text', tx, ty_start, 'GrC n=150', 'Color', grc_color_150, 'FontSize', 14, 'FontWeight', 'bold');

plot(tx - 0.02, ty_start - gap, 'v', 'MarkerFaceColor', grc_color_550, 'MarkerEdgeColor', grc_color_550, 'MarkerSize', 8);
builtin('text', tx, ty_start - gap, 'GrC n=550', 'Color', grc_color_550, 'FontSize', 14, 'FontWeight', 'bold');

plot(tx - 0.02, ty_start - 2*gap, 'v', 'MarkerFaceColor', mf_color, 'MarkerEdgeColor', mf_color, 'MarkerSize', 8);
builtin('text', tx, ty_start - 2*gap, 'MF n=150', 'Color', mf_color, 'FontSize', 14, 'FontWeight', 'bold');
hold off;

fprintf('GrC N150 slope: %.2f (%.1f axons/dim)\n', slope_grc150, 150/slope_grc150);
fprintf('GrC N550 slope: %.2f (%.1f axons/dim)\n', slope_grc550, 550/slope_grc550);
fprintf('MF  N150 slope: %.2f (%.1f axons/dim)\n', slope_mf, 150/slope_mf);


fileName = ['N150VS550'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
%%
slopes = [slope_mf, slope_grc150, slope_grc550];
labels = {'MFA (n=150)', 'GrC (n=150)', 'GrC (n=550)'};
colors = [mf_color; grc_color_150; grc_color_550];

figure('Position', [150, 150, 450, 450]);
hold on;

for i = 1:3
    b = bar(i, slopes(i), 'FaceColor', colors(i, :), 'EdgeColor', 'none', 'BarWidth', 0.7);
%     text(i, slopes(i) + 9, sprintf('%.1f', slopes(i)), ...
%         'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');
end

ylabel('Extrapolated dimensionality');
set(gca, 'XTick', 1:3, 'XTickLabel', labels, 'FontSize', 18, 'Box', 'off', 'LineWidth', 1.5);

ylim([0, max(slopes) * 1.2]);

ax = gca;
ax.XGrid = 'off';

fileName = ['last figure'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% firing rate mean

summary_folder = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\summaries';

% MF
mf_files = dir(fullfile(summary_folder, 'mf', '*.json'));
mf_rate_on = [];
mf_rate_off = [];
mf_rate_all = [];
mf_session_id = {};
for i = 1:length(mf_files)
    data = jsondecode(fileread(fullfile(summary_folder, 'mf', mf_files(i).name)));
    mf_rate_on  = [mf_rate_on;  data.sparseness.ALL.rate_ON_mean];
    mf_rate_off = [mf_rate_off; data.sparseness.ALL.rate_OFF_mean];
    mf_rate_all = [mf_rate_all; data.sparseness.ALL.rate_ALL_mean];
    [~, fname, ~] = fileparts(mf_files(i).name);
    parts = strsplit(fname, '_');
    mf_session_id{end+1} = strjoin(parts(1:4), '_');
end

% GrC
grc_files = dir(fullfile(summary_folder, 'grc', '*.json'));
grc_rate_on = [];
grc_rate_off = [];
grc_rate_all = [];
grc_session_id = {};
for i = 1:length(grc_files)
    data = jsondecode(fileread(fullfile(summary_folder, 'grc', grc_files(i).name)));
    grc_rate_on  = [grc_rate_on;  data.sparseness.ALL.rate_ON_mean];
    grc_rate_off = [grc_rate_off; data.sparseness.ALL.rate_OFF_mean];
    grc_rate_all = [grc_rate_all; data.sparseness.ALL.rate_ALL_mean];
    [~, fname, ~] = fileparts(grc_files(i).name);
    parts = strsplit(fname, '_');
    grc_session_id{end+1} = strjoin(parts(1:4), '_');
end

mf_scatter  = [62, 130, 65] / 255;
mf_bar      = mf_scatter * 0.5 + 0.5;
grc_scatter = [41, 28, 64] / 255;
grc_bar     = grc_scatter * 0.5 + 0.5;

% SEM
mf_on_sem   = nanstd(mf_rate_on)  / sqrt(sum(~isnan(mf_rate_on)));
grc_on_sem  = nanstd(grc_rate_on) / sqrt(sum(~isnan(grc_rate_on)));
mf_off_sem  = nanstd(mf_rate_off)  / sqrt(sum(~isnan(mf_rate_off)));
grc_off_sem = nanstd(grc_rate_off) / sqrt(sum(~isnan(grc_rate_off)));
mf_all_sem  = nanstd(mf_rate_all)  / sqrt(sum(~isnan(mf_rate_all)));
grc_all_sem = nanstd(grc_rate_all) / sqrt(sum(~isnan(grc_rate_all)));

figure('Position', [100, 100, 700, 350]);
hold on;

% === ON ===
bar(1, nanmean(mf_rate_on), 0.6, 'FaceColor', mf_bar, 'EdgeColor', 'none');
bar(2, nanmean(grc_rate_on), 0.6, 'FaceColor', grc_bar, 'EdgeColor', 'none');
errorbar(1, nanmean(mf_rate_on), mf_on_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(2, nanmean(grc_rate_on), grc_on_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        plot([1, 2], [mf_rate_on(i), grc_rate_on(match_idx(1))], '--', ...
            'Color', [0.2, 0.2,0.2,0.5], 'LineWidth', 1);
    end
end

scatter(ones(length(mf_rate_on), 1), mf_rate_on, 50, mf_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
scatter(2 * ones(length(grc_rate_on), 1), grc_rate_on, 50, grc_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');

% === OFF ===
bar(4, nanmean(mf_rate_off), 0.6, 'FaceColor', mf_bar, 'EdgeColor', 'none');
bar(5, nanmean(grc_rate_off), 0.6, 'FaceColor', grc_bar, 'EdgeColor', 'none');
errorbar(4, nanmean(mf_rate_off), mf_off_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(5, nanmean(grc_rate_off), grc_off_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        plot([4, 5], [mf_rate_off(i), grc_rate_off(match_idx(1))], '--', ...
            'Color',[0.2, 0.2,0.2,0.5], 'LineWidth', 1);
    end
end

scatter(4 * ones(length(mf_rate_off), 1), mf_rate_off, 50, mf_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
scatter(5 * ones(length(grc_rate_off), 1), grc_rate_off, 50, grc_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');

% === ALL ===
bar(7, nanmean(mf_rate_all), 0.6, 'FaceColor', mf_bar, 'EdgeColor', 'none');
bar(8, nanmean(grc_rate_all), 0.6, 'FaceColor', grc_bar, 'EdgeColor', 'none');
errorbar(7, nanmean(mf_rate_all), mf_all_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(8, nanmean(grc_rate_all), grc_all_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        plot([7, 8], [mf_rate_all(i), grc_rate_all(match_idx(1))], '--', ...
            'Color', [0.2, 0.2,0.2,0.5], 'LineWidth', 1);
    end
end

scatter(7 * ones(length(mf_rate_all), 1), mf_rate_all, 50, mf_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
scatter(8 * ones(length(grc_rate_all), 1), grc_rate_all, 50, grc_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');

axis auto;
yrange = ylim;

plot([3, 3], yrange, 'k--', 'LineWidth', 1);
plot([6, 6], yrange, 'k--', 'LineWidth', 1);

text(1.5, yrange(2) * 0.95, 'A-state', 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(4.5, yrange(2) * 0.95, 'Q-state', 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(7.5, yrange(2) * 0.95, 'All',     'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

set(gca, 'XTick', [1, 2, 4, 5, 7, 8], 'XTickLabel', {'MFA', 'GrC', 'MFA', 'GrC', 'MFA', 'GrC'});
xlim([0.3, 8.7]);
ylim(yrange);
ylabel('Average firing rate (Hz)');
set(gca, 'FontSize', 18, 'Box', 'off', 'LineWidth', 1);

hold off;

% === Statistical Tests (Paired t-test) ===
mf_on_paired  = [];
grc_on_paired = [];
mf_off_paired  = [];
grc_off_paired = [];
mf_all_paired  = [];
grc_all_paired = [];

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        mf_on_paired  = [mf_on_paired;  mf_rate_on(i)];
        grc_on_paired = [grc_on_paired; grc_rate_on(match_idx(1))];
        mf_off_paired  = [mf_off_paired;  mf_rate_off(i)];
        grc_off_paired = [grc_off_paired; grc_rate_off(match_idx(1))];
        mf_all_paired  = [mf_all_paired;  mf_rate_all(i)];
        grc_all_paired = [grc_all_paired; grc_rate_all(match_idx(1))];
    end
end

[p_on,  h_on]  = ttest(mf_on_paired,  grc_on_paired);
[p_off, h_off] = ttest(mf_off_paired, grc_off_paired);
[p_all, h_all] = ttest(mf_all_paired, grc_all_paired);

fprintf('Paired t-test (MF vs GrC):\n');
fprintf('  A-state (ON):  p = %.4f, h = %d\n', p_on,  h_on);
fprintf('  Q-state (OFF): p = %.4f, h = %d\n', p_off, h_off);
fprintf('  All:           p = %.4f, h = %d\n', p_all, h_all);

fileName = ['FiringRate_mean'];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% sparseness

summary_folder = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\summaries';

% MF
mf_files = dir(fullfile(summary_folder, 'mf', '*.json'));
mf_sp_on = [];
mf_sp_off = [];
mf_sp_all = [];
mf_session_id = {};
for i = 1:length(mf_files)
    data = jsondecode(fileread(fullfile(summary_folder, 'mf', mf_files(i).name)));
    mf_sp_on  = [mf_sp_on; data.sparseness.ALL.pop_sparseness_ON];
    mf_sp_off = [mf_sp_off; data.sparseness.ALL.pop_sparseness_OFF];
    mf_sp_all = [mf_sp_all; data.sparseness.ALL.pop_sparseness_ALL];
    [~, fname, ~] = fileparts(mf_files(i).name);
    parts = strsplit(fname, '_');
    mf_session_id{end+1} = strjoin(parts(1:4), '_');
end

% GrC
grc_files = dir(fullfile(summary_folder, 'grc', '*.json'));
grc_sp_on = [];
grc_sp_off = [];
grc_sp_all = [];
grc_session_id = {};
for i = 1:length(grc_files)
    data = jsondecode(fileread(fullfile(summary_folder, 'grc', grc_files(i).name)));
    grc_sp_on  = [grc_sp_on; data.sparseness.ALL.pop_sparseness_ON];
    grc_sp_off = [grc_sp_off; data.sparseness.ALL.pop_sparseness_OFF];
    grc_sp_all = [grc_sp_all; data.sparseness.ALL.pop_sparseness_ALL];
    [~, fname, ~] = fileparts(grc_files(i).name);
    parts = strsplit(fname, '_');
    grc_session_id{end+1} = strjoin(parts(1:4), '_');
end

mf_scatter  = [62, 130, 65] / 255;
mf_bar      = mf_scatter * 0.5 + 0.5;
grc_scatter = [41, 28, 64] / 255;
grc_bar     = grc_scatter * 0.5 + 0.5;

% SEM
mf_on_sem   = nanstd(mf_sp_on) / sqrt(sum(~isnan(mf_sp_on)));
grc_on_sem  = nanstd(grc_sp_on) / sqrt(sum(~isnan(grc_sp_on)));
mf_off_sem  = nanstd(mf_sp_off) / sqrt(sum(~isnan(mf_sp_off)));
grc_off_sem = nanstd(grc_sp_off) / sqrt(sum(~isnan(grc_sp_off)));
mf_all_sem  = nanstd(mf_sp_all) / sqrt(sum(~isnan(mf_sp_all)));
grc_all_sem = nanstd(grc_sp_all) / sqrt(sum(~isnan(grc_sp_all)));

figure('Position', [100, 100, 700, 350]);
hold on;

% === ON ===
bar(1, nanmean(mf_sp_on), 0.6, 'FaceColor', mf_bar, 'EdgeColor', 'none');
bar(2, nanmean(grc_sp_on), 0.6, 'FaceColor', grc_bar, 'EdgeColor', 'none');
errorbar(1, nanmean(mf_sp_on), mf_on_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(2, nanmean(grc_sp_on), grc_on_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        plot([1, 2], [mf_sp_on(i), grc_sp_on(match_idx(1))], '--', ...
            'Color', [0.2, 0.2,0.2,0.5], 'LineWidth', 1);
    end
end

scatter(ones(length(mf_sp_on), 1), mf_sp_on, 50, mf_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
scatter(2 * ones(length(grc_sp_on), 1), grc_sp_on, 50, grc_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');

% === OFF ===
bar(4, nanmean(mf_sp_off), 0.6, 'FaceColor', mf_bar, 'EdgeColor', 'none');
bar(5, nanmean(grc_sp_off), 0.6, 'FaceColor', grc_bar, 'EdgeColor', 'none');
errorbar(4, nanmean(mf_sp_off), mf_off_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(5, nanmean(grc_sp_off), grc_off_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        plot([4, 5], [mf_sp_off(i), grc_sp_off(match_idx(1))], '--', ...
            'Color', [0.2, 0.2,0.2,0.5], 'LineWidth', 1);
    end
end

scatter(4 * ones(length(mf_sp_off), 1), mf_sp_off, 50, mf_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
scatter(5 * ones(length(grc_sp_off), 1), grc_sp_off, 50, grc_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');

% === ALL ===
bar(7, nanmean(mf_sp_all), 0.6, 'FaceColor', mf_bar, 'EdgeColor', 'none');
bar(8, nanmean(grc_sp_all), 0.6, 'FaceColor', grc_bar, 'EdgeColor', 'none');
errorbar(7, nanmean(mf_sp_all), mf_all_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(8, nanmean(grc_sp_all), grc_all_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        plot([7, 8], [mf_sp_all(i), grc_sp_all(match_idx(1))], '--', ...
            'Color', [0.2, 0.2,0.2,0.5], 'LineWidth', 1);
    end
end

scatter(7 * ones(length(mf_sp_all), 1), mf_sp_all, 50, mf_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
scatter(8 * ones(length(grc_sp_all), 1), grc_sp_all, 50, grc_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');


axis auto;
yrange = ylim;

plot([3, 3], yrange, 'k--', 'LineWidth', 1);
plot([6, 6], yrange, 'k--', 'LineWidth', 1);

text(1.5, yrange(2) * 0.95, 'A-state', 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(4.5, yrange(2) * 0.95, 'Q-state', 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(7.5, yrange(2) * 0.95, 'All',     'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

set(gca, 'XTick', [1, 2, 4, 5, 7, 8], 'XTickLabel', {'MFA', 'GrC', 'MFA', 'GrC', 'MFA', 'GrC'});
xlim([0.3, 8.7]);
ylim(yrange);
ylabel('Population sparseness');
set(gca, 'FontSize', 18, 'Box', 'off', 'LineWidth', 1);

hold off;

% === Statistical Tests (Wilcoxon Signed-Rank Test) ===
mf_on_paired  = [];
grc_on_paired = [];
mf_off_paired  = [];
grc_off_paired = [];
mf_all_paired  = [];
grc_all_paired = [];

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        mf_on_paired  = [mf_on_paired;  mf_sp_on(i)];
        grc_on_paired = [grc_on_paired; grc_sp_on(match_idx(1))];
        mf_off_paired  = [mf_off_paired;  mf_sp_off(i)];
        grc_off_paired = [grc_off_paired; grc_sp_off(match_idx(1))];
        mf_all_paired  = [mf_all_paired;  mf_sp_all(i)];
        grc_all_paired = [grc_all_paired; grc_sp_all(match_idx(1))];
    end
end

[p_on,  h_on]  = ttest(mf_on_paired,  grc_on_paired);
[p_off, h_off] = ttest(mf_off_paired, grc_off_paired);
[p_all, h_all] = ttest(mf_all_paired, grc_all_paired);

fprintf('Signed-rank test (MF vs GrC):\n');
fprintf('  A-state (ON):  p = %.4f, h = %d\n', p_on,  h_on);
fprintf('  Q-state (OFF): p = %.4f, h = %d\n', p_off, h_off);
fprintf('  All:           p = %.4f, h = %d\n', p_all, h_all);


fileName = ['Sparseness'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% effdim

summary_folder = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\summaries';
% MF
mf_files = dir(fullfile(summary_folder, 'mf', '*.json'));
mf_effdim = [];
mf_session_id = {};
for i = 1:length(mf_files)
    data = jsondecode(fileread(fullfile(summary_folder, 'mf', mf_files(i).name)));
    mf_effdim = [mf_effdim; data.effdim_all];
    [~, fname, ~] = fileparts(mf_files(i).name);
    parts = strsplit(fname, '_');
    mf_session_id{end+1} = strjoin(parts(1:4), '_');  % e.g., '171212_16_19_37'
end
% GrC
grc_files = dir(fullfile(summary_folder, 'grc', '*.json'));
grc_effdim = [];
grc_session_id = {};
for i = 1:length(grc_files)
    data = jsondecode(fileread(fullfile(summary_folder, 'grc', grc_files(i).name)));
    grc_effdim = [grc_effdim; data.effdim_all];
    [~, fname, ~] = fileparts(grc_files(i).name);
    parts = strsplit(fname, '_');
    grc_session_id{end+1} = strjoin(parts(1:4), '_');
end

mf_scatter  = [62, 130, 65] / 255;
mf_bar      = mf_scatter * 0.5 + 0.5;
grc_scatter = [41, 28, 64] / 255;
grc_bar     = grc_scatter * 0.5 + 0.5;

mf_sem  = nanstd(mf_effdim) / sqrt(sum(~isnan(mf_effdim)));
grc_sem = nanstd(grc_effdim) / sqrt(sum(~isnan(grc_effdim)));

figure('Position', [100, 100, 350, 350]);
hold on;

% Bar
bar(1, nanmean(mf_effdim), 0.5, 'FaceColor', mf_bar, 'EdgeColor', 'none');
bar(2, nanmean(grc_effdim), 0.5, 'FaceColor', grc_bar, 'EdgeColor', 'none');

% Error bar
errorbar(1, nanmean(mf_effdim), mf_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(2, nanmean(grc_effdim), grc_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        plot([1, 2], [mf_effdim(i), grc_effdim(match_idx(1))], '--', ...
            'Color', [0.2, 0.2,0.2,0.5], 'LineWidth', 1);
    end
end

% Scatter
scatter(ones(length(mf_effdim), 1), mf_effdim, 50, mf_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
scatter(2 * ones(length(grc_effdim), 1), grc_effdim, 50, grc_scatter, 'filled', ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');

set(gca, 'XTick', [1, 2], 'XTickLabel', {'MFA', 'GrC'});
xlim([0.3, 2.7]);
ylabel('Effective dimensionality');
set(gca, 'FontSize', 18, 'Box', 'off', 'LineWidth', 1);
hold off;
mf_effdim_paired  = [];
grc_effdim_paired = [];

for i = 1:length(mf_session_id)
    match_idx = find(strcmp(grc_session_id, mf_session_id{i}));
    if ~isempty(match_idx)
        mf_effdim_paired  = [mf_effdim_paired;  mf_effdim(i)];
        grc_effdim_paired = [grc_effdim_paired; grc_effdim(match_idx(1))];
    end
end

% === Signed-rank test ===
[p_effdim, h_effdim] = signrank(mf_effdim_paired, grc_effdim_paired);
fprintf('Signed-rank test (MF vs GrC) - Effective Dimensionality:\n');
fprintf('  p = %.4f, h = %d\n', p_effdim, h_effdim);



fileName = ['effdim'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig 7f eigenvalue


folderPath = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\eigenspectrum\grc';
filePattern = fullfile(folderPath, '*.mat');
matFiles = dir(filePattern);
numFiles = length(matFiles);
eiggen_20 = zeros(20, numFiles);
fprintf('Processing %d files...\n', numFiles);
for ix = 1:numFiles
    baseFileName = matFiles(ix).name;
    fullFileName = fullfile(folderPath, baseFileName);
    data = load(fullFileName, 'eigvals_norm');
        currentEig = data.eigvals_norm;
        len = min(20, length(currentEig));
        eiggen_20(1:len, ix) = currentEig(1:len);
end

if contains(folderPath, 'mf', 'IgnoreCase', true)
    mean_color = [62, 130, 65] / 255;
else
    mean_color = [41, 28, 64] / 255;
end

figure('Position',[200 100 350 400]);
hold on;
for ix = 1:size(eiggen_20,2)

    face_color = [0.7,0.7,0.7];

    plot(eiggen_20(:,ix),'color',[face_color],'linewidth',1.5)
end
mean_eigen = nanmean(eiggen_20,2);
sem_mean_eigen = nanstd(eiggen_20,0,2) / sqrt(size(eiggen_20,2));
plot(1:20, mean_eigen, '-', 'Color', mean_color, 'LineWidth', 2);
errorbar(1:20, mean_eigen, sem_mean_eigen, 'Color', mean_color, 'linestyle', 'none', 'LineWidth', 1);
xlabel('PCs')
ylabel('Normalised eigenvalue')
set(gca, 'FontSize', 18, 'XTick', [1 5 10 15 20], ...
'LineWidth', 1, 'YScale', 'log');
ylim([0.005, 1])
yticks([0.01 0.1 1])
yticklabels({'0.01','0.1','1'})
fileName = ['eigenspectrum'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig 7e. pie PM,NM,NSM
redToWhite = [linspace(1, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)'];
whiteToBlue = [linspace(1, 0, 128)', linspace(1, 0, 128)', linspace(1, 1, 128)'];
redToBlueMap = [redToWhite; whiteToBlue];
mycm = flipud(redToBlueMap);
cl_off = mycm(40, :);
cl_on = mycm(256 - 40, :);
cl_ns = [0.7, 0.7, 0.7];


folder = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\summaries\mf';
files = dir(fullfile(folder, '*.json'));
nFiles = length(files);

all_pm = zeros(nFiles, 1);
all_ns = zeros(nFiles, 1);
all_nm = zeros(nFiles, 1);

for i = 1:nFiles
    filepath = fullfile(folder, files(i).name);
    text = fileread(filepath);
    data = jsondecode(text);
    
    all_pm(i) = data.n_PM;
    all_ns(i) = data.n_NS;
    all_nm(i) = data.n_NM;
end

% Bar + Scatter + SEM Error Bar
figure('Position', [100, 100, 350, 350]);
hold on;

means = [mean(all_pm),mean(all_nm), mean(all_ns) ];
sems  = [std(all_pm) / sqrt(nFiles), std(all_nm) / sqrt(nFiles), std(all_ns) / sqrt(nFiles)];
colors = [cl_on;cl_off; cl_ns];
labels = {'MFA⁺', 'MFA⁻', 'MFA*'};

% Bar
for k = 1:3
    bar(k, means(k), 0.6, 'FaceColor', colors(k, :), 'EdgeColor', 'none');
end

% Error Bar (SEM)
errorbar(1:3, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

% Scatter
jitter_width = 0;
data_all = [all_pm,all_nm, all_ns];
for k = 1:3
    jitter = (rand(nFiles, 1) - 0.5) * jitter_width;
    scatter(k + jitter, data_all(:, k), 30, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
end

set(gca, 'XTick', 1:3, 'XTickLabel', labels, 'FontSize', 18);
ylabel('Neuron Count', 'FontSize', 18);
xlim([0.3, 3.7]);
box off;
set(gca, 'TickDir', 'out');
hold off;

fileName = ['individual'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


% pie
total_pm = sum(all_pm);
total_ns = sum(all_ns);
total_nm = sum(all_nm);

figure('Position', [100, 100, 350, 300]);
hold on;
axis off;

h = pie([total_nm, total_ns, total_pm]);
colormap([cl_off; cl_ns; cl_on]);
set(findobj(h, 'type', 'text'), 'fontsize', 18);
axis image;
set(gca, 'Position', [0.13, 0.15, 0.6, 0.6]);
%title(sprintf('Total: PM=%d, NM=%d, NS=%d', total_pm, total_nm, total_ns), 'FontSize', 12);

fprintf('\n=== summary ===\n');
fprintf('Total PM: %d\n', total_pm);
fprintf('Total NM: %d\n', total_nm);
fprintf('Total NS: %d\n', total_ns);
fprintf('Total N:  %d\n', total_pm + total_ns + total_nm);

fileName = ['pie'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% (PM-NM)/(PM+NM) bar+scatter for MF and GrC across animals

mf_color  = [62, 130, 65] / 255 * 0.5 + 0.5;
grc_color = [41, 28, 64] / 255 * 0.5 + 0.5;

folder_mf  = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\summaries\mf';
folder_grc = 'C:\Users\Yizhou Xie\Downloads\yizhou_handoff_mainfig_cvpca_2026-02-17\summaries\grc';

files_mf  = dir(fullfile(folder_mf,  '*.json'));
files_grc = dir(fullfile(folder_grc, '*.json'));

nAnimals = 4;
ratio_mf  = zeros(nAnimals, 1);
ratio_grc = zeros(nAnimals, 1);

for i = 1:nAnimals
    d = jsondecode(fileread(fullfile(folder_mf, files_mf(i).name)));
    ratio_mf(i) = (d.n_PM - d.n_NM) / (d.n_PM + d.n_NM);
    %ratio_mf(i) = (d.n_NM ) / (d.n_PM);
    d = jsondecode(fileread(fullfile(folder_grc, files_grc(i).name)));
    ratio_grc(i) = (d.n_PM - d.n_NM) / (d.n_PM + d.n_NM);
    %ratio_grc(i) = (d.n_NM ) / (d.n_PM);
end

% Plot
figure('Position', [100, 100, 600, 375]);
hold on;

group_gap = 0.3;  % gap between animals
bar_width = 0.55;

for i = 1:nAnimals
    x_mf  = (i-1) * (2 + group_gap) + 1;
    x_grc = (i-1) * (2 + group_gap) + 2;
    bar(x_mf,  ratio_mf(i),  bar_width, 'FaceColor', mf_color,  'EdgeColor', 'none');
    bar(x_grc, ratio_grc(i), bar_width, 'FaceColor', grc_color, 'EdgeColor', 'none');
end

% x tick positions and labels
xtick_pos    = zeros(1, nAnimals);
xtick_labels = cell(1, nAnimals);
for i = 1:nAnimals
    xtick_pos(i)    = (i-1) * (2 + group_gap) + 1.5;
    xtick_labels{i} = sprintf('Animal %d', i);
end

%yline(0, '--k', 'LineWidth', 1);
set(gca, 'XTick', xtick_pos, 'XTickLabel', xtick_labels, 'FontSize', 18);
ylabel('(N+ + N-)/ (N+ - N-) ', 'FontSize', 18);

xlim([0.4, nAnimals * (2 + group_gap) + 0.6]);
ylim([-1, 1]);
box off;
set(gca, 'TickDir', 'out','FontSize',18);

% Legend
legend({'MFA', 'GrC'}, 'FontSize', 15, 'Box', 'off','Location','bestoutside');
hold off;

fileName = ['PMNMratio'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');




%%
