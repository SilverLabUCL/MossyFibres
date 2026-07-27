%
% This script displays sample z-scored dF/F traces, and generates 
% bootstrapped correlation, dertermines PM/NM/NSM
%
% Author: Yizhou
% Date: 2023/11/27

file = '171212_16_19_37'; 
quickAnalysis

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure2', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

%% plot CC
cc_MF_all = corr(dff_rz', 'rows', 'complete');
cc_MF_stat = corr(dff_rz(:, L_state==0)', 'rows', 'complete');
cc_MF_run = corr(dff_rz(:, L_state==1)', 'rows', 'complete');

for i = 1:3
    if i == 1
        cc_MF = cc_MF_all;
        suffix = 'all';
        ttl = 'All';
    elseif i == 2
        cc_MF = cc_MF_stat;
        suffix = 'stat';
        ttl = 'QW';
    elseif i == 3
        cc_MF = cc_MF_run;
        suffix = 'run';
        ttl = 'A-state';
    end
    
    % reorder by hierarchical clustering
    D = 1 - cc_MF;
    D(1:size(D,1)+1:end) = 0;
    Dv = squareform(D);
    Z = linkage(Dv, 'average');
    ord = optimalleaforder(Z, Dv);
    cc_MF = cc_MF(ord, ord);
    
    figure('Position', [100, 100, 300, 300])
    subplot(111)
    hold on
    title(ttl, 'FontWeight', 'normal');
    imagesc(cc_MF, [-1, 1])
    cb = colorbar();
    colormap('redblue')
    set(cb, 'position', [0.1, 0.2, 0.03, 0.15], 'Ticks', [-1, 0, 1])
    axis image
    xticks([1, Nmf])
    yticks([1, Nmf])
    xlabel('MFA#')
    ylabel('MFA#')
    set(gca, 'LineWidth', 0.1, 'FontSize', 17, 'TickDir', 'out')
    
    % print
    fileName = ['cc_matrix_' suffix];
    fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
end

%% - cc distribution
cc_MF_all(eye(length(cc_MF_all)) == 1) = nan;
cc_MF_stat(eye(length(cc_MF_stat)) == 1) = nan;
cc_MF_run(eye(length(cc_MF_run)) == 1) = nan;

figure('Position', [100, 100, 300, 250])
subplot(111)
hold on

[y, x] = histcounts(cc_MF_all, 'Normalization', 'count');
dx = diff(x(1:2));
xx = x(1:end-1) + dx/2;
N_tot = Nmf^2 - Nmf;
plot(xx, y / N_tot * 100, 'LineWidth', 2, 'color', [0, 0, 0])

set(gca, 'LineWidth', 1, 'FontSize', 20, 'TickDir', 'out')
xlabel('Pairwise corr.')
ylabel('% Pairs')
xlim([-1, 1]);

% print
fileName = ['cc_dist_allMF.png'];
fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

figure('Position', [100, 100, 300, 250])
subplot(111)
hold on

[y, x] = histcounts(cc_MF_stat, 'Normalization', 'probability');
dx = diff(x(1:2));
xx = x(1:end-1) + dx/2;
N_tot = Nmf^2 - Nmf;
plot(xx, y, 'LineWidth', 2, 'color', [50, 200, 200] / 255)

[y, x] = histcounts(cc_MF_run, 'Normalization', 'probability');
dx = diff(x(1:2));
xx = x(1:end-1) + dx/2;
N_tot = Nmf^2 - Nmf;
plot(xx, y, 'LineWidth', 2, 'color', [200, 50, 200] / 255)

legend('QW', 'AS', 'Location', 'best')
legend boxoff

xlabel('Pairwise corr.')
ylabel('Prob.')
set(gca, 'LineWidth', 1, 'FontSize', 17, 'TickDir', 'out')

% print
fileName = ['cc_dist_runVstat_' file];
fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% -- 2a. Heat Map of Activity
ccs = corr(spd_r', dff_r', 'rows', 'complete');
[a, b] = sort(ccs);

plot_heatMap = 1;

r2w = [linspace(1, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)'];
w2b = [linspace(1, 0, 128)', linspace(1, 0, 128)', linspace(1, 1, 128)'];
rbMap = [r2w; w2b];
mycm = flipud(rbMap);

if plot_heatMap
    suffix = 'changeSorted';
    my_b = b;
    Ntot = Nmf;
    dt = 10;
    
    xlm = [5000, 130 * 1e3 / dt];
    
    figure('Position', [360, 20, 480, 700]);
    sbplt1 = 1:8;
    sbplt2 = 9;
    sbplt3 = 10;
    sbp_no = 10;
    
    zz = dff_rz(my_b, :);
    
    subplot(sbp_no, 1, sbplt1); 
    hold on;
    
    sorted_zz = sort(zz(~isnan(zz)), 'ascend');
    n = numel(sorted_zz);
    v1 = round(sorted_zz(max(floor(n * 0.005), 1)), 1);
    v2 = round(sorted_zz(max(floor(n * 0.995), 1)), 1);
    
    imagesc(zz, [v1, v2]);
    ylb = ylabel('Mossy fibre axon number');
    set(ylb, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    yticks([1, Ntot]);
    ylim([0.5, Ntot + 0.5]);
    ylim([0, size(zz, 1) + 4]);
    
    for j = 1:Ntr
        tpuff = ((j - 1) * T(end) + 5e3) / dt;
        yTriTop = size(zz, 1) + 3;
        plot(tpuff, yTriTop, 'v', 'MarkerSize', 4, 'MarkerEdgeColor', [0.7, 0.7, 0.7], 'MarkerFaceColor', [0.7, 0.7, 0.7]);
    end
    
    xticklabels([]);
    xlim(xlm);
    cb = colorbar();
    caxis([-2, 5]);
    set(cb, 'position', [0.07, 0.82, 0.02, 0.05], 'Ticks', [v1, v2], 'FontSize', 11);
    title(cb, {'Z-scored', '\DeltaF/F'}, 'FontSize', 10, 'FontWeight', 'normal');
    
    pos = get(gca, 'Position');
    pos(1) = 0.15; pos(3) = 0.75;
    pos(2) = 0.26; pos(4) = 0.7;
    set(gca, 'Position', pos, 'LineWidth', 1, 'FontSize', 17, 'XColor', [0.6, 0.6, 0.6], 'YDir', 'normal', 'TickDir', 'out');
    
    subplot(sbp_no, 1, sbplt3); 
    hold on;
    
    AS_regions = (L_state == 1);
    AS_start = find(diff([0; AS_regions(:)]) == 1);
    AS_end = find(diff([AS_regions(:); 0]) == -1);
    time = 1:length(L_state);
    
    for i = 1:length(AS_start)
        patch([time(AS_start(i)), time(AS_end(i)), time(AS_end(i)), time(AS_start(i))], ...
              [-0.05, -0.05, 0.35, 0.35], ...
              [255, 231, 255] / 255, 'EdgeColor', 'none');
    end
    
    QW_regions = (Q_state == 1);
    QW_start = find(diff([0; QW_regions(:)]) == 1);
    QW_end = find(diff([QW_regions(:); 0]) == -1);
    
    for i = 1:length(QW_start)
        patch([time(QW_start(i)), time(QW_end(i)), time(QW_end(i)), time(QW_start(i))], ...
              [-0.05, -0.05, 0.35, 0.35], ...
              [192, 255, 255] / 255, 'EdgeColor', 'none');
    end
    
    wsk_r = MI_whisker_r;
    plot(time, wsk_r, 'LineWidth', 1, 'Color', [255, 163, 26] / 255);
    ylabel({'Whisk'});
    set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    yticks([0, 0.3]);
    
    pos = get(gca, 'Position');
    pos(1) = 0.15; pos(3) = 0.75;
    pos(2) = 0.04; pos(4) = 0.075;
    xlim(xlm);
    set(gca, 'Position', pos, 'LineWidth', 1, 'FontSize', 15, 'XColor', 'none', 'TickDir', 'out');
    
    x0 = max(xlim) - 250;
    y0 = min(wsk_r);
    plot(x0 - [0, 2000], [y0, y0], 'k-', 'LineWidth', 3);
    text(x0 - 500, y0 - 0.13, '20 s', 'FontSize', 15, 'HorizontalAlignment', 'center');
    
    subplot(sbp_no, 1, sbplt2); 
    hold on;
    for i = 1:length(AS_start)
        patch([time(AS_start(i)), time(AS_end(i)), time(AS_end(i)), time(AS_start(i))], ...
              [-0.05, -0.05, 1, 1], ...
              [255, 231, 255] / 255, 'EdgeColor', 'none');
    end
    
    for i = 1:length(QW_start)
        patch([time(QW_start(i)), time(QW_end(i)), time(QW_end(i)), time(QW_start(i))], ...
              [-0.05, -0.05, 1, 1], ...
              [192, 255, 255] / 255, 'EdgeColor', 'none');
    end
    
    plot(time, rescale(dis_R2), 'Color', [173, 210, 157] / 255, 'LineWidth', 1);
    ylabel({'Loco'});
    yticks([0, 0.8]);
    set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    
    pos = get(gca, 'Position');
    pos(1) = 0.15; pos(3) = 0.75;
    pos(2) = 0.16; pos(4) = 0.075;
    xlim(xlm);
    set(gca, 'Position', pos, 'LineWidth', 1, 'FontSize', 15, 'XColor', 'none', 'TickDir', 'out');
    
    fileName = ['HeatMap__', suffix, '.png'];
    fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
end

%% 2b. traces
plot_sampledff_separate = 1;
cl_nm = mycm(40, :);
cl_pm = mycm(256 - 40, :);
cl_ns = [0.7, 0.7, 0.7];
MF_id_down = my_b(1:3);
MF_id_up = my_b(end-2:end);
MF_id_ns = my_b(median(my_b) : median(my_b) + 2);
states = L_state;

if plot_sampledff_separate
    xlm = [0, 380 * 1e3 / dt];
    
    for ijk = 1:3
        if ijk == 1
            MF_id_sel = MF_id_down;
            cl = cl_nm;
            suffix = 'down';
        elseif ijk == 2
            MF_id_sel = MF_id_up;
            cl = cl_pm;
            suffix = 'up';
        elseif ijk == 3
            MF_id_sel = MF_id_ns;
            cl = cl_ns;
            suffix = 'ns';
        end
        
        MF_no_sel = length(MF_id_sel);
        h = figure('Position', [484, 119, 750, 600], 'visible', 'on');
        pre = 30;
    
        for ii = 1:MF_no_sel
            hsp = subplot(MF_no_sel, 1, ii); 
            hold on;
            if ii <= MF_no_sel
                mfid = MF_id_sel(ii);
                z2 = dff_rz(mfid, pre:end-pre);
                ymn = nanmin(z2);
                ymx = nanmean(z2) + 2 * nanstd(z2);
    
                h_group_AS = hggroup;
                for i = 1:length(AS_start)
                    patch([time(AS_start(i)), time(AS_end(i)), time(AS_end(i)), time(AS_start(i))], ...
                          [ymn - 3, ymn - 3, ymx + 3, ymx + 3], ...
                          [255, 231, 255] / 255, 'EdgeColor', 'none', 'Parent', h_group_AS);
                end
        
                h_group_QW = hggroup;
                QW_regions = Q_state == 1;
                QW_start = find(diff([0; QW_regions(:)]) == 1);
                QW_end = find(diff([QW_regions(:); 0]) == -1);
                for i = 1:length(QW_start)
                    patch([time(QW_start(i)), time(QW_end(i)), time(QW_end(i)), time(QW_start(i))], ...
                          [ymn - 3, ymn - 3, ymx + 3, ymx + 3], ...
                          [192, 255, 255] / 255, 'EdgeColor', 'none', 'Parent', h_group_QW);
                end
    
                trace = plot(fastsmooth(z2, 500, 3, 1), 'color', cl, 'LineWidth', 4);
                xlim([0, 380 * 1e3 / dt]);
                ylim([0.9 * ymn, 1.6 * ymx])
                pos = get(hsp, 'Position');
                pos(2) = pos(2) + 0.014 * (MF_no_sel - ii);
                set(hsp, 'Position', pos);
            end
            
            set(gca, 'LineWidth', 1, 'FontSize', 18, ...
                'FontName', 'Helvetica', ...
                'Box', 'off', ...
                'YTickLabel', [], ...
                'Ytick', [], ...
                'XMinorTick', 'off', 'YMinorTick', 'off', 'XColor', 'none')
        end
    
        if ijk == 3
            x0 = max(xlim) - 5000;
            y0 = min(z2) + 1;
            plot(x0 - [0, 2000], [y0, y0], 'k-', 'LineWidth', 3);
            text(x0 + 10, y0 - 0.05 * range(z2), '20 s', 'FontSize', 21, 'HorizontalAlignment', 'center');
            plot([x0, x0], [y0, y0 + 1], 'k-', 'LineWidth', 3)
            text(x0 + 400, y0 + 1.5, '1 \DeltaF/F', 'FontSize', 21)
            set(h_group_AS, 'visible', 'off')
            set(h_group_QW, 'visible', 'off')
            set(trace, 'visible', 'off')
        end
        
        % print
        fileName = ['sample_traces_separate__' suffix '.png'];
        fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
        exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    end
end

run_r_all = run_r;

filePath = fullfile('X:\MFB\MFB_AH_2023\Correlation_data', ['corr_', folder_name, '.mat']);
save(filePath, 'cc_MF_all', 'cc_MF_stat', 'cc_MF_run', 'folder_name', 'dff_r', ...
    'L_state', 'A_state', 'Q_state', 'MI_whisker_r', 'dff_rz', 'MI_wheel_r', 'spd_r', 'run_r_all', 'Nmf', 'whl_r', 'whl_rs', 'xyz');

%% pie
figure('Position', [100, 100, 320, 270])
subplot(111); 
hold on
axis off
block_size = 100;
sig_level = 2;
[cc_0, zc_0] = bootstrap_cc(dff_rz', L_state', block_size, 200);
pm_counts = nansum(zc_0 > sig_level);
ns_counts = nansum(abs(zc_0) < sig_level);
nm_counts = nansum(zc_0 < -sig_level);
h = pie([nm_counts, ns_counts, pm_counts], {'MFA⁻', 'MFA*', 'MFA⁺'});
colormap([cl_nm; cl_ns; cl_pm]);
set(findobj(h, 'type', 'text'), 'fontsize', 20);
axis image;
set(gca, 'LineWidth', 1, 'FontSize', 24, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.01, 0.01])
ht = title({'Modulation (AS)'}, 'fontsize', 20, 'fontweight', 'normal');
titlePos = get(ht, 'Position');
newTitlePos = titlePos + [0, 0.3, 0];
set(ht, 'Position', newTitlePos);
set(gca, 'Position', [0.13, 0.15, 0.6, 0.6])
fileName = ['modLoc_' file '__pie.png'];
fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

h = figure('Position', [100, 100, 300, 250]); 
hold on;
sig_level = 2;
db = 1;
zz = zc_0(zc_0 > sig_level);
bins = [min(zz(:)) - db/2 : db : max(zz(:)) + db/2];
histogram(zz, bins, 'FaceColor', cl_pm, 'FaceAlpha', 0.85)
zz = zc_0(zc_0 < -sig_level);
bins = [min(zz(:)) - db/2 : db : max(zz(:)) + db/2];
histogram(zz, bins, 'FaceColor', cl_nm, 'FaceAlpha', 0.85)
zz = zc_0(abs(zc_0) < sig_level);
bins = [min(zz(:)) - db/2 : db : max(zz(:)) + db/2];
histogram(zz, bins, 'FaceColor', cl_ns, 'FaceAlpha', 0.85)
xlabel('Corr. with AS')
ylabel('# MFAs')
set(gca, 'LineWidth', 1, 'FontSize', 17)
% print
fileName = ['Bootstrapped Corr_' file];
fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% second half
% 2e
corr_all_combined = 1;
% fig.2d
if corr_all_combined == 1
    load('X:\MFB\MFB_AH_2023\Correlation_data\CorrLoco.mat')
    
    figure('Position', [100, 100, 400, 350])
    subplot(111)
    hold on
    [y, x] = histcounts(all_cc_MF_stat, 'Normalization', 'probability');
    dx = diff(x(1:2));
    xx = x(1:end-1) + dx/2;
    plot(xx, y, 'LineWidth', 2, 'Color', [50, 200, 200] / 255);
    hold on;
    [y, x] = histcounts(all_cc_MF_run, 'Normalization', 'probability');
    dx = diff(x(1:2));
    xx = x(1:end-1) + dx/2;
    plot(xx, y, 'LineWidth', 2, 'Color', [200, 50, 200] / 255);
    
    legend('QW', 'AS', 'Location', 'best')
    legend boxoff
    
    xlabel('Pairwise corr.')
    ylabel('Probability')
    set(gca, 'LineWidth', 1, 'FontSize', 22, 'TickDir', 'out', 'box', 'off')
    
    mfbFolderPath = 'X:\MFB';
    currentDate = datestr(now, 'yyyy-mm-dd');
    folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure2', currentDate);
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    
    fileName = ['Pairwise corr_probability'];
    fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    
    % -
    figure('Position', [100, 100, 280, 230])
    subplot(111)
    hold on
    [y, x] = histcounts(all_cc_MF_stat, 'Normalization', 'cdf');
    dx = diff(x(1:2));
    xx = x(1:end-1) + dx/2;
    plot(xx, y, 'LineWidth', 3, 'Color', [50, 200, 200] / 255);
    hold on;
    [y, x] = histcounts(all_cc_MF_run, 'Normalization', 'cdf');
    dx = diff(x(1:2));
    xx = x(1:end-1) + dx/2;
    plot(xx, y, 'LineWidth', 3, 'Color', [200, 50, 200] / 255);
    
    xlabel('Pairwise corr.')
    ylabel('Cum prob')
    set(gca, 'LineWidth', 1, 'FontSize', 22, 'TickDir', 'out')
    
    [H, P, STAT] = kstest2(all_cc_MF_stat, all_cc_MF_run);
    
    % print
    fileName = ['Pairwise corr_cdf'];
    fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
end

%% - 2f('Bootstrapped Corr.')
filePath = fullfile('X:\MFB\MFB_AH_2023\Correlation_data', ['CorrLoco.mat']);
load(filePath); 
block_size = 100;

h = figure('Position', [100, 100, 300, 250]);
hold on;
sig_level = 2;
db = 1;
zz = zc_wl(zc_wl > sig_level);
bins = [min(zz(:)) - db/2 : db : max(zz(:)) + db/2];
histogram(zz, bins, 'FaceColor', cl_pm, 'FaceAlpha', 0.85)
zz = zc_wl(zc_wl < -sig_level);
bins = [min(zz(:)) - db/2 : db : max(zz(:)) + db/2];
histogram(zz, bins, 'FaceColor', cl_nm, 'FaceAlpha', 0.85)
zz = zc_wl(abs(zc_wl) < sig_level);
bins = [min(zz(:)) - db/2 : db : max(zz(:)) + db/2];
histogram(zz, bins, 'FaceColor', cl_ns, 'FaceAlpha', 0.85)
xlabel('Corr. with QW/AS')
ylabel('# MFAs')
set(gca, 'LineWidth', 1, 'FontSize', 15)
% print
fileName = ['Bootstrapped Corr'];
fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig.2g pie
figure('Position', [100, 100, 500, 450])
subplot(111);
hold on
axis off
pm_counts = nansum(zc_wl > sig_level);
ns_counts = nansum(abs(zc_wl) < sig_level);
nm_counts = nansum(zc_wl < -sig_level);
total_counts = pm_counts + ns_counts + nm_counts;
h = pie([nm_counts, ns_counts, pm_counts], {'NM', 'NSM', 'PM'});
colormap([cl_nm; cl_ns; cl_pm]);
set(findobj(h, 'Type', 'text'), 'fontsize', 24);
percentages = [nm_counts, ns_counts, pm_counts] / total_counts * 100;
textObjects = findobj(h, 'Type', 'text');

for i = 1:length(textObjects)
    Txt1 = textObjects(i).String;
    Txt2 = sprintf('%.1f%%', percentages(i));
    set(textObjects(i), 'String', {Txt1; Txt2});
    pos = get(textObjects(i), 'Position');
    pos(2) = pos(2) + 0.05; 
    set(textObjects(i), 'Position', pos);
end

axis image;

ht = title({'Modulation (QW/AS)'}, 'fontsize', 26, 'fontweight', 'normal');
titlePos = get(ht, 'Position');
newTitlePos = titlePos + [0, 0.3, 0];
set(ht, 'Position', newTitlePos);
set(gca, 'Position', [0.13, 0.15, 0.6, 0.6])

fileName = ['modLoc_allMice__pie'];
fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
