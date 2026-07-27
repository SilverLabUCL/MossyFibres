% This script performs alignment analysis on neural activity signals and behavior (running state) 
% in multiple sessions.
clear all;
ColorCodes;

folder_names = {
    '171212_16_19_37';
    '191018_13_39_41'; 
    %'191018_13_56_55';
    '191018_14_30_00';
    %'191018_14_11_33';
    '191209_13_44_12'; 
    %'191209_14_04_14'; % only one running state
    '191209_14_32_39'; 
    '191209_15_01_22'; 
    '191209_14_18_13';
    '191209_14_46_58';
    '200130_13_21_13 FunctAcq';
    %'200130_13_36_14 FunctAcq'; % only one running state
    '200130_13_49_09 FunctAcq';
    %'200130_14_02_12 FunctAcq';
    '200130_14_15_24 FunctAcq';
    '200130_14_29_30 FunctAcq';
};

zids_PM_all = [];
zids_NM_all = [];
all_ccd = [];
all_l_max_PM = [];
all_l_min_NM = [];

savepath ='X:\MFB';
datapath = 'X:\MFB\Processed\all_data';
file_list = dir(fullfile(datapath, '*.mat'));

for ij = [1,2,6,15]
  
   current_file = fullfile(datapath, file_list(ij).name);
    
   load(current_file);
   file = mf.info.file;
    
   dff_r = mf.c.dff_r;
    dff_rz = zscore(dff_r')';
    L_state = mf.b.L_state;

    prefix = [];

    if ~all(L_state == 0)
        [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', 100, 200);
    else
        continue
    end

    zth = 2;    
    zids_PM = find(zc_wL > zth);
    PM_val = zc_wL(zids_PM);
    [sorted_PM_val, sorted_PM] = sort(PM_val, 'descend');
    zids_PM = zids_PM(sorted_PM);
    zids_PM_all = [zids_PM_all zids_PM];

    zids_NM = find(zc_wL < -zth);
    NM_val = zc_wL(zids_NM);
    [sorted_NM_val, sorted_NM] = sort(NM_val, 'ascend');
    zids_NM = zids_NM(sorted_NM);
    zids_NM_all = [zids_NM_all zids_NM];

    plot_Tuning = 1;
    plot_cross_corr_example = 1;
    plot_cross_corr_allMFs = 1;

    if plot_Tuning
        dt = 10;

        t_ext_bef = 20e3; 
        t_ext_aft = 20e3;

        run_onset = find(diff(L_state) == 1);

        index = find(run_onset > 2000, 1, 'first');

        tids1 = [run_onset(index) - round(t_ext_bef/dt) : run_onset(index) + round(t_ext_aft/dt)];    
        tids2 = [run_onset(index+1) - round(t_ext_bef/dt) : run_onset(index+1) + round(t_ext_aft/dt)];

        tt1 = (tids1 - tids1(1)) * dt / 1e3; 
        tt2 = (tids2 - tids2(1)) * dt / 1e3;

        act1 = dff_rz(:, tids1); 
        act2 = dff_rz(:, tids2);

        lw = 1;
        try
            MF_id = [zids_PM(1:7), zids_NM(1:7)];
        catch
            MF_id = [zids_PM(1:7), zids_NM(1:length(zids_NM))];
        end

        MF_no = length(MF_id);
        MF_id_sorted = flip(MF_id);

        h = figure('Position', [400 48 400 800]);
        [ha, pos] = tight_subplot(15, 1, [.005 .0002], [.025 .025], [.1 .05]);

        for cnt = 1:MF_no+1
            axes(ha(cnt));
            hold on

            if cnt == MF_no+1
                whl_rs = fastsmooth(reshape(mf.b.dis_R2', 1, []), 100, 3, 1);
                z1 = whl_rs(tids1);
                z2 = whl_rs(tids2);

                cl = 'k';
                lw = 1;

                plot(tt1, z1, 'color', 'k', 'LineWidth', lw);
                plot(tt2, z2, 'color', 'k', 'LineWidth', lw);

                zsr = L_state; 
                zsr(L_state == 0) = nan;
                plot(tt1, zsr(tids1) * .3, 'r-', 'LineWidth', 1.5);
                plot(tt2, zsr(tids2) * .23, 'r-', 'LineWidth', 1.5);

y0_temp = nanmin([z1, z2]);
yf_temp = nanmax([z1, z2]);
y_range = yf_temp - y0_temp;
x_range = max(tt1) - min(tt1);

scale_y_start = y0_temp - 0.5;
scale_y_length = 0.5;
scale_x_pos = min(tt1) + 0.08 * x_range;

plot([scale_x_pos, scale_x_pos], ...
     [scale_y_start, scale_y_start + scale_y_length], ...
     'k-', 'LineWidth', 2)
text(-3.5, 7,'0.5 MI', 'FontSize', 12, 'VerticalAlignment', 'bottom')

scale_x_start = min(tt1) + 0.08 * x_range;
scale_x_length = 5;
scale_y_pos = y0_temp +10;

plot([scale_x_start, scale_x_start + scale_x_length], ...
     [scale_y_pos, scale_y_pos], ...
     'k-', 'LineWidth', 2);
text(scale_x_start + scale_x_length/2, ...
     scale_y_pos +15, ...
     '5 s', 'FontSize', 12, 'HorizontalAlignment', 'center')
            else
                mfid = MF_id(cnt);
                z1 = act1(mfid, :);
                z2 = act2(mfid, :);
                cl = my_colors(MF_no - cnt + 1, :);
                lw = 1.5;

                plot(tt1, z1, 'color', cl, 'LineWidth', lw);
                plot(tt2, z2, 'color', cl, 'LineWidth', lw);

                if cnt == 1 
                    text(-3.5, 7, '20%', 'fontsize', 12);
                    text(-4, 4, '\DeltaF/F', 'fontsize', 12);
                end
            end

            y0 = nanmin([z1, z2]); 
            yf = nanmax([z1, z2]);
            xlim([min(tt1), max(tt1)])
            ylim([y0 - .1 * abs(y0), 1.1 * yf])
            yticklabels([]) 
            y_20_percent = 0.2 * yf;
            plot([0, 0], [y0, y_20_percent], '-', 'color', 'k', 'LineWidth', 2);

            set(gca, 'LineWidth', 1, 'FontSize', 15, ...
                'FontName', 'Helvetica', ...
                'Box', 'off', ...
                'TickDir', 'out', ...
                'TickLength', [.02 .02], 'XColor', 'none', 'YColor', 'none')
        end

        currentDate = datestr(now, 'yyyy-mm-dd');
        savepath2 = fullfile(savepath, 'Figures/Supp.Figures', 'S4', currentDate);
        if ~exist(savepath2, 'dir')
            mkdir(savepath2);
        end

        % print
        fileName = [prefix 'runAligned_example' char(file)];
        fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
        exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    end

    % - plot cross-corr
    if plot_cross_corr_example
        fs = 15;
        figure('Position', [360 78 200 800]);
        cnt = 0;
        ccd = zeros(1, length(MF_id));
        for i = 1:length(MF_id)
            if i < 5 || i > 10
                cnt = cnt + 1;
                mfid = MF_id(i);
                zz = dff_r(mfid, :);
                whl_r = mf.b.dis_R2;
                nnids = logical(~isnan(zz) .* ~isnan(whl_r));
                [r, l] = xcorr(zz(nnids), whl_r(nnids), 'normalized');

                subplot(8, 1, cnt); 
                hold on
                cl = my_colors(MF_no - i + 1, :);
                plot(l, r, 'LineWidth', 2, 'Color', cl)

                delta_t = 20e3/dt;
                lids = logical((l > -delta_t) .* (l < delta_t));

                plot([0, 0], [min(r(lids)), max(r(lids))], 'k--') 
                xlim([-delta_t, delta_t])

                xticks([-10, 0, 10] * 1e3/dt)
                xticklabels([-10, 0, 10])

                if i ~= length(MF_id_sorted)
                    xticklabels([])
                else
                    xlabel('Time lag (s)')
                    ylabel('Corr.')
                end
                ccd(i) = r(l == 0) - nanmean(r(logical((l < 0) .* (l > -delta_t))));
                
                set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [.05, .05])
            end
        end

        % print
        fileName = [prefix 'runAligned_xcorr_example' char(file)];
        fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
        exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    end

    % - for all MFs
    if plot_cross_corr_allMFs
        calculate = 1;
        if calculate
            Nmf = size(dff_r,1);
            ccd = zeros(1, Nmf);
            l_max = zeros(1, Nmf);
            l_min = zeros(1, Nmf);
            N_sh = 100;
            ccd_sh = zeros(Nmf, N_sh);
            for i = 1:Nmf
                za = dff_r(i, :);
                zb = whl_r;
                delta_t = 20e3/dt;

                [ccd(i), l_max(i), l_min(i)] = my_xcorr(za, zb, delta_t, dt);
            end

            all_ccd = [all_ccd, ccd];
            all_l_max_PM = [all_l_max_PM, l_max(zids_PM)];
            all_l_min_NM = [all_l_min_NM, l_min(zids_NM)];
        end

        fs = 15;

        figure('Position', [100, 100, 500, 200])
        subplot(121)
        hold on

        bins = -0.25 : 0.025 : 0.25;
        histogram(all_ccd, bins, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'Normalization', 'probability', 'LineWidth', 2)

        xlabel('Rel. zero-lag corr.')
        ylabel('Fraction')
        set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [.02, .02])

        subplot(122)
        hold on

        bins = 0 : 1 : 20;

        all_l_max_PM = all_l_max_PM(all_l_max_PM ~= 0 & all_l_max_PM <= 19);
        all_l_min_NM = all_l_min_NM(all_l_min_NM ~= 0 & all_l_min_NM <= 19);

        histogram(all_l_max_PM, bins, 'EdgeColor', 'r', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'Normalization', 'probability')
        histogram(all_l_min_NM, bins, 'EdgeColor', 'b', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'Normalization', 'probability')

        legend('PM', 'NM')
        legend boxoff
        xlabel('Time to peak/trough')
        ylabel('Fraction')
        set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [.02, .02])

        % print
        fileName = [prefix 'runAligned_xcorr_quant_example' char(file)];
        fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
        exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    end
end

%% Supple 4d
ColorCodes

folder_names = {
    '171212_16_19_37'
    '191018_13_39_41';
    %'191018_13_56_55';
    '191018_14_30_00';
    %'191018_14_11_33';
    '200130_13_21_13 FunctAcq';
    %'200130_13_36_14 FunctAcq';
    '200130_13_49_09 FunctAcq';
    '200130_14_15_24 FunctAcq';
    '200130_14_29_30 FunctAcq';
    '191209_13_44_12';
    '191209_14_04_14';
    '191209_14_32_39';
    '191209_15_01_22';
    '191209_14_18_13';
    '191209_14_46_58';
};

t_ext_bef = 2000;
t_ext_aft = 20000;

dt = 10;

mean_x_pm_all = [];
mean_x_nm_all = [];

peak_time_pm = [];
peak_time_nm = [];

x_pm_all = [];
x_nm_all = [];
MI_wheel_r_all = [];

for file_i = 1:length(folder_names)
    file = folder_names{file_i};
    quickAnalysis;
    
    run_onset = find(diff(L_state) == 1);
    index = find(run_onset > 2500, 1, 'first');
    
    tids1 = [run_onset(index)-round(t_ext_bef/dt) : run_onset(index)+round(t_ext_aft/dt)];    
    run1 = run_onset(index);
    
    % for first movement
    tis_before = [run_onset(index)-round(t_ext_bef/dt) : run_onset(index)];  % /100 = s
    tis_after = [run_onset(index) : run_onset(index)+round(t_ext_aft/dt)];
    
    if ~all(L_state == 0)
        [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', 100, 100);
    else
        continue
    end
    
    zth = 2;    
    
    zids = find(abs(zc_wL) > zth);
    zids_PM = find(zc_wL > zth);
    PM_val = zc_wL(zids_PM);
    [sorted_PM_val, sorted_PM] = sort(PM_val, 'descend');
    zids_PM = zids_PM(sorted_PM);
    
    zids_NM = find(zc_wL < -zth);
    NM_val = zc_wL(zids_NM);
    [sorted_NM_val, sorted_NM] = sort(NM_val, 'ascend');
    zids_NM = zids_NM(sorted_NM);

    Ntot = size(dff_rz, 1);
    
    % get PM peak latency 
    x_pm = dff_rz(zids_PM, :);
    [peak_time, onset_time, rel_inc_pm, maxima_pm, minima_pm] = peak_latency(x_pm, tis_before, tis_after);
    peak_time_pm = [peak_time_pm peak_time/dt];
    
    x1 = run1 - ((t_ext_bef)/1000 + 8.5) * 100;
    x2 = run1 + (t_ext_aft/1000 + 0.5) * 100;
    
    [ai, bi] = sort(rel_inc_pm, 'descend');
    bi(isnan(ai)) = [];
    
    for ij = 1
        xlim([x1, x2])
        xlabel('Time')
        ylabel('\DeltaF/F')
        set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01])
    end
    
    % get NM peak latency
    x_nm = dff_rz(zids_NM, :);
    [peak_time, onset_time, rel_inc_nm, maxima_nm, minima_nm] = peak_latency(x_nm, tis_before, tis_after);
    peak_time_nm = [peak_time_nm peak_time/dt];
    
    x1 = run1 - ((t_ext_bef)/1000 + 8.5) * 100;
    x2 = run1 + (t_ext_aft/1000 + 0.5) * 100;  % before 10.5s, after 20.5s
    [ai, bi] = sort(rel_inc_nm, 'descend');
    bi(isnan(ai)) = [];
    
    % Locomotion index combine
    MI_wheel_r_normalized = (MI_wheel_r(x1:x2) - nanmin(MI_wheel_r(x1:x2))) / (nanmax(MI_wheel_r(x1:x2)) - nanmin(MI_wheel_r(x1:x2)));
    MI_wheel_r_all = [MI_wheel_r_all; MI_wheel_r_normalized];
    mean_wheel_r = nanmean(MI_wheel_r_all, 1);
    
    % smooth
    x_pm_s = zeros(size(x_pm));
    for i = 1:size(x_pm, 1)
        x_pm_s(i, :) = fastsmooth(x_pm(i, :), 10, 3, 1);
    end
    
    x_pm_all = [x_pm_all; x_pm_s(:, x1:x2)];
    mean_x_pm = nanmean(x_pm_all, 1);
    
    mean_x_pm_all = [mean_x_pm_all; mean_x_pm];
    mean_mean_x_pm_all = nanmean(mean_x_pm_all, 1);
    time_axis_pm = (x1:x2) / 100 - run1 / 100;
    
    [a, b] = max(mean_mean_x_pm_all);
    time_axis_pm(b)
    
    zero_i = find(time_axis_pm == 0, 1);
    
    max_positions = nan(size(x_pm_all, 1), 1);
    
    for iii = 1:size(x_pm_all, 1)
        end_index = size(x_pm_all, 2);
        [max_value, max_index] = max(x_pm_all(iii, zero_i:end));
        max_index = max_index + zero_i - 1;
        
        [~, idx1] = min(abs(time_axis_pm - (t_ext_aft/1000 - 1)));
        if max_index >= idx1 || max_index <= find(time_axis_pm == 0)
            continue
        else
            max_positions(iii) = time_axis_pm(max_index);
        end
    end
    
    x_nm_s = zeros(size(x_nm));
    for i = 1:size(x_nm, 1)
        x_nm_s(i, :) = fastsmooth(x_nm(i, :), 10, 3, 1);
    end
    
    x_nm_all = [x_nm_all; x_nm_s(:, x1:x2)];
    mean_x_nm = mean(x_nm_all, 1);
    mean_x_nm_all = [mean_x_nm_all; mean_x_nm];
    mean_mean_x_nm_all = nanmean(mean_x_nm_all, 1);
    
    time_axis_nm = (x1:x2) / 100 - run1 / 100;
    
    zero_i = find(time_axis_nm == 0, 1);
    
    min_positions = nan(size(x_nm_all, 1), 1);
    
    for iii = 1:size(x_nm_all, 1)
        end_index = size(x_nm_all, 2); 
        [min_value, min_index] = min(x_nm_all(iii, zero_i:end));
        min_index = min_index + zero_i - 1;
        [~, idx2] = min(abs(time_axis_nm - (t_ext_aft/1000 - 1)));
        if min_index >= idx2 || min_index <= find(time_axis_nm == 0)
            continue
        else
            min_positions(iii) = time_axis_nm(min_index);
        end
    end
    
    figure('Position', [100, 200, 350, 300]);
    subplot(2, 1, 1);
    hold on;
    ylim([-1.3, 1.5]);
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
    plot(time_axis_pm, mean_mean_x_pm_all, 'r');
    xlabel('Time (s) from run onset');
    ylabel('\DeltaF/F');
    hold off;
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01]);
    
    hold on;
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
    plot(time_axis_nm, mean_mean_x_nm_all, 'b');
    xlabel('Time (s) from run onset');
    ylabel('\DeltaF/F');
    hold off;
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01]);
    
    subplot(2, 1, 2);
    hold on;
    plot(x1:x2, mean_wheel_r, 'k');
    xlabel('Time (s)');
    ylabel('Loco');
    hold off;
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01]);
    
    % print
    fileName = ['Trace_PM_NM_average.png'];
    fullFilePath = fullfile(savepath2, fileName);
    print(fullFilePath, '-dpng', '-r300');
    
    figure(Position = [100 200 200 200]);
    data = [max_positions; min_positions];
    groups = [ones(size(max_positions)); 2 * ones(size(min_positions))];
    hBoxPlot = boxplot(data, groups, 'Colors', ['r', 'b'], 'Symbol', '', 'Widths', 0.6);
    
    ylabel('Time (s)');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'PM', 'NM'});
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01]);
    
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:length(h)
        set(h(j), 'LineWidth', 1);
    end
    
    hMed = findobj(gca, 'tag', 'Median');
    for j = 1:length(hMed)
        set(hMed(j), 'LineWidth', 1);
    end
    
    [p, h, stats] = ranksum(max_positions, min_positions);
    if p < 0.001
        num_stars = '***';
    elseif p < 0.01
        num_stars = '**';
    elseif p < 0.05
        num_stars = '*';
    end
    
    if exist('num_stars', 'var')
        hold on;
        star_x = 1.5;
        star_y = max([max_positions; min_positions]) + 0.1;
        text(star_x, star_y, num_stars, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15);
        hold off;
    end
    
    fileName = ['Boxplot'];
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    
end


%% -
function [ccd, l_max, l_min] = my_xcorr(za, zb, delta_t, dt)
    nnids = logical(~isnan(za) .* ~isnan(zb));
    [r, l] = xcorr(za(nnids), zb(nnids), 'normalized');
    lids = logical((l > -delta_t) .* (l < delta_t));
    ccd = r(logical(l == 0)) - nanmean(r(logical((l < 0) .* (l > -delta_t))));

    lp = l(l > 0);
    [~, zi] = max(r(logical((l > 0) .* (l < delta_t))));
    l_max = lp(zi) * dt / 1e3;
    [~, zi] = min(r(logical((l > 0) .* (l < delta_t))));
    l_min = lp(zi) * dt / 1e3;
end

function [x, y] = my_hist(zz, bins)
    [y, x] = histcounts(zz, bins, 'Normalization', 'count');
    dx = diff(x(1:2));
    xx = x(1:end-1) + dx/2;
    N_tot = length(zz);
    plot(xx, y/N_tot*100, 'LineWidth', 2, 'color', 'k')
end
