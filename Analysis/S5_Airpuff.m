%   This script performs analysis of peak latency in response 
%   to a puff stimulus across multiple experimental datasets.
%
%   Note:
%   - Parameters such as "tpuff", "t_before", and "t_after" are set based on experimental conditions.

% datasets with airpuff
folder_names = {
    '171212_16_19_37'; % 1
    '191018_13_39_41';
    '191018_13_56_55';
    '200130_13_21_13 FunctAcq'; %4
    '200130_14_29_30 FunctAcq';
    '191209_13_44_12';
    '191209_14_18_13';
    '191209_14_46_58';
};


%% S5.b-e
low_latency_t2p_all = [];
savepath = 'X:/MFB';
n1_all=0;
n2_all=0;
n3_all=0;

base_path = 'X:\MFB\Concatenation air puff\air puff';
file_list = dir(fullfile(base_path, 'concatenated_air_puff_*.mat'));

for file_i = 1:length(file_list)
    full_name = file_list(file_i).name;
    
    tokens = regexp(full_name, 'concatenated_air_puff_(.*)\.mat', 'tokens');
    
    if ~isempty(tokens)
        file = tokens{1}{1};
        fprintf('Processing file: %s -> Extracted identity: %s\n', full_name, file);
    else
        file = erase(full_name, 'concatenated_air_puff_');
        file = erase(file, '.mat');
    end
    
    load(fullfile(base_path, full_name));
    
    speed = all_speed;
    dff_gz = all_dff_gz;
    
    ColorCodes
    tpuff = 5e3;
    after_puff = 7; 
    
    Tb1 = 4e3;
    Tb2 = tpuff;
    Ta1 = tpuff;
    Ta2 = after_puff * 1e3;
    
    time_wind_before = Tb2 - Tb1;
    time_wind_after = Ta2 - Ta1;
    
    tis_before = find((T >= Tb1) .* (T < Tb2));
    tis_after = find((T >= Ta1) .* (T <= Ta2));
    
    tt = T(tis_after) / 1e3;
    
    spd_th = 0.05;
    spdm = nanmean(speed(:, [tis_before(1):500, 511:tis_after(end)]), 2);
    
    ROI_dff_mfa = dff_gz(:, :, :);
    all_dff_m = squeeze(nanmean(ROI_dff_mfa(:, abs(spdm) < spd_th, :), 2));
    disp(['fraction filtered by speed: ', num2str(sum(abs(spdm(:)) < spd_th) / numel(spdm) * 100), '%'])
    
    all_dff_m_tp = all_dff_m';
    zz = zeros(size(all_dff_m_tp));
    mu = nanmean(all_dff_m_tp, 1);
    sigma = nanstd(all_dff_m_tp, 0, 1);
    
    for i = 1:size(all_dff_m_tp, 2)
        zz(:, i) = (all_dff_m_tp(:, i) - mu(i)) / sigma(i);
    end
    zz = zz';
    
    [peak_time_temp, ~, ~, ~, ~] = peak_latency(zz, tis_before, tis_after);
    clean_peaks = isnan(peak_time_temp) | (peak_time_temp > 80);
    disp(['fraction kept after peak time filter: ', num2str(sum(clean_peaks) / length(clean_peaks) * 100), '%'])
   
    zz = zz(clean_peaks, :);
    Ntot = size(zz, 1);
    
    [peak_time, onset_time, rel_inc, maxima, minima] = peak_latency(zz, tis_before, tis_after);
    
    sig_inc_idx = find(~isnan(maxima));
    sig_dec_idx = find(~isnan(minima));
    all_idx = 1:Ntot;
    non_sig_idx = setdiff(all_idx, [sig_inc_idx, sig_dec_idx]);
    
    n1 = numel(sig_inc_idx);
    n2 = numel(sig_dec_idx);
    n3 = numel(non_sig_idx);
    
    n1_all = n1_all + n1;
    n2_all = n2_all + n2;
    n3_all = n3_all + n3;
    
    currentDate = datestr(now, 'yyyy-mm-dd');
    savepath2 = fullfile(savepath, 'Figures', 'Supp.Figures', 'S5', currentDate);
    if ~exist(savepath2, 'dir')
        mkdir(savepath2);
    end
    
    % ====================================================================
    % ====================================================================
    
    sorting_configs = {
        struct('ref_time', onset_time, 'fig_suffix', 'onsetSortedSplit_', 'line_color', [0, 0.8, 0.2], 'title_str', 'Onset-sorted')
        struct('ref_time', peak_time,  'fig_suffix', 'peakSortedSplit_',  'line_color', [1, 0.8, 0.1], 'title_str', 'Peak-sorted')
    };
    
    for ijk = 1:2
        config = sorting_configs{ijk};
        ref_time = config.ref_time;
        
        [~, sort_inc] = sort(ref_time(sig_inc_idx));
        [~, sort_dec] = sort(ref_time(sig_dec_idx));
        [~, sort_none] = sort(ref_time(non_sig_idx));
        
        inc_sorted = sig_inc_idx(sort_inc);
        dec_sorted = sig_dec_idx(sort_dec);
        none_sorted = non_sig_idx(sort_none);
        
        sorted_idx = [inc_sorted, dec_sorted, none_sorted];
        
        fig = figure('Position', [360 20 369 600]); 
        hold on;
        imagesc(zz(sorted_idx, :), [0 5]);
        colormap(gray);
        
        dt = 10;
        
        plot([tpuff/dt, tpuff/dt], [0, Ntot], '--', 'LineWidth', 2, 'color', 'w');
        
        jitter_range = 0.25;
        time_vec = (tpuff + ref_time(sorted_idx)) / dt;
        x_jitter = (rand(1, numel(time_vec)) - 0.5) * 2 * jitter_range;
        time_jittered = time_vec + x_jitter;
        
        plot(time_jittered, 1:Ntot, '-', 'Color', config.line_color, 'LineWidth', 2);
        
        y1 = 0.5;
        y2 = y1 + n1;
        y3 = y2 + n2;
        
        h1 = n1 - 0.5;
        h2 = n2 - 1;
        h3 = n3 - 0.5;
        
        box_width = diff([1, 7000]) / dt;
        
        rectangle('Position', [3e3/dt, y1, box_width, h1], 'EdgeColor', [1 0 0], 'LineWidth', 2);
        rectangle('Position', [3e3/dt, y2, box_width, h2], 'EdgeColor', [0 0 1], 'LineWidth', 2);
        rectangle('Position', [3e3/dt, y3, box_width, h3], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
        
        xticks((3:10) * 1e3 / dt);
        xticklabels({'-2'; '-1'; '0'; '1'; '2'; '3'; '4'; '5'});
        xlabel('Time relative to puff (s)');
        ylabel('Sorted MFA #');
        xlim([3e3, 10e3] / dt);
        ylim([.5, Ntot + .5]);
        
        title(config.title_str);
        
        set(gca, 'LineWidth', 1, 'FontSize', 18, ...
            'FontName', 'Helvetica', ...
            'Box', 'off', ...
            'TickDir', 'out', ...
            'TickLength', [.02 .02], ...
            'XMinorTick', 'off', ...
            'YMinorTick', 'off', ...
            'YDir', 'reverse');
        
        fileName = [config.fig_suffix 'HeatMap_dff_' num2str(after_puff) 's_' char(file)];
        fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
        exportgraphics(fig, fullFilePathPDF, 'ContentType', 'vector');
    end
    
    % ====================================================================
    % ====================================================================
    
    t_before = 0.5e3;
    t_after = 2e3;
    dt = 10;
    
    x1 = t_before / dt;
    x2 = (tpuff + t_after) / dt;
    
    low_latency = 600;
    mid_latency = 1200;
    
    colors = { [0.7, 0.7, 0.7], ...
               [0.7, 0.7, 0.7], ...
               [0.7, 0.7, 0.7] };
    
    fig_summary = figure('Position', [200 100 350 600]);
    
    ymin = -6; 
    ymax = 8;
    
    subplot(3,1,1); hold on;
    title(sprintf('MFA+ (n = %d)', numel(sig_inc_idx)));
    
    for ij = 1:numel(sig_inc_idx)
        ii = sig_inc_idx(ij);
        if peak_time(ii) <= low_latency
            color = colors{1};
        elseif peak_time(ii) <= mid_latency
            color = colors{2};
        else
            color = colors{3};
        end
        plot(T/1000, zz(ii, :), 'Color', color, 'LineWidth', 1);
    end
    
    plot([tpuff, tpuff]/1000, [ymin, ymax], 'k--', 'LineWidth', 0.6);
    plot([tpuff + t_after, tpuff + t_after]/1000, [ymin, ymax], 'k-', 'LineWidth', 1.2);
    
    avg_inc = nanmean(zz(sig_inc_idx, :), 1);
    plot(T/1000, avg_inc, 'k-', 'LineWidth', 2);
    
    xlim([T(x1), T(x2)*1.5]/1000);
    ylim([ymin, ymax]);
    set(gca, 'XColor', 'none', 'Box', 'off', 'FontSize', 18);
    ylabel('\DeltaF/F');
    
    subplot(3,1,2); hold on;
    title(sprintf('MFA- (n = %d)', numel(sig_dec_idx)));
    
    for ij = 1:numel(sig_dec_idx)
        ii = sig_dec_idx(ij);
        if peak_time(ii) <= low_latency
            color = colors{1};
        elseif peak_time(ii) <= mid_latency
            color = colors{2};
        else
            color = colors{3};
        end
        plot(T/1000, zz(ii, :), 'Color', color, 'LineWidth', 1);
    end
    
    plot([tpuff, tpuff]/1000, [ymin, ymax], 'k--', 'LineWidth', 0.6);
    plot([tpuff + t_after, tpuff + t_after]/1000, [ymin, ymax], 'k-', 'LineWidth', 1.2);
    
    avg_dec = nanmean(zz(sig_dec_idx, :), 1);
    plot(T/1000, avg_dec, 'k-', 'LineWidth', 2);
    
    xlim([T(x1), T(x2)*1.5]/1000);
    ylim([ymin, ymax]);
    set(gca, 'XColor', 'none', 'Box', 'off', 'FontSize', 18);
    ylabel('\DeltaF/F');
    
    subplot(3,1,3); hold on;
    title('Speed (across trials)');
    
    speed_crop = speed(abs(spdm) < spd_th, :);
    avg_speed = mean(speed_crop, 1, 'omitnan');
    
    plot(T/1000, avg_speed, 'Color', [0.3 0.3 0.9], 'LineWidth', 1.5);
    
    plot([tpuff, tpuff]/1000, [-5 5], 'k--', 'LineWidth', 0.6);
    plot([tpuff + t_after, tpuff + t_after]/1000, [-5 5], 'k-', 'LineWidth', 1.2);
    
    xlim([T(x1), T(x2)*1.5]/1000);
    ylim([-5 5]);
    xlabel('Time (s)');
    ylabel('Speed (cm/s)');
    set(gca, 'Box', 'off', 'FontSize', 18);
    
    fileName = ['summary_PM_NM_Speed_' char(file)];
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(fig_summary, fullFilePathPDF, 'ContentType', 'vector')
    
    low_latency_idx = peak_time <= low_latency;
    low_latency_t2p = peak_time(low_latency_idx);
    low_latency_t2p_all = [low_latency_t2p_all low_latency_t2p];
end

fprintf('\n=== Total ===\n');
fprintf('Total PM (sig_inc): %d\n', n1_all);
fprintf('Total NM (sig_dec): %d\n', n2_all);
fprintf('Total Non-sig: %d\n', n3_all);
fprintf('Total low-latency: %d\n', numel(low_latency_t2p_all));


%% Supp 5f
time_to_peak = struct( ...
    'Animal1', struct('inc', [], 'dec', []), ...
    'Animal2', struct('inc', [], 'dec', []), ...
    'Animal3', struct('inc', [], 'dec', []), ...
    'Animal4', struct('inc', [], 'dec', []));

Animal1 = {'171212_16_19_37'};
Animal4 = {'200130_13_21_13 FunctAcq', '200130_14_29_30 FunctAcq'};
Animal2 = {'191018_13_39_41', '191018_13_56_55'};
Animal3 = {'191209_13_44_12',  '191209_14_18_13', '191209_14_46_58'};


base_path = 'X:\MFB\Concatenation data\air puff';
file_list = dir(fullfile(base_path, 'concatenated_air_puff_*.mat'));

for file_i = 1:length(file_list)
    full_name = file_list(file_i).name;
    
    tokens = regexp(full_name, 'concatenated_air_puff_(.*)\.mat', 'tokens');
    
    if ~isempty(tokens)
        file = tokens{1}{1};
        fprintf('Processing file: %s -> Extracted identity: %s\n', full_name, file);
    else
        file = erase(full_name, 'concatenated_air_puff_');
        file = erase(file, '.mat');
    end
    
    load(fullfile(base_path, full_name));

    speed = all_speed;
    dff_gz = all_dff_gz;
    
    tpuff = 5e3;
    t_before = 0.5 * 1e3;
    t_after = 2 * 1e3;
    
    Tb1 = tpuff - t_before;
    Tb2 = tpuff;
    Ta1 = tpuff;
    Ta2 = tpuff + t_after;
    
    tis_before = find((T >= Tb1) .* (T < Tb2));
    tis_after = find((T >= Ta1) .* (T <= Ta2));
    

    spd_th = 0.05;
    spdm = nanmean(speed(:, [tis_before(1):500, 511:tis_after(end)]), 2); % mean speed during the window
    
    ROI_dff_mfa = dff_gz(:, :, :);
    ['fraction filtered: ', num2str(sum(spdm(:) < spd_th) / numel(spdm) * 100), '%']
    x = squeeze(nanmean(ROI_dff_mfa(:, abs(spdm) < spd_th, :), 2));
    [peak_time, ~, ~, ~, ~,~] = peak_latency(x, tis_before, tis_after);
    clean_peaks = isnan(peak_time) | (peak_time > 300);
    x= x(clean_peaks,:);
     [peak_time, onset_time, rel_inc, maxima, minima,peak_type] = peak_latency(x, tis_before, tis_after);


    idx_inc = peak_type == 1;
    idx_dec = peak_type == -1;

    
    pt_inc = peak_time(idx_inc);
    pt_dec = peak_time(idx_dec);
    
    if strcmp(file, 'animal 1')
        time_to_peak.Animal1.inc = [time_to_peak.Animal1.inc, pt_inc];
        time_to_peak.Animal1.dec = [time_to_peak.Animal1.dec, pt_dec];
    elseif strcmp(file, 'animal 2')
        time_to_peak.Animal2.inc = [time_to_peak.Animal2.inc, pt_inc];
        time_to_peak.Animal2.dec = [time_to_peak.Animal2.dec, pt_dec];
    elseif strcmp(file, 'animal 3')
        time_to_peak.Animal3.inc = [time_to_peak.Animal3.inc, pt_inc];
        time_to_peak.Animal3.dec = [time_to_peak.Animal3.dec, pt_dec];
    elseif strcmp(file, 'animal 4')
        time_to_peak.Animal4.inc = [time_to_peak.Animal4.inc, pt_inc];
        time_to_peak.Animal4.dec = [time_to_peak.Animal4.dec, pt_dec];
    end

end

% mean & std for sig_inc
mean_peak_time_Animal1_inc = nanmean(time_to_peak.Animal1.inc);
mean_peak_time_Animal2_inc = nanmean(time_to_peak.Animal2.inc);
mean_peak_time_Animal3_inc = nanmean(time_to_peak.Animal3.inc);
mean_peak_time_Animal4_inc = nanmean(time_to_peak.Animal4.inc);

std_peak_time_Animal1_inc = nanstd(time_to_peak.Animal1.inc);
std_peak_time_Animal2_inc = nanstd(time_to_peak.Animal2.inc);
std_peak_time_Animal3_inc = nanstd(time_to_peak.Animal3.inc);
std_peak_time_Animal4_inc = nanstd(time_to_peak.Animal4.inc);

% mean & std for sig_dec
mean_peak_time_Animal1_dec = nanmean(time_to_peak.Animal1.dec);
mean_peak_time_Animal2_dec = nanmean(time_to_peak.Animal2.dec);
mean_peak_time_Animal3_dec = nanmean(time_to_peak.Animal3.dec);
mean_peak_time_Animal4_dec = nanmean(time_to_peak.Animal4.dec);

std_peak_time_Animal1_dec = nanstd(time_to_peak.Animal1.dec);
std_peak_time_Animal2_dec = nanstd(time_to_peak.Animal2.dec);
std_peak_time_Animal3_dec = nanstd(time_to_peak.Animal3.dec);
std_peak_time_Animal4_dec = nanstd(time_to_peak.Animal4.dec);


fprintf('--- sig_inc ---\n');
fprintf('A1: %.2f ± %.2f ms\n', mean_peak_time_Animal1_inc, std_peak_time_Animal1_inc);
fprintf('A2: %.2f ± %.2f ms\n', mean_peak_time_Animal2_inc, std_peak_time_Animal2_inc);
fprintf('A3: %.2f ± %.2f ms\n', mean_peak_time_Animal3_inc, std_peak_time_Animal3_inc);
fprintf('A4: %.2f ± %.2f ms\n', mean_peak_time_Animal4_inc, std_peak_time_Animal4_inc);

fprintf('--- sig_dec ---\n');
fprintf('A1: %.2f ± %.2f ms\n', mean_peak_time_Animal1_dec, std_peak_time_Animal1_dec);
fprintf('A2: %.2f ± %.2f ms\n', mean_peak_time_Animal2_dec, std_peak_time_Animal2_dec);
fprintf('A3: %.2f ± %.2f ms\n', mean_peak_time_Animal3_dec, std_peak_time_Animal3_dec);
fprintf('A4: %.2f ± %.2f ms\n', mean_peak_time_Animal4_dec, std_peak_time_Animal4_dec);

% boxplot

figure(Position=[100 100 600 400]);
data_inc = [time_to_peak.Animal1.inc'; time_to_peak.Animal2.inc'; ...
            time_to_peak.Animal3.inc'; time_to_peak.Animal4.inc'];
group_labels_inc = [repmat({'Animal 1'}, length(time_to_peak.Animal1.inc), 1); ...
                    repmat({'Animal 2'}, length(time_to_peak.Animal2.inc), 1); ...
                    repmat({'Animal 3'}, length(time_to_peak.Animal3.inc), 1); ...
                    repmat({'Animal 4'}, length(time_to_peak.Animal4.inc), 1)];
boxplot(data_inc, group_labels_inc, 'colors', [0.5,0.5,0.5]);
ylabel('Time to Peak (ms)');
title('Time to Peak');
colors = [1, 0.4, 0.4;
          0.1, 0.5, 0.8;
          1, 0.7, 0;
          0.7, 0.1, 0.7;
          ];
hBox = findobj(gca, 'Tag', 'Box');
hBox = flipud(hBox);
for j = 1:length(hBox)
    patch(get(hBox(j), 'XData'), get(hBox(j), 'YData'), colors(j, :), ...
          'FaceAlpha', 0.5, 'EdgeColor', colors(j, :),'linewidth',2);
end

set(gca, 'LineWidth', 1, 'FontSize', 18, ...
    'FontName', 'Helvetica', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'ZMinorTick', 'off', ...
    'YLim', [0 2100], ...
    'YTick', 0:1000:2000);
fileName = ['Puff_latency_box_inc'];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


figure(Position=[100 100 600 400]);
data_dec = [time_to_peak.Animal1.dec'; time_to_peak.Animal2.dec'; ...
            time_to_peak.Animal3.dec'; time_to_peak.Animal4.dec'];
group_labels_dec = [repmat({'Animal 1'}, length(time_to_peak.Animal1.dec), 1); ...
                    repmat({'Animal 2'}, length(time_to_peak.Animal2.dec), 1); ...
                    repmat({'Animal 3'}, length(time_to_peak.Animal3.dec), 1); ...
                    repmat({'Animal 4'}, length(time_to_peak.Animal4.dec), 1)];
boxplot(data_dec, group_labels_dec,  'colors', [0.5,0.5,0.5]);
ylabel('Time to Trough (ms)');
title('Time to Trough');



colors = [1, 0.4, 0.4;
          0.1, 0.5, 0.8;
          1, 0.7, 0;
          0.7, 0.1, 0.7;
          ];
hBox = findobj(gca, 'Tag', 'Box');
hBox = flipud(hBox);
for j = 1:length(hBox)
    patch(get(hBox(j), 'XData'), get(hBox(j), 'YData'), colors(j, :), ...
          'FaceAlpha', 0.5, 'EdgeColor', colors(j, :),'linewidth',2);
end

set(gca, 'LineWidth', 1.5, 'FontSize', 18, ...
    'FontName', 'Helvetica', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'ZMinorTick', 'off', ...
    'YLim', [0 2100], ...
    'YTick', 0:1000:2000);



fileName = ['Puff_latency_box_dec'];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% pie

data = [n1_all, n2_all, n3_all];
labels = {'PM', 'NM', 'NSM'};

h = pie(data, labels);

colors = [1 0 0;
          0 0 1;
          0.7 0.7 0.7];

for i = 1:2:length(h)
    set(h(i), 'FaceColor', colors((i + 1) / 2, :));
end
set(findobj(h, 'Type', 'text'), 'fontsize', 24);
percentages = [n1_all, n2_all, n3_all] / sum([n2_all, n3_all, n1_all]) * 100;
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

ht = title({'Modulation (air puff)'}, 'fontsize', 26, 'fontweight', 'normal');
titlePos = get(ht, 'Position');
newTitlePos = titlePos + [0, 0.5, 0];
set(ht, 'Position', newTitlePos);
set(gca, 'Position', [0.13, 0.15, 0.6, 0.6])

fileName = ['airpuff_allMice__pie'];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% Time to 50% increase / decrease

time_to_50_inc = struct('Animal1', [], 'Animal2', [], 'Animal3', [], 'Animal4', []);
time_to_50_dec = struct('Animal1', [], 'Animal2', [], 'Animal3', [], 'Animal4', []);

base_path = 'X:\MFB\Concatenation air puff\air puff';
file_list = dir(fullfile(base_path, 'concatenated_air_puff_*.mat'));

t50_limit = 2000;
dt = 10;

for file_i = 1:length(file_list)
    full_name = file_list(file_i).name;
    tokens = regexp(full_name, 'concatenated_air_puff_(.*)\.mat', 'tokens');
    if ~isempty(tokens)
        file = tokens{1}{1};
        fprintf('Processing: %s\n', file);
    else
        file = erase(erase(full_name, 'concatenated_air_puff_'), '.mat');
    end

    load(fullfile(base_path, full_name));
    speed = all_speed;
    dff_gz = all_dff_gz;

    tpuff = 5e3;
    t_before = 0.5e3;
    t_after = 2e3;

    Tb1 = tpuff - t_before; Tb2 = tpuff;
    Ta1 = tpuff;             Ta2 = tpuff + t_after;

    tis_before = find((T >= Tb1) .* (T < Tb2));
    tis_after  = find((T >= Ta1) .* (T <= Ta2));

    spd_th = 0.05;
    spdm = nanmean(speed(:, [tis_before(1):500, 511:tis_after(end)]), 2);
    ROI_dff_mfa = dff_gz(:, :, :);
    x = squeeze(nanmean(ROI_dff_mfa(:, abs(spdm) < spd_th, :), 2));

    Ncells = size(x, 1);
    t50_inc = nan(1, Ncells);
    t50_dec = nan(1, Ncells);

    for i = 1:Ncells
        x_s = fastsmooth(x(i, :), 10, 3, 1);
        zz_before = x_s(tis_before);
        zz_after  = x_s(tis_after);

        bl_mean    = nanmean(zz_before);
        peak_val   = nanmax(zz_after);
        trough_val = nanmin(zz_after);

        is_inc = false;
        if peak_val > bl_mean
            thresh_50 = bl_mean + 0.5 * (peak_val - bl_mean);
            idx_cross = find(zz_after >= thresh_50, 1, 'first');
            if ~isempty(idx_cross) && idx_cross * dt <= t50_limit
                t50_inc(i) = idx_cross * dt;
                is_inc = true;
            end
        end

        if ~is_inc && trough_val < bl_mean
            thresh_50 = bl_mean + 0.5 * (trough_val - bl_mean);
            idx_cross = find(zz_after <= thresh_50, 1, 'first');
            if ~isempty(idx_cross) && idx_cross * dt <= t50_limit
                t50_dec(i) = idx_cross * dt;
            end
        end
    end

    fprintf('  file = "%s" | Ncells=%d | inc=%d | dec=%d\n', ...
        file, Ncells, sum(~isnan(t50_inc)), sum(~isnan(t50_dec)));

    if strcmp(file, 'animal 1')
        time_to_50_inc.Animal1 = [time_to_50_inc.Animal1, t50_inc];
        time_to_50_dec.Animal1 = [time_to_50_dec.Animal1, t50_dec];
    elseif strcmp(file, 'animal 2')
        time_to_50_inc.Animal2 = [time_to_50_inc.Animal2, t50_inc];
        time_to_50_dec.Animal2 = [time_to_50_dec.Animal2, t50_dec];
    elseif strcmp(file, 'animal 3')
        time_to_50_inc.Animal3 = [time_to_50_inc.Animal3, t50_inc];
        time_to_50_dec.Animal3 = [time_to_50_dec.Animal3, t50_dec];
    elseif strcmp(file, 'animal 4')
        time_to_50_inc.Animal4 = [time_to_50_inc.Animal4, t50_inc];
        time_to_50_dec.Animal4 = [time_to_50_dec.Animal4, t50_dec];
    else
        fprintf('  *** WARNING: "%s" did not match any animal! ***\n', file);
    end
end

fprintf('\n=== Time to 50%% increase ===\n');
fprintf('A1: %.2f +/- %.2f ms (n=%d)\n', nanmean(time_to_50_inc.Animal1), nanstd(time_to_50_inc.Animal1), sum(~isnan(time_to_50_inc.Animal1)));
fprintf('A2: %.2f +/- %.2f ms (n=%d)\n', nanmean(time_to_50_inc.Animal2), nanstd(time_to_50_inc.Animal2), sum(~isnan(time_to_50_inc.Animal2)));
fprintf('A3: %.2f +/- %.2f ms (n=%d)\n', nanmean(time_to_50_inc.Animal3), nanstd(time_to_50_inc.Animal3), sum(~isnan(time_to_50_inc.Animal3)));
fprintf('A4: %.2f +/- %.2f ms (n=%d)\n', nanmean(time_to_50_inc.Animal4), nanstd(time_to_50_inc.Animal4), sum(~isnan(time_to_50_inc.Animal4)));

fprintf('\n=== Time to 50%% decrease ===\n');
fprintf('A1: %.2f +/- %.2f ms (n=%d)\n', nanmean(time_to_50_dec.Animal1), nanstd(time_to_50_dec.Animal1), sum(~isnan(time_to_50_dec.Animal1)));
fprintf('A2: %.2f +/- %.2f ms (n=%d)\n', nanmean(time_to_50_dec.Animal2), nanstd(time_to_50_dec.Animal2), sum(~isnan(time_to_50_dec.Animal2)));
fprintf('A3: %.2f +/- %.2f ms (n=%d)\n', nanmean(time_to_50_dec.Animal3), nanstd(time_to_50_dec.Animal3), sum(~isnan(time_to_50_dec.Animal3)));
fprintf('A4: %.2f +/- %.2f ms (n=%d)\n', nanmean(time_to_50_dec.Animal4), nanstd(time_to_50_dec.Animal4), sum(~isnan(time_to_50_dec.Animal4)));

% --- Across animals ---
all_t50_inc = [time_to_50_inc.Animal1, time_to_50_inc.Animal2, time_to_50_inc.Animal3, time_to_50_inc.Animal4];
all_t50_dec = [time_to_50_dec.Animal1, time_to_50_dec.Animal2, time_to_50_dec.Animal3, time_to_50_dec.Animal4];

fprintf('\n=== Across animals ===\n');
fprintf('Time to 50%% inc: %.2f +/- %.2f ms (n=%d)\n', nanmean(all_t50_inc), nanstd(all_t50_inc), sum(~isnan(all_t50_inc)));
fprintf('Time to 50%% dec: %.2f +/- %.2f ms (n=%d)\n', nanmean(all_t50_dec), nanstd(all_t50_dec), sum(~isnan(all_t50_dec)));

% --- Boxplot: Time to 50% increase ---
colors_box = [1, 0.4, 0.4; 0.1, 0.5, 0.8; 1, 0.7, 0; 0.7, 0.1, 0.7];

figure('Position', [100 100 600 400]);
data_50 = [time_to_50_inc.Animal1'; time_to_50_inc.Animal2'; ...
           time_to_50_inc.Animal3'; time_to_50_inc.Animal4'];
group_labels_50 = [repmat({'Animal 1'}, length(time_to_50_inc.Animal1), 1); ...
                   repmat({'Animal 2'}, length(time_to_50_inc.Animal2), 1); ...
                   repmat({'Animal 3'}, length(time_to_50_inc.Animal3), 1); ...
                   repmat({'Animal 4'}, length(time_to_50_inc.Animal4), 1)];
boxplot(data_50, group_labels_50, 'colors', [0.5, 0.5, 0.5]);
ylabel('Time to 50% increase (ms)');
title('Time to 50% Increase');

hBox = findobj(gca, 'Tag', 'Box');
hBox = flipud(hBox);
for j = 1:length(hBox)
    patch(get(hBox(j), 'XData'), get(hBox(j), 'YData'), colors_box(j, :), ...
          'FaceAlpha', 0.5, 'EdgeColor', colors_box(j, :), 'LineWidth', 2);
end

set(gca, 'LineWidth', 1, 'FontSize', 18, ...
    'FontName', 'Helvetica', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'YLim', [0 2100], ...
    'YTick', 0:500:2000);

exportgraphics(gcf, fullfile(savepath2, 'Time_to_50pct_inc.pdf'), 'ContentType', 'vector');

% --- Boxplot: Time to 50% decrease ---
figure('Position', [100 100 600 400]);
data_50d = [time_to_50_dec.Animal1'; time_to_50_dec.Animal2'; ...
            time_to_50_dec.Animal3'; time_to_50_dec.Animal4'];
group_labels_50d = [repmat({'Animal 1'}, length(time_to_50_dec.Animal1), 1); ...
                    repmat({'Animal 2'}, length(time_to_50_dec.Animal2), 1); ...
                    repmat({'Animal 3'}, length(time_to_50_dec.Animal3), 1); ...
                    repmat({'Animal 4'}, length(time_to_50_dec.Animal4), 1)];
boxplot(data_50d, group_labels_50d, 'colors', [0.5, 0.5, 0.5]);
ylabel('Time to 50% decrease (ms)');
title('Time to 50% Decrease');

hBox = findobj(gca, 'Tag', 'Box');
hBox = flipud(hBox);
for j = 1:length(hBox)
    patch(get(hBox(j), 'XData'), get(hBox(j), 'YData'), colors_box(j, :), ...
          'FaceAlpha', 0.5, 'EdgeColor', colors_box(j, :), 'LineWidth', 2); 
end

set(gca, 'LineWidth', 1, 'FontSize', 18, ...
    'FontName', 'Helvetica', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'YLim', [0 2100], ...
    'YTick', 0:500:2000);

exportgraphics(gcf, fullfile(savepath2, 'Time_to_50pct_dec.pdf'), 'ContentType', 'vector');

