%% load data
file ='200130_13_36_14 FunctAcq';
quickAnalysis
dff0 = ROI_dff_all;
Nmf0 = size(dff0,1);
dff0_r = reshape(permute(dff0,[1,3,2]), Nmf0, []);
dff0_rz = (dff0_r - nanmean(dff0_r')') ./ nanstd(dff0_r')';

%% S2a
cc0 = corr(dff0_r', 'rows', 'complete','type','Spearman'); 
cc0(eye(length(cc0))==1)=nan;
cc_th = prctile(cc0(~isnan(cc0)), 99); % 99th percentile
figure(position = [100 100 350 300])
histogram(cc0(~isnan(cc0)),'FaceColor', [0 0 0],'BinWidth',0.02,'normalization','count')
xline(cc_th, '--r', 'LineWidth', 2, 'Label', '99th percentile');
xlabel('Spearman CC');
ylabel('Count');
yticks([0 6000 12000])
set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'ZMinorTick'  , 'off'  ...
    )
mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Supp.Figures/S2/', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
fileName = ['example correlation distribution'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

groups_cc = {};
pregroups1 = {};

for i = 1:Nmf0
    pregroups1{i} = i;

    for j = 1:Nmf0
        if i ~= j && cc0(i, j) > cc_th
            pregroups1{i} = [pregroups1{i}, j]; 
        end
    end
end

for i = 1:Nmf0
    if numel(pregroups1{i}) <= 2 
        groups_cc{i} = pregroups1{i};
        continue;
    end

    current_group = pregroups1{i};
    primary = current_group(1);
    cc_members = primary;

    for j = 2:length(current_group)
        for k = j+1:length(current_group)
            if cc0(current_group(j), current_group(k)) < cc_th
                if cc0(primary, current_group(j)) >= cc0(primary, current_group(k)) 
                    current_group(k) = [];
                else
                    current_group(j) = [];
                end
                break; 
            end
        end
        if numel(current_group) < 3
            break; 
        end
    end
    groups_cc{i} = current_group;
end


%% LDR
load(['X:\MFB\Smooth_data\', file, '.mat']);
load(['X:\MFB\LDR\', file, '.mat']);

ldr_th = 1.5; 
pregroups2 = cell(1, Nmf0);

for i = 1:Nmf0
    primary = i;
    pregroups2{i} = primary;
    for j = 2:length(groups_cc{i})
        if ldr(primary, groups_cc{i}(j)) <= ldr_th
            pregroups2{i} = [pregroups2{i}, groups_cc{i}(j)];
        end
    end
end

groups_ldr = cell(1, Nmf0);

for i = 1:Nmf0
    if numel(pregroups2{i}) <= 1
        groups_ldr{i} = pregroups2{i};
        continue;
    end

    current_group = pregroups2{i};
    primary = current_group(1); 
    ldr_members = primary;

    for j = 2:length(current_group)
        isgroup = true;
        for k = 2:length(current_group)
            if j ~= k && ldr(current_group(j), current_group(k)) > ldr_th
                isgroup = false;
                if ldr(primary, current_group(j)) < ldr(primary, current_group(k))
                    ldr_members = [ldr_members, current_group(j)];
                else
                    ldr_members = [ldr_members, current_group(k)];
                end
                break;
            end
        end
    
        if isgroup && ~ismember(current_group(j), ldr_members)
            ldr_members = [ldr_members, current_group(j)];
        end
    end
    groups_ldr{i} = ldr_members; 
end

figure(position = [100 100 350 300])
plot(cc0,ldr,'o','color',[0.4 0.4 0.4])
yline(ldr_th, '--r', 'LineWidth', 2, 'Label', 'LDR = 1.5');
xlabel('Spearman CC');
ylabel('LDR');
ylim([0 10])
set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'ZMinorTick'  , 'off'  ...
    )
fileName = ['LDR threshold'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% S2c
m = [92,88]% [226 244]%[92,88]%;
fig = figure('visible', 'on', 'position', [100 100 1000 250]);

offset = 1* range(dff0_r(m(1),:)) + 1.4*range(dff0_r(m(2),:));

trace1 = dff0_r(m(1), :);
trace2 = dff0_r(m(2), :) - offset;

x_limits = [1, size(dff0_r, 2)];
y_limits = [min(trace2), max(trace1)];

ax = axes(fig);
hold on;
plot(trace1, 'k');
plot(trace2, 'k');
xlim(x_limits);
ylim(y_limits);
axis off;

% Scale bar
x_start = x_limits(2) - 1100;
y_start = y_limits(1)+3.8;
line([x_start, x_start+1000], [y_start, y_start],   'Color','k','LineWidth',2,'Clipping','off');
line([x_start, x_start],      [y_start, y_start+0.5],  'Color','k','LineWidth',2,'Clipping','off');
text(x_start+500, y_start, '10 s',        'FontSize',15,'HorizontalAlignment','center','VerticalAlignment','top');
text(x_start-50,  y_start+0.5, '0.5 \DeltaF/F', 'FontSize',15,'HorizontalAlignment','right');

hold off;


inset_range = 38957:39617;
y_box = [min(min(trace1(inset_range)), min(trace2(inset_range))), ...
         max(max(trace1(inset_range)), max(trace2(inset_range)))];
rectangle('Position', [inset_range(1), y_box(1), length(inset_range), diff(y_box)], ...
          'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--');


exportgraphics(fig, fullfile(folderPath, 'example_with_scale_bar.pdf'), 'ContentType','vector');

fig_inset = figure('visible', 'on', 'position', [100 100 70 80]);
hold on;
plot(inset_range, trace1(inset_range), 'k', 'LineWidth', 0.5);
plot(inset_range, trace2(inset_range) + offset, 'r', 'LineWidth', 0.5);
xlim([inset_range(1), inset_range(end)]);
axis off;

inset_xl = xlim; inset_yl = ylim;
sb_x = inset_xl(2) - 120;
sb_y = inset_yl(2) - 0.05 * diff(inset_yl);
sb_t = 100;    % 100 points = 1 s
sb_df = 0.2;   % 0.2 dF/F
line([sb_x, sb_x + sb_t], [sb_y, sb_y], 'Color', 'k', 'LineWidth', 1, 'Clipping', 'off');
line([sb_x, sb_x], [sb_y - sb_df, sb_y], 'Color', 'k', 'LineWidth', 1, 'Clipping', 'off');
text(sb_x + sb_t/2, sb_y, '1 s', 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(sb_x - 5, sb_y - sb_df, '0.2 \DeltaF/F', 'FontSize', 10, 'HorizontalAlignment', 'right');
hold off;

exportgraphics(fig_inset, fullfile(folderPath, 'mismatch_zoom.pdf'), 'ContentType', 'vector');



zz0 = dff0_r_smooth;
zz0(isnan(zz0)) = 0;
i = m(1); j = m(2);
a = get_var_ratio(zz0(i,:), zz0(j,:), 1);

figs_before = findobj('type', 'figure');
a = get_var_ratio(zz0(i,:), zz0(j,:), 1);
figs_after = findobj('type', 'figure');

new_figs = setdiff(figs_after, figs_before);




for k = 1:length(new_figs)
    exportgraphics(new_figs(k), fullfile(folderPath, ['var_ratio_fig', num2str(k), '.pdf']), 'ContentType', 'vector');
end


