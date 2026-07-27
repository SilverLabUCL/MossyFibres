
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   %'191018_13_56_55'; % sessions are removed due to no beh
   %'191018_14_30_00';
   %'191018_14_11_33';
   '191209_13_44_12';
   %'191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22'; % No loco, state, but good whisker  
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   %'200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };

savepath = 'X:\MFB';
fields = fieldnames(results);
numFields = numel(fields)-1;
numLoops =length(folder_names);
num_files = length(folder_names);
N = 20;
all_best_idx = struct();
all_CVEV = struct();
all_XYZ = cell(1, num_files);
all_groups = cell(1, num_files);

for i = 1:numFields
    field = fields{i};
    all_best_idx.(field) = zeros(numLoops, N);
    all_CVEV.(field) = zeros(numLoops, N);
end


for file_i = 1:12
    file = char(folder_names(file_i));
    quickAnalysis
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


    for i = 1:numFields
        field = fields{i};
        best_idx = results.(field).BestMFidx(1:N);
        each_CVEV = results.(field).EachMFCVEV(1:N);

        all_best_idx.(field)(file_i, :) = best_idx;
        all_CVEV.(field)(file_i, :) = each_CVEV;
        
    end
    all_XYZ{file_i} = xyz;
    all_groups{file_i} = group_ids;

end


%% feedback 1
fields = fieldnames(results);

fields = fields(1:end-1);
num_fields = numel(fields);
num_files = size(all_best_idx.(fields{1}), 1);
N = size(all_best_idx.(fields{1}), 2);  % top N MFA

x_all = [];
y_all = [];  % CVEV
g_all = [];

for f = 1:num_fields
    field = fields{f};
    for file_i = 1:num_files
        best_idx = all_best_idx.(field)(file_i, :);
        cvev_vals = all_CVEV.(field)(file_i, :);
        for k = 1:N
            mfa_idx = best_idx(k);
            mfb_list = all_groups{file_i}{mfa_idx};
            len = numel(mfb_list);
            if len > 0
                x_all(end+1) = len;
                y_all(end+1) = cvev_vals(k);
                g_all(end+1) = f;
            end
        end
    end
end

trace_colors = [
    [0 153 76]/255
    [0 204 102]/255
    [0 204 204]/255
    [1.0 0.5 0.1]
    [0.8 0.1 0.6]
    [0.5 0.2 0.7]
    [0.1 0.2 0.7]
    [0.1 0.2 0.4]
];
colors = trace_colors(1:num_fields, :);

figure(Position=[100 100 700 400]); hold on;
for f = 1:num_fields
    idx = (g_all == f);
    scatter(x_all(idx), y_all(idx), 10, colors(f,:), 'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
end

p = polyfit(x_all, y_all, 1);
x_fit = linspace(min(x_all), max(x_all), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '-','color',[0.7,0.7,0.7], 'LineWidth', 2);

text(1.1*mean(x_fit), max(y_all)*0.8, ...
     sprintf('y = %.2fx + %.2f', p(1), p(2)), ...
     'FontSize', 16);

%
xlabel('MFB count per MFA');
ylabel('CVEV');
set(gca, 'FontSize', 18);
box off;
xlim([0 12])


currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Supp.Figures/S8/', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end

fileName = ['Significant MF and MFB counts'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% feedback 2
numFields=8;

overlap_ratio_matrix = zeros(numFields, numFields);

for i = 1:numFields
    for j = 1:numFields
        overlap_ratios = zeros(numLoops, 1);
        for k = 1:numLoops
            idx_i = all_best_idx.(fields{i})(k, :);
            idx_j = all_best_idx.(fields{j})(k, :);
            num_overlap = numel(intersect(idx_i, idx_j));
            overlap_ratios(k) = num_overlap / N;
        end
        overlap_ratio_matrix(i, j) = mean(overlap_ratios);
    end
end

figure('Position', [0.5 0.5 500 500]);
imagesc(overlap_ratio_matrix);colormap('parula')
colorbar;
xticks(1:numFields);
yticks(1:numFields);
xticklabels(fields);
yticklabels(fields);

axis square;
set(gca, 'FontSize', 18);

fileName = ['Mean Overlap Ratio'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% feedback 3 
    trace_colors = [
        [0 153 76]/255
        [0, 204, 102] / 255;
        [0 204 204]/255;
        1.0 0.5 0.1;
        0.8 0.1 0.6;
        0.5 0.2 0.7;
        0.1 0.2 0.7;
        0.1 0.2 0.4;  
    ];


fields = fieldnames(all_best_idx);
num_fields = numel(fields);
num_files = size(all_best_idx.(fields{1}), 1);
N = 20; 

xyz_per_field_per_file = cell(num_fields, num_files);
single_MFA_count = zeros(num_fields, num_files);

for file_i = 1:num_files
    group = all_groups{file_i};
    for f = 1:num_fields
        field = fields{f};
        best_MFA = all_best_idx.(field)(file_i, :);
        xyz_list = [];
        count = 0;
        for k = 1:length(best_MFA)
            mfa_idx = best_MFA(k);
            mfb_list = group{mfa_idx};
            if numel(mfb_list) == 1
                mfb_idx = mfb_list(1);
                xyz_list = [xyz_list; all_XYZ{file_i}(mfb_idx, :)];
                count = count + 1;
            end
        end
        xyz_per_field_per_file{f, file_i} = xyz_list;
        single_MFA_count(f, file_i) = count;
        fprintf('File %d - %s: %d single-MFA found.\n', file_i, field, count);
    end
end



for file_i = 12%:num_files
    figure('Position',[360, 278, 1200, 800]*.5); hold on; view(3);
    all_mfb_xyz = all_XYZ{file_i};
    rotate3d on;
    grid on;

    scatter3(all_mfb_xyz(:,1), all_mfb_xyz(:,2), all_mfb_xyz(:,3), 10, [0.5 0.5 0.5], ...
        'filled', 'MarkerFaceAlpha', 0.2,'DisplayName', 'All MFBs'); 

    colors = trace_colors(1:num_fields-1, :);
    jitter_scale = 0.5;
    
    for f = 1:num_fields-1
        xyz_pts = xyz_per_field_per_file{f, file_i};
        if ~isempty(xyz_pts)
            jitter = (rand(size(xyz_pts)) - 0.5) * 2 * jitter_scale;
            xyz_offset = xyz_pts + jitter;
            scatter3(xyz_offset(:,1), xyz_offset(:,2), xyz_offset(:,3), ...
                     12, colors(f,:), 'filled', 'DisplayName', fields{f});
        end
    end

    legend(Box="off");
    %title(sprintf('File %d - Single-MFA MFB Positions per Behavior', file_i));
    xlabel('X'); ylabel('Y'); zlabel('Z');

    xyz_mx = [250,250, max(all_mfb_xyz(:,3))];
    xyz_mn = [0,0,  max(all_mfb_xyz(:,3))];
    xlim([0,xyz_mx(1)]);
    ylim([0,xyz_mx(2)]);
    %zlim([xyz_mn(3)-10,xyz_mx(3)+10]);
    %zticks([-100, -80, -60, -40, -20, 0]);
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
    grid on

end



%%

B = 500;       
min_sig = 3;   
sig_color = [0.9 0.3 0.3];
band_color = [0.8 0.8 0.8];

include_fields = 1:(num_fields-1);  
special_xyz_all_by_file = cell(num_files,1);

for file_i = 1:num_files
    xyz_all = [];
    for f = include_fields
        xyz_all = [xyz_all; xyz_per_field_per_file{f, file_i}];
    end
    special_xyz_all_by_file{file_i} = unique(xyz_all, 'rows');
end

L_all_files = [];   
L_perm_all_files = [];

for file_i = 1:num_files
    xyz_all = all_XYZ{file_i};                     
    sig_xyz = special_xyz_all_by_file{file_i};     
    
    if isempty(sig_xyz)
        fprintf('File %d: no Sig_MF, skip\n', file_i);
        continue;
    end
    
    z_planes = unique(xyz_all(:,3));
    L_plane_all = []; 
    n_sig_plane = [];   
    L_perm_all_planes = [];
    
    for p = 1:length(z_planes)
        z_val = z_planes(p);
        all_xy = xyz_all(xyz_all(:,3)==z_val, 1:2);
        sig_xy = sig_xyz(sig_xyz(:,3)==z_val, 1:2);
        n_sig = size(sig_xy,1);
        
        if n_sig < min_sig
            continue;
        end
        
        [r_vals, ~, L_obs] = ripleyKL2D(sig_xy, 0);
        L_plane_all = [L_plane_all; L_obs];
        n_sig_plane = [n_sig_plane; n_sig];
        
        L_perm = zeros(B, length(r_vals));
        for b = 1:B
            sel = randperm(size(all_xy,1), n_sig);
            [~,~,L_perm(b,:)] = ripleyKL2D(all_xy(sel,:), 0);
        end
        L_perm_all_planes = cat(1, L_perm_all_planes, L_perm);
    end
    
    if isempty(L_plane_all)
        fprintf('File %d: no valid plane, skip\n', file_i);
        continue;
    end
    
    weights = n_sig_plane / sum(n_sig_plane);
    L_file = sum(L_plane_all .* weights, 1);
    L_all_files = [L_all_files; L_file];
    L_perm_all_files = cat(1, L_perm_all_files, L_perm_all_planes);
end

meanL = mean(L_all_files, 1);
seL = std(L_all_files, 0, 1) / sqrt(size(L_all_files,1));

L_low = prctile(L_perm_all_files, 2.5, 1);
L_high = prctile(L_perm_all_files, 97.5, 1);

figure('Color','w','Position',[100 100 350 350]);
fill([r_vals fliplr(r_vals)], [L_low fliplr(L_high)], ...
    [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold on;
plot(r_vals, meanL, '-', 'Color', sig_color, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 1.5);

xlabel('r (μm)'); ylabel('L(r)');
title('Sig. MF');
set(gca,'FontSize',18,'LineWidth',1,'TickDir','out','Box','off');
