%% Axons Grouping
% -input:
% ROI_dff_all, your raw dff data, size should be #MF x trial x time for each trial
% -output: 
% merged_groups, the group index after manual check
% merged_ldr, the group index after LDR
% merged_cc, the group index after correlation
% 'TP', 'TN', 'FP', 'FN' are evaluation parameters, TP = True positive, 
% TN = True negative, FP = False positive, FN= False negative

dff0 = ROI_dff_all; % #MF x trial x time
Nmf0 = size(dff0,1);
dff0_r = reshape(permute(dff0,[1,3,2]), Nmf0, []); % size of dff0_r = #MF x time
dff0_rz = (dff0_r - nanmean(dff0_r')') ./ nanstd(dff0_r')';

%% Correlation
cc0 = corr(dff0_r', 'rows', 'complete','type','Spearman'); % you can change to pearson
cc0(eye(length(cc0))==1)=nan;
cc_th = prctile(cc0(~isnan(cc0)), 99); % 99th percentile

histogram(cc0(~isnan(cc0)),'FaceColor', [0 0 0],'BinWidth',0.02,'normalization','count')
xline(cc_th, '--r', 'LineWidth', 2, 'Label', '99th percentile');
xlabel('CC');
ylabel('Count');

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

%%
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

%% manual check from results after ldr

manual_check = 1;
manualcheckplot2(dff0_r(:,1:5000), groups_cc, groups_ldr); % dff activity, group after cc, group after ldr
% make sure maunal_groups is in the workspace as the output

%% to merge
merged_cc = mergegroups(groups_cc);
merged_ldr = mergegroups(groups_ldr);

if manual_check == 1
    merged_groups = mergegroups(manual_groups);
else
    merged_groups = mergegroups(groups_ldr);
end

disp('all merged')

%% check if duplicate elements and update
% cc dup
d_cc=report_duplicates(merged_cc);
for i = 1:length(d_cc)
    member = d_cc{i}.Member;
    in_groups = d_cc{i}.Groups;
    group_avg_cc = zeros(1, length(in_groups));
    for j = 1:length(in_groups)
        group_idx = in_groups(j);
        group_members = merged_cc{group_idx};
        others = group_members(group_members ~= member);
        group_avg_cc(j) = mean(cc0(member, others));
    end
    [~, best_group_idx] = max(group_avg_cc);
    best_group = in_groups(best_group_idx);
    
    for j = 1:length(in_groups)
        if in_groups(j) ~= best_group
            merged_cc{in_groups(j)} = merged_cc{in_groups(j)}(merged_cc{in_groups(j)} ~= member);
        end
    end
end
merged_cc = merged_cc(cellfun('prodofsize',merged_cc) ~= 0);


% ldr dup
d_ldr=report_duplicates(merged_ldr);
for i = 1:length(d_ldr)
    member = d_ldr{i}.Member;
    in_groups = d_ldr{i}.Groups;
    group_avg_ldr = zeros(1, length(in_groups));
    for j = 1:length(in_groups)
        group_idx = in_groups(j);
        group_members = merged_ldr{group_idx};
        others = group_members(group_members ~= member);
        group_avg_ldr(j) = mean(ldr(member, others));
    end

    [~, best_group_idx] = min(group_avg_ldr);
    best_group = in_groups(best_group_idx);
    
    for j = 1:length(in_groups)
        if in_groups(j) ~= best_group
            merged_ldr{in_groups(j)} = merged_ldr{in_groups(j)}(merged_ldr{in_groups(j)} ~= member);
        end
    end
end
merged_ldr = merged_ldr(cellfun('prodofsize',merged_ldr) ~= 0);

% manual dup
d_groups=report_duplicates(merged_groups);
for i = 1:length(d_groups)
    member = d_groups{i}.Member;
    in_groups = d_groups{i}.Groups;
    group_avg_groups = zeros(1, length(in_groups));
    for j = 1:length(in_groups)
        group_idx = in_groups(j);
        group_members = merged_groups{group_idx};
        others = group_members(group_members ~= member);
        group_avg_groups(j) = mean(ldr(member, others));
    end

    [~, best_group_idx] = min(group_avg_groups);
    best_group = in_groups(best_group_idx);

    for j = 1:length(in_groups)
        if in_groups(j) ~= best_group
            merged_groups{in_groups(j)} = merged_groups{in_groups(j)}(merged_groups{in_groups(j)} ~= member);
        end
    end
end
merged_groups = merged_groups(cellfun('prodofsize',merged_groups) ~= 0);

%% Manual vs auto
[TP, TN, FP, FN] = manual_vs_auto(merged_groups, merged_ldr, merged_cc);
accuracy = (TP + TN) / (TP + TN + FP + FN);

%% save
full_filename = ['axon_idx', '_', file];
save_path = fullfile('X:\MFB\Processed\Axon_index', full_filename); % change to your dir
save(save_path, 'merged_groups', 'merged_ldr', 'merged_cc', 'TP', 'TN', 'FP', 'FN','accuracy'); 
