% Select list
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_13_56_55';
   '191018_14_11_33';
   '191018_14_30_00';
   '191209_13_44_12';
   '191209_14_04_14';
   '191209_14_18_13';
   '191209_14_32_39';
   '191209_14_46_58';
   '191209_15_01_22';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   %'200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };


%% cross validated explained variance 

Nigel_c  = [1,.4,.4];
Jeremy_c = [.1,.5,.8];
Bernie_c = [1,.7,0];
Bill_c   = [.7,.1,.7];

ve = cell(length(folder_names),1);
varexp = cell(length(folder_names),1);
varexp_all = cell(length(folder_names),1);
varmax = [];
dimmax = [];
numcomp = nan(length(folder_names),1);

Nsub = 100;

for dataset_i = 1:length(folder_names)

    file = char(folder_names(dataset_i));
    fprintf('\n=== Processing %d/%d: %s ===\n', dataset_i, length(folder_names), file);

    quickAnalysis;
    Mice = getMice(folder_names{dataset_i});

    if strcmp(Mice,'Nigel')
        face_color = Nigel_c;
    elseif strcmp(Mice,'Jeremy')
        face_color = Jeremy_c;
    elseif strcmp(Mice,'Bernie')
        face_color = Bernie_c;
    elseif strcmp(Mice,'Bill')
        face_color = Bill_c;
    else
        face_color = [.5 .5 .5];
    end

    zz = dff_rz;
    nSamples = size(zz,1);

    [ve{dataset_i}, dm, vm] = get_dim(zz, Nsub, Nsub,100);
    fprintf('Shared dimension = %d, peak var = %.3f\n', dm, vm);

    ve_mean = nanmean(ve{dataset_i}, 1);
    ve_sem = nanstd(ve{dataset_i}, 1) / sqrt(size(ve{dataset_i},1));

    varexp{dataset_i} = ve_mean;
    numcomp(dataset_i) = dm;
    varmax = [varmax vm];
    dimmax = [dimmax dm];
    varexp_all{dataset_i} = ve{dataset_i};
end

save([savepath '\varexp.mat'], 've', 'varexp', 'varexp_all', 'varmax', 'dimmax', 'numcomp', 'Nsub', 'folder_names');


%% fig. 4c plot shared PC 
load([savepath '\varexp.mat'])
figure('Position',[100 100 400 350]); hold on;

for dataset_i = 1:length(folder_names)
    if ismember(dataset_i,[3,4,7])
        continue
    end
    file = char(folder_names(dataset_i));
    Mice = getMice(folder_names{dataset_i});

    if strcmp(Mice,'Nigel')
        face_color = Nigel_c;
    elseif strcmp(Mice,'Jeremy')
        face_color = Jeremy_c;
    elseif strcmp(Mice,'Bernie')
        face_color = Bernie_c;
    elseif strcmp(Mice,'Bill')
        face_color = Bill_c;
    else
        face_color = [.5 .5 .5];
    end

    y_mean = varexp{dataset_i};
    y_sem = nanstd(ve{dataset_i}, 1) / sqrt(size(ve{dataset_i},1));
    x = 1:numel(y_mean);
    y_top = y_mean + y_sem;
    y_bot = y_mean - y_sem;
    ix = ~isnan(y_mean);

    fill([x(ix), fliplr(x(ix))], [y_top(ix), fliplr(y_bot(ix))], ...
        face_color, 'LineStyle', 'none', 'FaceAlpha', 0.3);
    plot(x(ix), y_mean(ix), 'k', 'LineWidth', 2);

    [vm, dm] = max(y_mean);
    plot(dm, vm + 0.03, 'v', 'MarkerFaceColor', face_color, 'MarkerEdgeColor', face_color,'MarkerSize',8);
end

xlabel('Num. components');
ylabel('CVEV');
ylim([0, 0.6]);
yticks([0 0.3 0.6]);
xlim([0 40]);
%title('MFA population activity')
set(gca, 'FontSize', 18, 'Box', 'off','LineWidth',1);

currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure4', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end

fileName = ['PCA_cvev'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% fig 4d estimate

varmax_all = []; dimmax_all = []; col_all = [];
for dataset_ix = 1:length(varexp)
    if ~isempty(varexp{dataset_ix})
        for k = 1:70
            [varmax_,dimmax_] = max(varexp_all{dataset_ix}(k,:));
            varmax_all = [varmax_all; varmax_];
            dimmax_all = [dimmax_all; dimmax_];
        end
    end
end

figure('Position', [30 20 400 350])
ix = find(~isnan(varmax_all) & varmax_all>0);
plot(varmax_all(varmax_all>0),dimmax_all(varmax_all>0),'vk','MarkerFaceColor',[.6,.6,.6],'MarkerEdgeColor',[.6,.6,.6])
slope = (varmax_all(ix)' * varmax_all(ix)) \ (varmax_all(ix)' * dimmax_all(ix));
hold on
plot([0, 1], [0, 1] * slope, 'k')
y_fit = slope;
plot([0, 1], [1, 1] * y_fit, ':k')


xlabel('Max variance explained')
ylabel('Lower bound of dim.')
set(gca,'FontSize',18,'LineWidth',1)
box('off')
yticks([0 25 50])
ylim([0 50])
xticks([0 0.5 1])


for dataset_ix = 1:13
    file = folder_names{dataset_ix};
    Mice = getMice(file);

    if strcmp(Mice,'Nigel')
        face_color = Nigel_c;
    elseif strcmp(Mice,'Jeremy')
        face_color = Jeremy_c;
    elseif strcmp(Mice,'Bernie')
        face_color = Bernie_c;
    elseif strcmp(Mice,'Bill')
        face_color = Bill_c;
    end

    if ~isnan(varmax(dataset_ix))
        plot(varmax(dataset_ix),dimmax(dataset_ix),'vk','MarkerFaceColor',face_color,'MarkerEdgeColor',face_color,'MarkerSize',8)
    end
end

fileName = ['PCA_low_bound'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig. 4ab
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_14_30_00';
   '191209_13_44_12';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
};

sig_level = 2;
all_loadings = [];
all_loadings_pm = [];
all_loadings_nm = [];
all_loadings_ns = [];
eigen_20 = [];


concat_folder = 'X:\MFB\MFB_AH_2023\Correlation_data\4mice';
confiles = dir(fullfile(concat_folder, 'concat_*.mat'));

for dataset_ix = 1:length(folder_names)

%     filePath = fullfile(concat_folder, confiles(dataset_ix).name);
%     load(filePath);
    file = char(folder_names(dataset_ix));
    quickAnalysis;

    zz = zscore(dff_r')';
    nan_times = any(isnan(dff_r), 1);  
    L_state = L_state(~nan_times);
    zz = zz(:,~any(isnan(zz), 1));%dff_rz;

    if ~all(L_state==0)
        [cc_wL, zc_wL] = bootstrap_cc(zz', L_state', 100, 300);
    else
        continue
    end

    pm = zc_wL > sig_level;
    [~, pm_ids] = find(pm == 1);
    
    nm = zc_wL < -sig_level;
    [~, nm_ids] = find(nm == 1);

    ns = abs(zc_wL) < sig_level;
    [~, ns_ids] = find(ns == 1);

   % asd = [pm_ids nm_ids];

    [coeff, score, latent, tsquared, explained] = pca(zz(:,:)');

    if size(coeff, 2) < 3
        continue;
    else
        total_variance = sum(latent);
        ne = latent/ total_variance;
        eigen_20 = [eigen_20, ne(1:20)];
        loadings = coeff(:, 1:3) .* sqrt(latent(1:3))';
       

        all_loadings = [loadings; all_loadings];

        if ~isempty(pm_ids)
            all_loadings_pm = [all_loadings_pm; loadings(pm_ids, :)];
        end

        if ~isempty(nm_ids)
            all_loadings_nm = [all_loadings_nm; loadings(nm_ids, :)];
        end
    end

fprintf('Dataset %d: mean NM loading PC1 = %.3f, n_pm=%d, n_nm=%d\n', ...
    dataset_ix, mean(loadings(nm_ids,1)), length(pm_ids), length(nm_ids));

end


%%
figure('Position',[200 100 550 380]); hold on;
for ix = 1:size(eigen_20,2)
    
    Mice = getMice(folder_names{ix});

    if strcmp(Mice,'Nigel')
        face_color = Nigel_c;
    elseif strcmp(Mice,'Jeremy')
        face_color = Jeremy_c;
    elseif strcmp(Mice,'Bernie')
        face_color = Bernie_c;
    elseif strcmp(Mice,'Bill')
        face_color = Bill_c;
    end
    plot(eigen_20(:,ix),'color',face_color,'linewidth',2)
end
mean_eigen = nanmean(eigen_20,2);
sem_mean_eigen = nanstd(eigen_20,0,2) / sqrt(size(eigen_20,2));
errorbar(1:20, mean_eigen, sem_mean_eigen, 'k', 'linestyle', 'none', 'LineWidth', 2);
bar(1:20, mean_eigen, 'FaceAlpha',0, 'LineWidth',1)
xlabel('Number of components')
ylabel('Normalised eigenvalue')
set(gca, 'FontSize', 18, 'XTick', [1 5 10 15 20],'LineWidth',1)
ylim([0, .4])
yticks([0 0.2 0.4])

fileName = ['PCA_normalised_eigenvalue'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


figure('Position',[200 100 430 400]); 
hold on;

for ix = 1:size(eigen_20,2)
    Mice = getMice(folder_names{ix});
    if strcmp(Mice,'Nigel')
        face_color = Nigel_c;
    elseif strcmp(Mice,'Jeremy')
        face_color = Jeremy_c;
    elseif strcmp(Mice,'Bernie')
        face_color = Bernie_c;
    elseif strcmp(Mice,'Bill')
        face_color = Bill_c;
    end

    plot(eigen_20(:,ix),'color',[face_color, 0.5],'linewidth',1.5)
end

mean_eigen = nanmean(eigen_20,2);
sem_mean_eigen = nanstd(eigen_20,0,2) / sqrt(size(eigen_20,2));

plot(1:20, mean_eigen, 'k-', 'LineWidth', 2); 

errorbar(1:20, mean_eigen, sem_mean_eigen, 'k', 'linestyle', 'none', 'LineWidth', 1);

xlabel('Num. components')
ylabel('Normalised eigenvalue')

set(gca, 'FontSize', 18, 'XTick', [1 5 10 15 20], ...
    'LineWidth', 1, 'YScale', 'log'); 

ylim([0.005, 1])
yticks([0.01 0.1 1])
yticklabels({'0.01','0.1','1'})

fileName = ['PCA_normalised_eigenvalue_log'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
%%
% loadings
group_names = {'All', 'MFA^{+}', 'MFA^{-}', 'MFA^{*}'};
colors = {
    [0.85, 0.33, 0.98], ...
    [1, 0, 0], ...
    [0, 0, 1], ...
    [0.5, 0.5, 0.5]
};
figure('Position', [100, 100, 800, 450]);
for pc = 1:3
    data_all = all_loadings(:, pc);
    data_pm  = all_loadings_pm(:, pc);
    data_nm  = all_loadings_nm(:, pc);
    for grp = 1:3
        subplot(3, 3, (pc-1)*3 + grp);
        switch grp
            case 1
                histogram(data_all, 'BinWidth', 0.05, 'Normalization', 'probability', ...
                    'FaceColor', colors{1},'LineWidth',0.01);
            case 2
                histogram(data_pm, 'BinWidth', 0.05, 'Normalization', 'probability', ...
                    'FaceColor', colors{2},'LineWidth',0.01);
            case 3
                histogram(data_nm, 'BinWidth', 0.05, 'Normalization', 'probability', ...
                    'FaceColor', colors{3},'LineWidth',0.01);
        end
        xlabel(sprintf('PC%d', pc));
        ylabel('Frac.');
        title(group_names{grp}, 'Interpreter', 'tex');
        box off
        set(gca, 'FontSize', 12)
        xlim([-1,1])
        ylim([0 0.3])
        yticks([0 0.15 0.3])
    end
end

fileName = ['PCA_loadings_all_pm_nm'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% Calculate dimensionality for chosen subpopulation size (~10 min for each subpopulation size)

savepath = 'X:\MFB';
datapath = 'X:\MFB\Processed\all_data';
file_list = dir(fullfile(datapath, '*.mat'));
N_file = length(file_list);

for sub_i= [75 100 125 150 175]

N_sub = sub_i

acquisition_rate = 100;
dimmax = nan(N_file,1);
varexp = cell(N_file,1);

for dataset_i = 1:N_file
    current_file = fullfile(datapath, file_list(dataset_i).name);
    
    data = load(current_file);

    zz = zscore(data.mf.c.dff_r')';
    Ls = data.mf.b.L_state;
    zz(isnan(zz)) = 0;

    %zz0 = zz(:,1:30000);
    if size(zz,1) >= N_sub
        % Calculate dimensionality for grouped axons
        [varexp{dataset_i},dimmax(dataset_i),~] = get_dim(zz,N_sub,N_sub,100);
    end

end

save([char(savepath),'/Processed/A4dimensionality_N',num2str(N_sub)],'dimmax','varexp')
end


%% fig 4e

valid_datasets = [1,2,4,6,8:16];
Animal3_ids = [2,3,4,5];
Animal4_ids = [6:11];
Animal2_ids = [12:16];
Animal1_ids = [1];

c_Animal1 = [1,.4,.4];
c_Animal2 = [.1,.5,.8];
c_Animal3 = [1,.7,0];
c_Animal4 = [.7,.1,.7];
mice_colors = containers.Map({'Animal1','Animal2', 'Animal3', 'Animal4'}, ...
                                 {c_Animal1, c_Animal2, c_Animal3, c_Animal4});

N_sub = [75:25:175];

%N_sub = [200:100:600];
slope = zeros(size(N_sub));
sem=[];
N_file = 16%length(folder_names);
for k = 1:length(N_sub)
    % Recalculate max variance and dimensionality
    varmax = nan(N_file,1);
    dimmax = nan(N_file,1);
    slope_val = nan(N_file,1);
    for dataset_i =valid_datasets
        load(fullfile(savepath, 'Processed\', ['A4dimensionality_N',num2str(N_sub(k)),'.mat']))
        if ~isempty(varexp{dataset_i}) 
            
            [pks, locs] = findpeaks(nanmean(varexp{dataset_i}), 'MinPeakProminence', 0.01);
            varmax(dataset_i) = pks(1);
            dimmax(dataset_i) = locs(1);


            %[varmax(dataset_i),dimmax(dataset_i)] = max(nanmean(varexp{dataset_i},1));
            val = varmax(dataset_i);

            if isnan(val) || isinf(val) || val == 0
            
                varmax(dataset_i) = nan;
                dimmax(dataset_i) = nan;
            
            end

        end

        slope_val(dataset_i) = dimmax(dataset_i) ./ varmax(dataset_i);
       
    end
    slope_ind = N_sub(k) ./ slope_val;
    ix = find(~isnan(varmax));
    slope(k) = (varmax(ix)'*varmax(ix))\(varmax(ix)'*dimmax(ix));
    sem(k) = nanstd(slope_ind) / sqrt(sum(~isnan(slope_ind)));
end

figure('Position', [30 20 350 350])
bar(N_sub', N_sub./slope, 'FaceColor', [.6, .6, .6], 'EdgeColor', 'w', 'LineWidth', 1)
hold on;
%errorbar(N_sub', N_sub./slope, sem, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
set(gca, 'FontSize', 18, 'Box', 'off','LineWidth',1)
xlim([min(N_sub)-20, max(N_sub)+20])
ylim([0 4])
xlabel('Number of MFAs')
ylabel('MFAs per dimension')

fileName = ['subgroups_' char("errorbar") '_corr'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% Effective dimensionality
% Select
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_14_30_00';
   '191209_13_44_12';
   '191209_14_18_13';
   '191209_14_32_39';
   '191209_14_46_58';
   '191209_15_01_22';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
};

D_eff_all=[];

for file_i  =1:length(folder_names)
    file = char(folder_names(file_i));
    quickAnalysis;

    %dff_p2_z = (dff_p_2 - nanmean(dff_p_2')') ./ nanstd(dff_p_2')';
    
%     beh = [rescale(dis_L2);rescale(dis_R2);rescale(MI_wheel_r);
%         rescale(MI_whisker_r);rescale(wam_r);rescale(wa_r);rescale(wsp_r)];


    valid_time_idx = all(~isnan(dff_rz), 1);
    X_valid = dff_r(:, valid_time_idx);
    X_valid = X_valid - mean(X_valid, 2);

    C = cov(X_valid');
    eigvals = eig(C);
    eigvals = eigvals(eigvals > eps);
  
    lambda_norm = eigvals / sum(eigvals);

    D_eff = 1 / sum(lambda_norm.^2);

    fprintf('Effective dimensionality (normalized PR): %.2f\n', D_eff);

    D_eff_all = [D_eff_all D_eff];
end

mean(D_eff_all)

%% effective dim for motion 
data_dir = 'E:\Hana Video\Loco\infer';
mat_files = dir(fullfile(data_dir, '*.mat'));

D_eff_all=[];

for data_i = 1:length(mat_files)
    file_path = fullfile(data_dir, mat_files(data_i).name);
    fprintf('\nProcessing file: %s\n', mat_files(data_i).name);
    data = load(file_path);

    if isfield(data, 'proc') && isfield(data.proc, 'motSVD') && numel(data.proc.motSVD) >= 2
        Wbeh = data.proc.motSVD{2};

        if isnumeric(Wbeh)
            C = cov(Wbeh);
            eigvals = eig(C);
            eigvals = eigvals(eigvals > eps);
            lambda_norm = eigvals / sum(eigvals);
            D_eff = 1 / sum(lambda_norm.^2);

            fprintf('Effective dimensionality (normalized PR): %.2f\n', D_eff);
        else
            warning('Wbeh is not a numeric matrix, skipping this file.');
        end
    else
        warning('No valid proc.motSVD{2} in file %s, skipping.', mat_files(data_i).name);
    end
       D_eff_all = [D_eff_all D_eff];
end


%% shared dimension

data_dir   = 'E:\Hana Video\Loco\infer';
mat_files  = dir(fullfile(data_dir, '*.mat'));

N_sub   = 125;
D_max   = N_sub;
skip_idx = [3 4 7]; 

varexp = cell(numel(mat_files),1);
varexp_all = cell(numel(mat_files),1);
dimmax = nan(numel(mat_files),1);
varmax = nan(numel(mat_files),1);
filesinfo = struct([]);

for data_i = 1:numel(mat_files)
    if ismember(data_i, skip_idx)
        fprintf('Skip %s (no motion)\n', mat_files(data_i).name);
        filesinfo(data_i).file = [];
        continue
    end

    file_path = fullfile(data_dir, mat_files(data_i).name);
    data = load(file_path);

    fprintf('\nProcessing %s\n', mat_files(data_i).name);
    filesinfo(data_i).file = folder_names{data_i};
    if isfield(data, 'proc') && isfield(data.proc, 'motSVD') && numel(data.proc.motSVD) >= 2
        Wbeh = data.proc.motSVD{2};
        Wbeh = Wbeh(:,1:200);
        Wbeh_r = Wbeh';

        % === shared dimension ===
        [ve, d, v] = get_dim(Wbeh_r, N_sub, D_max,25);

        varexp{data_i}     = nanmean(ve,1);
        varexp_all{data_i} = ve;
        dimmax(data_i)     = d;
        varmax(data_i)     = v;

        fprintf('→ Shared dim: %d | Var explained: %.3f\n', d, v);
    else
        fprintf('No valid motSVD found in file %s\n', mat_files(data_i).name);
    end
end

results.varexp = varexp;
results.varexp_all = varexp_all;
results.dimmax = dimmax;
results.varmax = varmax;
results.filesinfo = filesinfo;
results.N_sub = N_sub;

save(fullfile(savepath, 'motionSVDvarexp.mat'), '-struct', 'results');

%% plot 4f motion energy shared dimension
S = load(fullfile(savepath, 'motionSVDvarexp.mat'));

varexp     = S.varexp;
varexp_all = S.varexp_all;
dimmax     = S.dimmax;
varmax     = S.varmax;
filesinfo  = S.filesinfo;

Nigel_c  = [1,.4,.4];
Jeremy_c = [.1,.5,.8];
Bernie_c = [1,.7,0];
Bill_c   = [.7,.1,.7];

skip_idx = [3 4 7];

figure('Position',[100 100 500 350]); hold on;
for data_i = 1:numel(varexp)
    if isempty(varexp{data_i}) || ismember(data_i, skip_idx), continue; end

    if exist('getMice','file')
        file = filesinfo(data_i).file;
        Mice = getMice(file);
    else
        Mice = 'unknown';
    end
    switch Mice
        case 'Nigel',  face_color = Nigel_c;
        case 'Jeremy', face_color = Jeremy_c;
        case 'Bernie', face_color = Bernie_c;
        case 'Bill',   face_color = Bill_c;
        otherwise,     face_color = [.5 .5 .5];
    end

    ve_mean = varexp{data_i};
    ve_sem  = nanstd(varexp_all{data_i}, 0, 1) ./ sqrt(size(varexp_all{data_i},1));

    x = 1:numel(ve_mean);
    good = isfinite(ve_mean);

    fill([x(good), fliplr(x(good))], ...
         [ve_mean(good)+ve_sem(good), fliplr(ve_mean(good)-ve_sem(good))], ...
         face_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(x(good), ve_mean(good), 'k', 'LineWidth', 2);

    [vmm,dmm] = max(ve_mean(good));
    plot(dmm, vmm+0.03, 'v', 'MarkerFaceColor', face_color, 'MarkerEdgeColor', face_color);
end

set(gca, 'Box','off', 'FontSize', 18,'LineWidth',1);
xlim([0 100]);
ylim([0 0.6]);
yticks([0 0.3 0.6])
xlabel('Number of components');
ylabel('CVEV');
title('Motion energy PC');

fileName = 'motionsharedPC';
exportgraphics(gcf, fullfile(savepath2, [fileName,'.pdf']), 'ContentType', 'vector');


%% shared dimension estimate

c_Animal1 = [1,.4,.4];
c_Animal2 = [.1,.5,.8];
c_Animal3 = [1,.7,0];
c_Animal4 = [.7,.1,.7];

varmax_all = [];
dimmax_all = [];

for dataset_ix = 1:length(varexp)
    if ~isempty(varexp{dataset_ix})
        for k = 1:70
            [varmax_, dimmax_] = max(varexp_all{dataset_ix}(k,:));
            if varmax_ > 0
                varmax_all = [varmax_all; varmax_];
                dimmax_all = [dimmax_all; dimmax_];
            end
        end
    end
end

figure('Position', [30 20 500 350])
ix = find(~isnan(varmax_all));

plot(varmax_all, dimmax_all, 'vk', 'MarkerFaceColor', [.6,.6,.6], 'MarkerEdgeColor', [.6,.6,.6])
hold on

% regression line
slope = (varmax_all(ix)' * varmax_all(ix)) \ (varmax_all(ix)' * dimmax_all(ix));
plot([0, 1], [0, 1] * slope, 'k')

y_fit = slope;
plot([-1, 1], [1, 1] * y_fit, ':k')

xlabel('Max variance explained')
ylabel('Lower bound of dim.')
set(gca,'FontSize',18,'LineWidth',1)
box('off')
yticks([0 50 100])
ylim([0 100])
xticks([0 0.5 1])
xlim([0 1])


for dataset_ix = 1:length(filesinfo)
    file = filesinfo(dataset_ix).file;
    if isempty(file)
        continue
    end
    Mice = getMice(file);

    if strcmp(Mice,'Nigel')
        face_color = c_Animal1;
    elseif strcmp(Mice,'Jeremy')
        face_color = c_Animal2;
    elseif strcmp(Mice,'Bernie')
        face_color = c_Animal3;
    elseif strcmp(Mice,'Bill')
        face_color = c_Animal4;
    end

    if ~isnan(varmax(dataset_ix)) && varmax(dataset_ix) > 0
        plot(varmax(dataset_ix), dimmax(dataset_ix), 'vk', ...
            'MarkerFaceColor', face_color, 'MarkerEdgeColor', face_color)
    end
end

fileName = ['estimate_motionsharedPC'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');



%% PLS on 
rng shuffle

data_dir   = 'E:\Hana Video\Loco\infer';
mat_files  = dir(fullfile(data_dir, '*.mat'));

NCOMP_MAX  = 50;
Kx_use     = 80;
Lblock     = 375;          % interleaved block
fractrain  = 0.8;
skip_idx   = [3 4 7];
tmax       = 380;
shuffle_block = 300;

EV_all   = nan(numel(mat_files), NCOMP_MAX);
EV_shuf  = nan(numel(mat_files), NCOMP_MAX);

for data_i = 1:numel(mat_files)
    if ismember(data_i, skip_idx)
        fprintf('Skip no motion: %s\n', mat_files(data_i).name);
        continue
    end

    file_path = fullfile(data_dir, mat_files(data_i).name);
    fprintf('\n=== Processing file %d/%d: %s ===\n', data_i, numel(mat_files), mat_files(data_i).name);
    load(file_path);

    file = folder_names{data_i};
    quickAnalysis;
    get_forelimb;

    use_beh = (tt >= tmin & tt <= tmax);
    use_neu = (tr >= tmin & tr <= tmax);
    tt2 = tt(use_beh);

    Wbeh  = proc.motSVD{2};
    if size(Wbeh,2) > size(Wbeh,1), Wbeh = Wbeh'; end
    Wbeh2 = Wbeh(use_beh, :);
    wbeh_rs = interp1(tt2, Wbeh2, tr(use_neu), 'linear', 'extrap');   % [T' × K]
    X_full = wbeh_rs(:, 1:min(Kx_use, size(wbeh_rs,2)));

    Y_all = (dff_rz(:, use_neu)).';   % [T' × N]
    NT = size(X_full,1);

    [indtrain, indtest] = splitInterleaved(NT, Lblock, fractrain, 1);

    muX = mean(X_full(indtrain,:));  sdX = std(X_full(indtrain,:));  sdX(sdX==0)=1;
    muY = mean(Y_all(indtrain,:));   sdY = std(Y_all(indtrain,:));   sdY(sdY==0)=1;

    Xtr = (X_full(indtrain,:) - muX)./sdX;
    Ytr = (Y_all(indtrain,:)  - muY)./sdY;
    Xte = (X_full(indtest,:)  - muX)./sdX;
    Yte = (Y_all(indtest,:)   - muY)./sdY;

    meanEV_all = nan(1, NCOMP_MAX);
    for ncomp = 1:NCOMP_MAX
        ncomp_eff = min([ncomp, size(Xtr,2), size(Ytr,2)]);
        [~,~,~,~,BETA] = plsregress(Xtr, Ytr, ncomp_eff);
        Yhat = [ones(sum(indtest),1) Xte] * BETA;

        den = var(Yte, 0, 1); den(den==0)=eps;
        EV_neuron = 1 - var(Yte - Yhat, 0, 1) ./ den;
        meanEV_all(ncomp) = mean(EV_neuron, 'omitnan');
    end
    EV_all(data_i,:) = meanEV_all;

    X_shuf = zeros(size(X_full));
    for col = 1:size(X_full,2)
        X_shuf(:,col) = blockShuffle(X_full(:,col), shuffle_block);
    end

    muXs = mean(X_shuf(indtrain,:));  sdXs = std(X_shuf(indtrain,:));  sdXs(sdXs==0)=1;
    Xtrs = (X_shuf(indtrain,:) - muXs)./sdXs;
    Xtes = (X_shuf(indtest,:)  - muXs)./sdXs;

    meanEV_shuf = nan(1, NCOMP_MAX);
    for ncomp = 1:NCOMP_MAX
        ncomp_eff = min([ncomp, size(Xtrs,2), size(Ytr,2)]);
        [~,~,~,~,BETA] = plsregress(Xtrs, Ytr, ncomp_eff);
        Yhat = [ones(sum(indtest),1) Xtes] * BETA;

        den = var(Yte, 0, 1); den(den==0)=eps;
        EV_neuron = 1 - var(Yte - Yhat, 0, 1) ./ den;
        meanEV_shuf(ncomp) = mean(EV_neuron, 'omitnan');
    end
    EV_shuf(data_i,:) = meanEV_shuf;

    fprintf('File %d done | real peak EV=%.3f, shuffle=%.3f\n', ...
            data_i, max(meanEV_all), max(meanEV_shuf));
end

%%
x = 1:NCOMP_MAX;

EV_mean   = nanmean(EV_all,   1);
EV_sem    = nanstd(EV_all,   0, 1) ./ sqrt(sum(~isnan(EV_all),1));
EVsh_mean = nanmean(EV_shuf,  1);
EVsh_sem  = nanstd(EV_shuf,  0, 1) ./ sqrt(sum(~isnan(EV_shuf),1));

figure(position = [ 100 100 400 350]); hold on;

fill([x, fliplr(x)], [EV_mean+EV_sem, fliplr(EV_mean-EV_sem)], ...
     [0.2 0.5 1.0], 'FaceAlpha', 0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(x, EV_mean, '-', 'Color', [0.2 0.5 1.0], 'LineWidth', 2, 'DisplayName','Real');

[maxEV, idx_max] = max(EV_mean);
plot(x(idx_max), maxEV + 0.008, 'v', ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 6, ...
    'HandleVisibility','off');

% Shuffle
fill([x, fliplr(x)], [EVsh_mean+EVsh_sem, fliplr(EVsh_mean-EVsh_sem)], ...
     [0.6 0.6 0.6], 'FaceAlpha', 0.25, 'EdgeColor','none', 'HandleVisibility','off');
plot(x, EVsh_mean, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'DisplayName','Shuffle');

xlabel('PLS components');
ylabel('CVEV');
box off; xlim([1 NCOMP_MAX]);
legend('Location','best','Box','off');

xlim([0 50])
ylim([ -0.1 0.1])
yticks([-0.1 -0.05 0 0.05 0.1])
set(gca, 'FontSize', 18, 'Box', 'off','LineWidth',1)

fileName = ['PLS_motionPC'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');