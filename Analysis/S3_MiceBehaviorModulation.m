%This code is used to calculate the proportions of PM, NM, and NSM for each concatenated mouse, 
%based on the bootstrapped correlation between activity and different behavior parameters.

clear all;
folderPath =  'X:\MFB\MFB_AH_2023\Correlation_data\4mice';

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
fp = fullfile(mfbFolderPath, 'Figures\Supp.Figures\S3\', currentDate);
if ~exist(fp, 'dir')
    mkdir(fp);
end

% Files
files = dir(fullfile(folderPath, 'concat_*.mat'));

% Color maps
redToWhite = [linspace(1, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)'];
whiteToBlue = [linspace(1, 0, 128)', linspace(1, 0, 128)', linspace(1, 1, 128)'];
redToBlueMap = [redToWhite; whiteToBlue];
mycm = flipud(redToBlueMap);

cl_off = mycm(40, :);
cl_on = mycm(256 - 40, :);
cl_ns = [0.7, 0.7, 0.7];

all_cc_MF_stat = [];
all_cc_MF_run = [];
Nmf = 0;
cc_wl = [];
zc_wl = [];
cc_pw = [];
zc_pw = [];
block_size = 100;
for idx = 1:length(files)
    filePath = fullfile(folderPath, files(idx).name);
    load(filePath);
    dff_rz = (dff_r - nanmean(dff_r, 2)) ./ nanstd(dff_r, 0, 2);
    cc_MF_all = corr(dff_rz', 'rows', 'complete');
    cc_MF_stat = corr(dff_rz(:, L_state==0)', 'rows', 'complete');
    cc_MF_run = corr(dff_rz(:, L_state==1)', 'rows', 'complete');
    Nmf = Nmf + size(dff_rz, 1);
    all_cc_MF_stat = [all_cc_MF_stat; cc_MF_stat(tril(true(size(cc_MF_stat)), -1))];
    all_cc_MF_run = [all_cc_MF_run; cc_MF_run(tril(true(size(cc_MF_run)), -1))];
    %[cc_0, zc_0] = bootstrap_cc(dff_rz', dis_R2_all', block_size, 100);
    %[cc_0, zc_0] = bootstrap_cc(dff_rz', MI_wheel_r', block_size, 100);
    %[cc_0, zc_0] = bootstrap_cc(dff_rz', MI_whisker_r', block_size, 100);
    [cc_0, zc_0] = bootstrap_cc(dff_rz', L_state', block_size, 200);
    
    cc_wl=[cc_wl,cc_0];
    zc_wl=[zc_wl,zc_0];
    
    figure;
    subplot(111);
    hold on;
    [y, x] = histcounts(cc_MF_stat(tril(true(size(cc_MF_stat)), -1)), 'Normalization', 'probability');
    dx = diff(x(1:2));
    xx = x(1:end-1) + dx / 2;
    plot(xx, y, 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7]);
    [y, x] = histcounts(cc_MF_run(tril(true(size(cc_MF_run)), -1)), 'Normalization', 'probability');
    dx = diff(x(1:2));
    xx = x(1:end-1) + dx / 2;
    plot(xx, y, 'LineWidth', 2, 'Color', 'r');
    legend('QW', 'AS', 'Location', 'best');
    legend boxoff;
    xlabel('Pairwise corr.');
    ylabel('Probability');
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out');
    fileName = ['Pairwise_corr_probability_' files(idx).name];
    fullFilePathPDF = fullfile(fp, [fileName,'.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    
    figure('Position', [100, 100, 350, 300]);
    subplot(111); 
    hold on;
    axis off;
    
    sig_level = 2;
    pm_counts = nansum(zc_0 > sig_level);
    ns_counts = nansum(abs(zc_0) < sig_level);
    nm_counts = nansum(zc_0 < -sig_level);
    
    total_counts = pm_counts + ns_counts + nm_counts;
    percentages = [nm_counts, ns_counts, pm_counts] / total_counts * 100;
    
    %labels = {'MFA⁻', 'MFA*', 'MFA⁺'};
    
    h = pie([nm_counts, ns_counts, pm_counts]);
    colormap([cl_off; cl_ns; cl_on]);
    set(findobj(h, 'type', 'text'), 'fontsize', 18);
    
   
    axis image;

    set(gca, 'Position', [0.13, 0.15, 0.6, 0.6]);
    
    fileName = ['modLoc_' files(idx).name '__pie'];
    fullFilePathPDF = fullfile(fp, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
end


%% pie
figure('Position',[100,100,400,350])
subplot(111); hold on
axis off

pm_counts = nansum(zc_wl>sig_level);
ns_counts = nansum(abs(zc_wl)<sig_level);
nm_counts = nansum(zc_wl<-sig_level);
total_counts = pm_counts + ns_counts + nm_counts;
percentages = [nm_counts, ns_counts, pm_counts] / total_counts * 100;


labels = {sprintf('MFA⁻\n%.1f%%', percentages(1)), ...
          sprintf('MFA*\n%.1f%%', percentages(2)), ...
          sprintf('MFA⁺\n%.1f%%', percentages(3))};

h = pie([nm_counts, ns_counts, pm_counts], labels);
colormap([cl_off;cl_ns;cl_on]);
set(findobj(h,'type','text'),'fontsize',18);

textObjects = findobj(h,'Type','text');
for i = 1:length(textObjects)
    pos = get(textObjects(i), 'Position');
    pos(2) = pos(2)+0.05;
    set(textObjects(i), 'Position', pos);
end

axis image;

set(gca, 'Position', [0.13, 0.15, 0.6, 0.6]);

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
fp = fullfile(mfbFolderPath, 'Figures\Supp.Figures\S3\', currentDate);
if ~exist(fp, 'dir')
    mkdir(fp);
end

fileName = ['modLoc_allMice__pie'];
fullFilePathPDF = fullfile(fp, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


% -('Corr. with AS') 

db = .1/2;
bins = [-1-db/2:db:1+db/2];

h = figure('Position',[100,100,300,250]);

hold on
histogram(cc_wl, bins,'FaceColor', [.7,.7,.7], 'FaceAlpha',.5)
histogram(cc_wl(abs(zc_wl)>sig_level),bins, 'DisplayStyle','stairs','edgeColor', 'k','LineWidth',2)

xlabel('Corr. with Loco')
ylabel('# MFAs')

xlim([-1,1]);

set(gca, 'LineWidth', 1, 'FontSize', 15)

fileName = ['Corr with AS.png'];
fullFilePath = fullfile(fp, fileName);
print(fullFilePath, '-dpng', '-r300');

h = figure('Position',[100,100,300,250]); hold on;
sig_level=2;
db = 1;

zz = zc_0(zc_0>sig_level);
bins = [min(zz(:))-db/2:db:max(zz(:))+db/2];
histogram(zz, bins, 'FaceColor', cl_on, 'FaceAlpha',.85)

zz = zc_0(zc_0<-sig_level);
bins = [min(zz(:))-db/2:db:max(zz(:))+db/2];
histogram(zz, bins, 'FaceColor', cl_off, 'FaceAlpha',.85)

zz = zc_0(abs(zc_0)<sig_level);
bins = [min(zz(:))-db/2:db:max(zz(:))+db/2];
histogram(zz, bins, 'FaceColor', cl_ns, 'FaceAlpha',.85)

xlabel('Corr. with Loco')
ylabel('# MFAs')

set(gca, 'LineWidth', 1, 'FontSize', 15)

fileName = ['Bootstrapped Corr_all'];
fullFilePathPDF = fullfile(fp, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


figure('Position',[100,100,300,250])
subplot(111)
hold on

[y, x] = histcounts(all_cc_MF_stat, 'Normalization', 'probability');
dx = diff(x(1:2));
xx = x(1:end-1) + dx/2;
plot(xx, y, 'LineWidth', 2, 'Color', [50, 200, 200]/255);
hold on;

[y, x] = histcounts(all_cc_MF_run, 'Normalization', 'probability');
dx = diff(x(1:2));
xx = x(1:end-1) + dx/2;
plot(xx, y, 'LineWidth', 2, 'Color', [200 50 200]/255);

legend('QW', 'AS', 'Location', 'best')
legend boxoff

xlabel('Pairwise corr.')
ylabel('Probability')
set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out','box','off')

fileName = ['Pairwise corr_probability'];
fullFilePathPDF = fullfile(fp, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
