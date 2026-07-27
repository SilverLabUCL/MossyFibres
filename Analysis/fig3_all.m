clear all;close all;clc
path_home = 'X:\MFB\MFB_AH_2023';
addpath([path_home '/Init']);
addpath([path_home '/Analysis/Analysis_pipeline']);
addpath([path_home '/Analysis/Analysis_pipeline/Utilities']);
addpath 'X:\MFB\MFB_AH_2023\Init'
run 'X:\MFB\MFB_AH_2023\Init\Initialize';

%% 3a. heat map of activity
load('X:\MFB\MFB_AH_2023\Correlation_data\4mice\concat_Animal 4.mat')

redToWhite = [linspace(1, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)'];
whiteToBlue = [linspace(1, 0, 128)', linspace(1, 0, 128)', linspace(1, 1, 128)'];
redToBlueMap = [redToWhite; whiteToBlue];

mycm=flipud(redToBlueMap);
suffix = 'changeSorted';

Nmf = size(dff_r,1);
Ntot = Nmf;
dt = 10;
%
figure('Position', [30 20 1200 550])
sbplt1 = [1:8];
sbplt2 = [9];
sbplt3 = [10];
sbp_no = 10;    
dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';  
ccs = corr(L_state', dff_rz', 'rows', 'complete');
[a,b]=sort(ccs);
zz = dff_rz(b,:);

% clean
nans_t = any(isnan(dff_rz), 1);
dff_rz(:, nans_t) = [];
MI_whisker_r(:, nans_t) = [];
L_state(:, nans_t) = [];
A_state(:, nans_t) = [];
Q_state(:, nans_t) = [];
MI_wheel_r(:, nans_t) = [];


zz = dff_rz(b,:);
sorted_zz = sort(zz(~isnan(zz)), 'ascend');
n = numel(sorted_zz);
lower_i = max(floor(n * 0.01), 1);
upper_i = max(floor(n * 0.99), 1);
lv = round(sorted_zz(lower_i), 1);
uv = round(sorted_zz(upper_i), 1);

xlm = [0,size(zz,2)];
subplot(sbp_no,1,sbplt1); hold on
imagesc(zz, [lv,uv]);


ylb = ylabel('Mossy fibre axon number');
set(ylb, 'Units', 'Normalized', 'Position', [-.04, 0.5, 0]);
yticks([1, Ntot])
ylim([0+.5,Ntot+.5])
xticklabels([]);
xlim(xlm)
cb = colorbar();
set(cb, 'Position', [0.92, 0.65, 0.012, 0.05], 'Ticks', [lv, uv], 'FontSize', 10);
titleStr = {'Z-scored', '\DeltaF/F'};
title(cb, titleStr, 'FontSize', 10, 'FontWeight', 'normal');

pos = get(gca, 'Position');
pos(1) = 0.15; pos(3) = 0.75;
pos(2) = .26; pos(4) = .4;
set(gca, 'Position', pos, 'LineWidth', 1, 'FontSize', 12, 'XColor', 'none', 'YDir', 'normal','TickDir', 'out')
time = (1:length(L_state)) ; % 100 = 1s

% Subplot for MI_whisker_r
subplot(sbp_no,1,sbplt3) 
hold on
wsk_r = rescale(MI_whisker_r);
ylb = ylabel({'WMI'});
set(ylb, 'Units', 'Normalized', 'Position', [-0.04, 0.5, 0]);
pos = get(gca, 'Position');
pos(1) = 0.15; pos(3) = 0.75;
pos(2) = .033; pos(4) = .075;

set(gca, 'Position', pos,'LineWidth', 1, 'FontSize', 12, 'XColor', 'none','TickDir', 'out')
xlim(xlm)
ylim([-0.05 1])
yticks([0 1])
AS_regions = A_state == 1;
AS_start = find(diff([0; AS_regions(:)]) == 1);
AS_end = find(diff([AS_regions(:); 0]) == -1);
for i = 1:length(AS_start)
    patch([time(AS_start(i)), time(AS_end(i)), time(AS_end(i)), time(AS_start(i))], ...
          [-0.05, -0.05, 1, 1], ...
          [255 205 255]/255, 'EdgeColor', 'none');
end
QW_regions = Q_state == 1;
QW_start = find(diff([0; QW_regions(:)]) == 1);
QW_end = find(diff([QW_regions(:); 0]) == -1);
for i = 1:length(QW_start)
    patch([time(QW_start(i)), time(QW_end(i)), time(QW_end(i)), time(QW_start(i))], ...
          [-0.05, -0.05, 1, 1], ...
          [192 255 255]/255, 'EdgeColor', 'none');
end

plot(time, wsk_r, 'LineWidth',1, 'color', [255,163,26]/255)

hold on;
x0 = max(xlim)/4.5; 
y0 = min(wsk_r) ; 
plot(x0-[0,6000], [y0-0.05,y0-0.05], 'k-', 'LineWidth',3);
text(x0-3000, y0-0.35, '60 s', 'FontSize', 12, 'HorizontalAlignment', 'center'); 

% Subplot for Right limb
subplot(sbp_no,1,sbplt2)
hold on
ylb = ylabel({'RFL'});
set(ylb, 'Units', 'Normalized', 'Position', [-.04, 0.5, 0]);

pos = get(gca, 'Position');
pos(1) = 0.15; pos(3) = 0.75;
pos(2) = .14; pos(4) = .075;
set(gca, 'Position', pos)
xlim(xlm)
ylim([-0.05 1])
yticks([0 1])
for i = 1:length(AS_start)
    patch([time(AS_start(i)), time(AS_end(i)), time(AS_end(i)), time(AS_start(i))], ...
          [-0.05, -0.05, 1, 1], ...
          [255 205 255]/255, 'EdgeColor', 'none');
end

for i = 1:length(QW_start)
    patch([time(QW_start(i)), time(QW_end(i)), time(QW_end(i)), time(QW_start(i))], ...
          [-0.05, -0.05, 1, 1], ...
          [192 255 255]/255, 'EdgeColor', 'none');
end
plot(time, rescale(dis_R2_all), 'color', [173,210,157]/255, 'LineWidth',1);
set(gca, 'LineWidth', 1, 'FontSize', 12, 'XColor', 'none','TickDir', 'out')

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
fP = fullfile(mfbFolderPath, 'Figures', 'Figure3', currentDate);
if ~exist(fP, 'dir')
    mkdir(fP);
end

fileName = ['HeatMap__' suffix];
fullFilePathPDF = fullfile(fP, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% 3b.pc

folderpath= 'X:\MFB\MFB_AH_2023\Correlation_data\pnp';
files = dir(fullfile(folderpath, 'concat_*.mat'));

idx = 4;
filePath = fullfile(folderpath, files(idx).name);
load(filePath);

dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')'; 
zz = dff_rz;

zz = zscore(dff_r')';
dt = 10;
nsh = 10;
fs = 18;
zzf = fillmissing(zz, 'spline', 5,'MaxGap',10);

% if idx ==1
%     zzf = zzf(:,100:47120);
%     A_state = A_state(100:47120);
%     Q_state = Q_state(100:47120);
% end

[coeff, score, LATENT, TSQUARED, EXPLAINED] = pca(zzf', 'NumComponents',3);
dFF_L = zzf(:, A_state==1);
dFF_Q = zzf(:, A_state==0);              

mu_L = nanmean(dFF_L,2);
mu_Q = nanmean(dFF_Q,2);

w = mu_L - mu_Q;
w = w / norm(w);

A_or_Q = w'*(zzf - nanmean(zzf,2));
A_or_Q = (A_or_Q - min(A_or_Q)) / (max(A_or_Q) - min(A_or_Q)); 

c = nan(size(A_or_Q,2),3);
% ls = 1:nsh:length(A_state); % lower sampling
% score = score(ls, :);
% A_or_Q = A_or_Q(ls);

win = nsh;
score_sm = movmean(score, win, 1);
A_or_Q_sm = movmean(A_or_Q, win);

ds = nsh;
ls = 1:ds:size(score_sm, 1);
score = score_sm(ls, :);
A_or_Q = A_or_Q_sm(ls);


figure('Position',[100,100,350,300]); hold on;
for k = 1:1:length(score)-1
    c_ = A_or_Q(k);
    if c_ < 0
        c_ = 0;
    elseif c_ > 1
       c_ = 1;
    end
    c(k,:) = c_*[1,0,1] + (1-c_)*[0,1,1];

    if ~isnan(c(k,:))
        plot3(score(k:k+1,1),score(k:k+1,2),score(k:k+1,3),'-','Color',c(k,:),'LineWidth', 0.75);
    end
end

if idx == 1
    view(119,-48)
elseif idx == 2
    view(109,-67)
elseif idx==3  
    view(-152,69)
elseif idx==4
    view(-152,69)
elseif idx== 5 || idx==6||idx ==7
    view(74,68)
elseif idx==8
    view(115,-58)
elseif idx==9
    view(-72,26)
elseif idx==10
    view(-73,-51)
end

xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength',[.025,.025],'XTick',[],'YTick',[],'ZTick',[])


savepath = 'X:/MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure3', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end


fileName = ['pc_activitiy_proj_' files(idx).name(1:end-4)];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%%
set(gcf, 'Position', [100, 100, 800, 600]);
daspect([0.5 1 0.5]);
axis vis3d;

videoFile = ['X:\MFB\Figures\Figure3\Video_' files(idx).name(1:end-4)];
v = VideoWriter(videoFile, "MPEG-4");
v.FrameRate = 30;
v.Quality = 100;
open(v);

numFrames = 180;
[az0, ~] = view;
el_fixed = 35;
angle1 = az0 + (0:numFrames-1) * 360 / numFrames;

for t = 1:numFrames
    view(angle1(t), el_fixed);
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Video saved: %s.mp4\n', videoFile);
%% 3c. manifold_distance


%Jeremy
Jeremy_names = {
    '191018_13_39_41';
    '191018_13_56_55';
    '191018_14_30_00';
    '191018_14_11_33';
    };

%Bernie
Bernie_names = {
   '191209_13_44_12';
   '191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
    };

%Nigel
Nigel_names = {
    '171212_16_19_37';
    };

%Bill
folder_names = {
    '200130_13_21_13 FunctAcq';
    '200130_13_36_14 FunctAcq';
    '200130_13_49_09 FunctAcq';
    %'200130_14_02_12 FunctAcq';
    '200130_14_15_24 FunctAcq';
    '200130_14_29_30 FunctAcq';
    };

%Select
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   %'191018_13_56_55';
   '191018_14_30_00';
   %'191018_14_11_33';
   '191209_13_44_12';
   %'191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   %'200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };




angle_LQ_val = [];
angle_LQ_sh = [];
angle_LQ_sh_bl = [];
Ns = [];

dist_L = [];
dist_Q = [];
dist_LQ = [];

ijks = [];

for ii = 1:length(folder_names)
    %
    file=char(folder_names(ii));
    quickAnalysis;
    block_size = 100;     
    dt = 10;
    
    [N,T] = size(dff_r);
    Ns = [Ns, N];
    
    dFF_L = dff_rz(:, A_state==1);
    dFF_Q = dff_rz(:, Q_state==1);  
    dsr = 100;
    xL = dFF_L(:,1:dsr:end)';
    xQ = dFF_Q(:,1:dsr:end)';
    
    if size(dFF_L,2)<100
        xL = dFF_L(:,1:end)';
    end

    x = pdist2(xL, xL);
    dist_L = [dist_L, nanmean(x(:))];
    
    x = pdist2(xQ, xQ);
    dist_Q = [dist_Q, nanmean(x(:))];
    
    x = pdist2(xL, xQ);
    dist_LQ = [dist_LQ, nanmean(x(:))];

end
%%
figure('Position', [100,100,350,350])
avg_L = nanmean(dist_L);
sem_L = nanstd(dist_L) / sqrt(sum(~isnan(dist_L)));
avg_Q = nanmean(dist_Q);
sem_Q = nanstd(dist_Q) / sqrt(sum(~isnan(dist_Q)));
avg_LQ = nanmean(dist_LQ);
sem_LQ = nanstd(dist_LQ) / sqrt(sum(~isnan(dist_LQ)));
bar_positions = 1:3;

b = bar(bar_positions, [avg_L, avg_Q, avg_LQ], 'FaceColor', 'none', 'EdgeColor', 'black','LineWidth', 1);

hold on;

for iii = 1:length(dist_L)
    plot(bar_positions, [dist_L(iii), dist_Q(iii), dist_LQ(iii)], 'o-', 'Color', [0.7, 0.7, 0.7]);
    hold on
end

scatter(ones(size(dist_L)) * bar_positions(1), dist_L, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'magenta','LineWidth', 1.5);
scatter(ones(size(dist_Q)) * bar_positions(2), dist_Q, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'cyan','LineWidth', 1.5);
scatter(ones(size(dist_LQ)) * bar_positions(3), dist_LQ, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.5,0.5,0.5],'LineWidth', 1.5);

hErrorbar = errorbar(bar_positions, [avg_L, avg_Q, avg_LQ], [sem_L, sem_Q, sem_LQ], 'k', 'linestyle', 'none');
uistack(hErrorbar, 'top');
fs = 15;
xlim([0,4])
ylim([10 30])
yticks([10 20 30])
box('off')
ylabel('Avg. distance')
xticks([1,2,3])
xticklabels({'AS', 'QW', 'AS-QW'})
xlabel({'Within/between manifolds'})
set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength',[.025,.025])

box("off")


[p_LQ, h_LQ] = signrank(dist_L, dist_Q);
[p_L_LQ, h_L_LQ] = signrank(dist_L, dist_LQ);
[p_Q_LQ, h_Q_LQ] = signrank(dist_Q, dist_LQ);

disp(['AA and QQ: p value = ', num2str(p_LQ)]);
disp(['AA and AQ: p value = ', num2str(p_L_LQ)]);
disp(['QQ and AQ: p value  = ', num2str(p_Q_LQ)]);

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
fP = fullfile(mfbFolderPath, 'Figures', 'Figure3', currentDate);
if ~exist(fP, 'dir')
    mkdir(fP);
end

fileName = ['Manifold_distance_all_np'];
fullFilePathPDF = fullfile(fP, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% 3d Manifold Angle

savepath ='X:\MFB';
datapath = 'X:\MFB\Processed\all_data';
file_list = dir(fullfile(datapath, '*.mat'));
N_file = length(file_list);

num_PCs = 3;
angle_LQ_val = [];
angle_LQ_min = [];
angle_LQ_sh_bl_max = [];
angle_LQ_sh_bl_min = [];
Ns = [];
dt = 10;
good_sessions = [1,2,4,6,8:16];
for i = good_sessions

    current_file = fullfile(datapath, file_list(i).name);
    
    data = load(current_file);

    dff_rz = zscore(data.mf.c.dff_r')';
    L_state = data.mf.b.L_state;

    file = data.mf.info.file;
    [N, T] = size(dff_rz);
    Ns = [Ns, N];

    if all(L_state == 0)
        continue
    end

    dFF_L = dff_rz(:, L_state == 1);
    dFF_Q = dff_rz(:, L_state == 0);
    
    coeff_L = pca(dFF_L');
    coeff_Q = pca(dFF_Q');
    
    A_metric = eye(size(coeff_Q,1));
    [theta, ~, ~] = subspacea(coeff_Q(:,1:num_PCs), coeff_L(:,1:num_PCs), A_metric);
    angle_LQ_val = [angle_LQ_val, max(theta)];
    angle_LQ_min = [angle_LQ_min, min(theta)];

    ashes_max = [];
    ashes_min = [];
    N_itr = 20;
    for j = 1:N_itr
        bl_sz = 100;
        dff_sh_bl = dff_rz(:, randblock(1:floor(T / bl_sz) * bl_sz, [1, bl_sz]));
        assert(mod(size(dff_sh_bl,2), 2) == 0, 'Shuffled data must have even columns.');
        
        dff_1 = dff_sh_bl(:, 1:size(dff_sh_bl, 2)/2);
        dff_2 = dff_sh_bl(:, size(dff_sh_bl, 2)/2+1:end);
        
        coeff_1 = pca(dff_1');
        coeff_2 = pca(dff_2');
        
        [theta_sh, ~, ~] = subspacea(coeff_1(:,1:num_PCs), coeff_2(:,1:num_PCs), A_metric);
        ashes_max = [ashes_max, max(theta_sh)];
        ashes_min = [ashes_min, min(theta_sh)];
    end

    angle_LQ_sh_bl_max = [angle_LQ_sh_bl_max, mean(ashes_max, 'omitnan')];
    angle_LQ_sh_bl_min = [angle_LQ_sh_bl_min, mean(ashes_min, 'omitnan')];
end

[p_max, ~, ~] = signrank(angle_LQ_val, angle_LQ_sh_bl_max);
[p_min, ~, ~] = signrank(angle_LQ_min, angle_LQ_sh_bl_min);
disp(['Max Angle p = ', num2str(p_max)]);
disp(['Min Angle p = ', num2str(p_min)]);
%%

figure('Position', [100, 100, 300, 350])
hold on

h1 = plot(nan, nan, 'o', 'LineWidth', 2, 'Color', [0,0.5,0.4], 'MarkerFaceColor', [0,0.5,0.4]);
h2 = plot(nan, nan, 'o', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.7 0.7 0.7], 'LineWidth', 1);

for i = 1:length(angle_LQ_val)
    plot([1, 2], [angle_LQ_val(i), angle_LQ_sh_bl_max(i)], '-o', ...
        'LineWidth', 2, 'Color', [0,0.5,0.4], 'MarkerFaceColor', [0,0.5,0.4])
    plot([1, 2], [angle_LQ_min(i), angle_LQ_sh_bl_min(i)], ...
        'Color', [0.7 0.7 0.7], 'Marker', 'o', 'LineWidth', 0.5, ...
        'MarkerFaceColor', [0.7 0.7 0.7])
end

legend([h1, h2], {'Max Angle', 'Min Angle'}, 'Location', 'best','Box','off');


xlim([0.5, 2.5])
ylim([0, 2])
xticks([1, 2])
xticklabels({'QW-AS', 'Shuffle'})
xtickangle(45)
fs = 15;
set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [.025, .025])
ylabel('Angle (rad.)')

fileName = 'manifold_angle_all';
fullFilePathPDF = fullfile(fP, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% 3e manifold_angle_percentile both PM/NM exclude 
num_PCs = 3;
Ns = [];
prcnts = 10:10:90;
angles_LQ_min = [];
angles_LQ_max = [];
angles_LQ_shuffle_min = [];
angles_LQ_shuffle_max = [];

fprintf('Total files to process: %d\n', length(good_sessions));

good_sessions = [1,2,4,6,8:16];
for file_i = 1:length(good_sessions)
    fprintf('\n=== Processing file %d/%d ===\n', file_i, length(good_sessions));
    
    current_file = fullfile(datapath, file_list(file_i).name);
    fprintf('File: %s\n', file_list(file_i).name);
    
    data = load(current_file);
    dff_rz = zscore(data.mf.c.dff_r')';
    L_state = data.mf.b.L_state;
    dt = 10;
    block_size = 100;
    
    if ~all(L_state==0)
        [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', block_size, 200);
    else
        fprintf('Skipping file (all L_state==0)\n');
        continue
    end
    
    [N,T] = size(dff_rz);
    Ns = [Ns, N];
    vals = prctile(abs(zc_wL), 100 - prcnts);
    
    angle_LQ_min_current = nan(1, length(prcnts));
    angle_LQ_max_current = nan(1, length(prcnts));
    angle_LQ_sh_min = nan(1, length(prcnts));
    angle_LQ_sh_max = nan(1, length(prcnts));
    
    for i = 1:length(prcnts)
        fprintf('  Percentile %d%% (%d/%d)...', prcnts(i), i, length(prcnts));
        
        ids = find(abs(zc_wL) < vals(i));
        dFF_L = dff_rz(ids, L_state==1);
        dFF_Q = dff_rz(ids, L_state==0);
        coeff_L = pca(dFF_L');
        coeff_Q = pca(dFF_Q');
        A_metric = eye(size(coeff_Q,1));
        [theta, ~, ~] = subspacea(coeff_Q(:,1:num_PCs), coeff_L(:,1:num_PCs), A_metric);
        
        angle_LQ_min_current(i) = min(theta);
        angle_LQ_max_current(i) = max(theta);
        
        ashes_min = [];
        ashes_max = [];
        N_itr = 20;
        
        for kk = 1:N_itr
            bl_sz = 100;
            dff_sh_bl = dff_rz(:, randblock(1:floor(T/bl_sz)*bl_sz, [1, bl_sz]));
            dff_1 = dff_sh_bl(ids, 1:size(dff_sh_bl,2)/2);
            dff_2 = dff_sh_bl(ids, size(dff_sh_bl,2)/2+1:end);
            coeff_1 = pca(dff_1');
            coeff_2 = pca(dff_2');
            [theta_sh, ~, ~] = subspacea(coeff_1(:,1:num_PCs), coeff_2(:,1:num_PCs), A_metric);
            
            ashes_min = [ashes_min, min(theta_sh)];
            ashes_max = [ashes_max, max(theta_sh)];
        end
        
        angle_LQ_sh_min(i) = nanmean(ashes_min);
        angle_LQ_sh_max(i) = nanmean(ashes_max);
        
        fprintf(' Done\n');
    end
    
    angles_LQ_min = [angles_LQ_min; angle_LQ_min_current];
    angles_LQ_max = [angles_LQ_max; angle_LQ_max_current];
    angles_LQ_shuffle_min = [angles_LQ_shuffle_min; angle_LQ_sh_min];
    angles_LQ_shuffle_max = [angles_LQ_shuffle_max; angle_LQ_sh_max];
    
    fprintf('File %d/%d completed\n', file_i, length(file_list));
end

fprintf('\n=== All files processed ===\n');
fprintf('Total valid files: %d\n', size(angles_LQ_min, 1));

figure('Position', [100,100,400,300])
hold on

meanangle_min = nanmean(angles_LQ_min, 1);
meanangle_max = nanmean(angles_LQ_max, 1);
meanshuffle_min = nanmean(angles_LQ_shuffle_min, 1);
meanshuffle_max = nanmean(angles_LQ_shuffle_max, 1);

SEM_min = nanstd(angles_LQ_min, 0, 1) / sqrt(size(angles_LQ_min, 1));
SEM_max = nanstd(angles_LQ_max, 0, 1) / sqrt(size(angles_LQ_max, 1));
SEM_sh_min = nanstd(angles_LQ_shuffle_min, 0, 1) / sqrt(size(angles_LQ_shuffle_min, 1));
SEM_sh_max = nanstd(angles_LQ_shuffle_max, 0, 1) / sqrt(size(angles_LQ_shuffle_max, 1));

color_max = [0, 0.5, 0.4];
color_max_light = [0.6, 0.8, 0.75];
color_min = [0, 0, 0];
color_min_light = [0.7, 0.7, 0.7];

h1 = errorbar(prcnts, meanangle_max, SEM_max, '-', 'LineWidth', 2, 'Color', color_max);
h2 = errorbar(prcnts, meanshuffle_max, SEM_sh_max, '-', 'LineWidth', 2, 'Color', color_max_light);
h3 = errorbar(prcnts, meanangle_min, SEM_min, 'k-', 'LineWidth', 2, 'Color', color_min);
h4 = errorbar(prcnts, meanshuffle_min, SEM_sh_min, '-', 'LineWidth', 2, 'Color', color_min_light);

legend([h1, h2, h3, h4], {'Max angle', 'Max shuffle', 'Min angle', 'Min shuffle'}, ...
    'Location', 'Best')
legend boxoff

ylim([0, 1.7])
fs = 15;
set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [.025, .025])
ylabel('Angle (rad.)')
xlabel('MFA+/MFA- excl. percentile (%)')

fileName = 'manifold_angle_percentiles_combined';
fullFilePathPDF = fullfile(fP, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
fprintf('Saved: %s\n', fullFilePathPDF);

fprintf('\n=== Analysis complete ===\n');

%%
datasavep = 'X:\MFB\Processed\Figure 3';
if ~exist(datasavep, 'dir')
    mkdir(datasavep);
end
fileName = 'manifold_angle_percentiles.mat';
fullFilePath = fullfile(datasavep, fileName);

save(fullFilePath, ...
    'angles_LQ_min', 'angles_LQ_max', ...           
    'angles_LQ_shuffle_min', 'angles_LQ_shuffle_max', ...
    'prcnts', ...                               
    'num_PCs', ...                           
    'Ns', ...                                      
    'file_list');                 

fprintf('Data saved to: %s\n', fullFilePath);

%% 3fg pm/nm exclude only
num_PCs = 3;
Ns = [];
prcnts = 10:10:90;
sig_level = 2;


color_max = [0, 0.5, 0.4]; 
color_max_light = [0.6, 0.8, 0.75];
color_min = [0, 0, 0];
color_min_light = [0.7, 0.7, 0.7];

for md = 1:2
    angles_LQ_min = [];
    angles_LQ_max = [];
    angles_LQ_shuffle_min = [];
    angles_LQ_shuffle_max = [];
    
    fprintf('\n=== Processing mode %d ===\n', md);
    
    good_sessions = [1,2,4,6,8:16];
    for file_i = 1:length(good_sessions)
        fprintf('\n=== Processing file %d/%d ===\n', file_i, length(good_sessions));
        
        current_file = fullfile(datapath, file_list(file_i).name);
        fprintf('File: %s\n', file_list(file_i).name);
        
        data = load(current_file);
        dff_rz = zscore(data.mf.c.dff_r')';
        L_state = data.mf.b.L_state;
        dt = 10;
        [N, T] = size(dff_rz);
        Ns = [Ns, N];


        pm_ids = data.mf.a.pm_ids;
        nm_ids = data.mf.a.nm_ids;
        zc_wL = data.mf.a.zc_wL;

        vals_pm = prctile(zc_wL(pm_ids), 100 - prcnts);
        vals_nm = prctile(zc_wL(nm_ids), 100 - prcnts);
        
        angle_LQ_min_val = nan(1, length(prcnts));
        angle_LQ_max_val = nan(1, length(prcnts));
        angle_LQ_sh_min = nan(1, length(prcnts));
        angle_LQ_sh_max = nan(1, length(prcnts));
        
        for i = 1:length(prcnts)
            fprintf('  Percentile %d%% (%d/%d)...', prcnts(i), i, length(prcnts));
            
            if md == 1
                suffix = 'PM';
                exclude_ids = pm_ids(abs(zc_wL(pm_ids)) >= abs(vals_pm(i)));
                keep_mask = true(N, 1);
                keep_mask(exclude_ids) = false;
                ids = find(keep_mask);
                zz = dff_rz;
                
            elseif md == 2
                suffix = 'NM';
                exclude_ids = nm_ids(abs(zc_wL(nm_ids)) >= abs(vals_nm(i)));
                keep_mask = true(N, 1);
                keep_mask(exclude_ids) = false;
                ids = find(keep_mask);
                zz = dff_rz;
            end


            if length(ids) < num_PCs
                fprintf(' Skipped (N=%d < num_PCs=%d)\n', length(ids), num_PCs);
                continue;
            end
            
            dFF_L = zz(ids, L_state == 1);
            dFF_Q = zz(ids, L_state == 0);
            coeff_L = pca(dFF_L');
            coeff_Q = pca(dFF_Q');
            
            if sum(all(coeff_Q == 1)) == 1 || sum(all(coeff_L == 1)) == 1
                fprintf(' Skipped (PCA issue)\n');
                break;
            end
            
            A_metric = eye(size(coeff_Q, 1));
            [theta, ~, ~] = subspacea(coeff_Q(:,1:num_PCs), coeff_L(:,1:num_PCs), A_metric);
            
            angle_LQ_min_val(i) = min(theta);
            angle_LQ_max_val(i) = max(theta);
            
            ashes_min = [];
            ashes_max = [];
            N_itr = 20;
            
            for kk = 1:N_itr
                bl_sz = 100;
                dff_sh_bl = zz(:, randblock(1:floor(T/bl_sz)*bl_sz, [1, bl_sz]));
                dff_1 = dff_sh_bl(ids, 1:size(dff_sh_bl,2)/2);
                dff_2 = dff_sh_bl(ids, size(dff_sh_bl,2)/2+1:end);
                coeff_1 = pca(dff_1');
                coeff_2 = pca(dff_2');
                A_metric = eye(size(coeff_1, 1));
                [theta_sh, ~, ~] = subspacea(coeff_1(:,1:num_PCs), coeff_2(:,1:num_PCs), A_metric);
                
                ashes_min = [ashes_min, min(theta_sh)];
                ashes_max = [ashes_max, max(theta_sh)];
            end
            
            angle_LQ_sh_min(i) = nanmean(ashes_min);
            angle_LQ_sh_max(i) = nanmean(ashes_max);
            
            fprintf(' Done\n');
        end
        
        angles_LQ_min = [angles_LQ_min; angle_LQ_min_val];
        angles_LQ_max = [angles_LQ_max; angle_LQ_max_val];
        angles_LQ_shuffle_min = [angles_LQ_shuffle_min; angle_LQ_sh_min];
        angles_LQ_shuffle_max = [angles_LQ_shuffle_max; angle_LQ_sh_max];
    end
    

    figure('Position', [100, 100, 400, 300]);
    hold on;
    
    meanangle_min = nanmean(angles_LQ_min, 1);
    meanangle_max = nanmean(angles_LQ_max, 1);
    meanshuffle_min = nanmean(angles_LQ_shuffle_min, 1);
    meanshuffle_max = nanmean(angles_LQ_shuffle_max, 1);
    
    SEM_min = nanstd(angles_LQ_min, 0, 1) / sqrt(size(angles_LQ_min, 1));
    SEM_max = nanstd(angles_LQ_max, 0, 1) / sqrt(size(angles_LQ_max, 1));
    SEM_sh_min = nanstd(angles_LQ_shuffle_min, 0, 1) / sqrt(size(angles_LQ_shuffle_min, 1));
    SEM_sh_max = nanstd(angles_LQ_shuffle_max, 0, 1) / sqrt(size(angles_LQ_shuffle_max, 1));
    
    h1 = errorbar(prcnts, meanangle_max, SEM_max, '-', 'LineWidth', 2, 'Color', color_max);
    h2 = errorbar(prcnts, meanshuffle_max, SEM_sh_max, '-', 'LineWidth', 2, 'Color', color_max_light);
    h3 = errorbar(prcnts, meanangle_min, SEM_min, '-', 'LineWidth', 2, 'Color', color_min);
    h4 = errorbar(prcnts, meanshuffle_min, SEM_sh_min, '-', 'LineWidth', 2, 'Color', color_min_light);
    
    legend([h1, h2, h3, h4], {'Max angle', 'Max shuffle', 'Min angle', 'Min shuffle'}, ...
        'Location', 'Best');
    legend boxoff;
    
    ylim([0, 1.7]);
    xlim([0, 100]);
    fs = 15;
    set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [0.025, 0.025]);
    ylabel('Angle (rad.)');
    xlabel([suffix, ' excl. percentile (%)']);
    
    mfbFolderPath = 'X:\MFB';
    currentDate = datestr(now, 'yyyy-mm-dd');
    fP = fullfile(mfbFolderPath, 'Figures', 'Figure3', currentDate);
    if ~exist(fP, 'dir')
        mkdir(fP);
    end
    
    fileName = ['manifold_angle_percentiles_', suffix, '_combined'];
    fullFilePathPDF = fullfile(fP, [fileName,'.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    fprintf('Saved: %s\n', fullFilePathPDF);
    
    datasavep = 'X:\MFB\Processed\Figure 3';
    if ~exist(datasavep, 'dir')
        mkdir(datasavep);
    end
    
    fileName_mat = ['manifold_angle_percentiles_', suffix];
    fullFilePath = fullfile(datasavep, [fileName_mat, '.mat']);
    save(fullFilePath, ...
        'angles_LQ_min', 'angles_LQ_max', ...
        'angles_LQ_shuffle_min', 'angles_LQ_shuffle_max', ...
        'prcnts', 'num_PCs', 'Ns');
    fprintf('Data saved: %s\n', fullFilePath);
end

fprintf('\n=== All analysis complete ===\n');

