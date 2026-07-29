clear all;close all;clc
%% Select a dataset
file = '191209_13_44_12';
quickAnalysis;
%% fig.1d SpatialPosition3d

MF_id = [];
figure('Position',[360, 278, 800, 800]*.6)
hold on; rotate3d on; grid on;

xyz_mx = [250,250, 0];
xyz_mn = [0,0, -100];
scatter3(Xall(incl_ids), Yall(incl_ids), Zall(incl_ids), 30, Zall(incl_ids), 'LineWidth', 2);
grayLevels = linspace(0.1, 0.8, 256);
customGray = repmat(grayLevels', 1, 3);
colormap([0.5 0.5 0.5]);
xlim([0,xyz_mx(1)]);
ylim([0,xyz_mx(2)]);
zlim([xyz_mn(3),xyz_mx(3)]);
xlabel('X (μm)'); ylabel('Y (μm)'); zlabel('Z (μm)');
view(-28,28);

set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'   , ...
    'ZMinorTick'  , 'off'  )

Date = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure1', Date);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end

fileName = ['SpatialPosition3d_' file];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% fig.1f subgroups correlation
cc0 = corr(dff0_r', 'rows', 'complete');

% choose # subgroups
subgroups_n = 4;
ssgroups = subgroups(randperm(length(subgroups), subgroups_n));
sids = cellfun(@(sg) sg(randperm(numel(sg), min(2, numel(sg)))), ssgroups, 'UniformOutput', false);
allMembers = horzcat(sids{:});
ccS = cc0(allMembers, allMembers);

% color
nColors = 256;
halfN = nColors / 2;
r = [linspace(0, 1, halfN), linspace(1, 1, halfN)];
g = [linspace(0, 1, halfN), linspace(1, 0, halfN)];
b = [linspace(1, 1, halfN), linspace(1, 0, halfN)];
colorMap = [r' g' b'];

fig = figure;
set(fig, 'Position', [400, 100, 410, 380]);

imagesc(ccS);
grid on;
set(gca, 'GridLineStyle', '-', 'GridColor', 'k', 'GridAlpha', 0.2, 'LineWidth', 1.5);
colormap(colorMap);
caxis([-1 1]);
axis off;

axPos = get(gca, 'Position');

cb = colorbar('Position', [axPos(1), axPos(2)- 0.04, axPos(3)*0.3, 0.03], 'Orientation', 'horizontal');
set(cb, 'Ticks', [-1, 0, 1],'FontSize',14);
cbPos = get(cb, 'Position');
cbRightEdge = cbPos(1) + cbPos(3);
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
txt = text(cbRightEdge + max(XLim)*0.34, max(YLim)*1.047, 'Correlation', 'VerticalAlignment', 'middle','FontSize',18);

xticks([]);
yticks([]);

numSubgroups = numel(sids);
numColors = size(my_colors, 1);
rng shuffle
randIndices = randperm(numColors);
colors = [
       127, 88, 175;
       100, 197, 235;
       232, 77, 138;
       254, 179, 38]/255;

totalElements = sum(cellfun(@(x) numel(x), sids));

% Calculate the boundary indices for each subgroup (x,y)
x_Sum = min(YLim)+[0, cumsum(cellfun(@length, sids))];
totalLength = max(x_Sum);
groupItv = diff(x_Sum);
y_Sum = zeros(size(x_Sum));
currentStart = totalLength;
for i = 1:length(groupItv)
    y_Sum(i) = currentStart;
    currentStart = currentStart - groupItv(i);
end
y_Sum(end) = 0.5;

axPos = get(gca, 'Position');

% left side color
leftAx = axes('Position', [axPos(1)-0.04, axPos(2), 0.02, axPos(4)], 'Units', 'normalized');
set(leftAx, 'YLim', [0, size(ccS, 1)]+0.5, 'XLim', [0, 1], 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
axes(leftAx);
for i = length(ssgroups):-1:1
    startIdx = y_Sum(i);
    endIdx = y_Sum(i+1);
    patch([0, 1, 1, 0], [startIdx, startIdx, endIdx, endIdx], colors(i,:), 'EdgeColor', 'none');
    midLineY = (startIdx + endIdx) / 2;
    line([0, 1], [midLineY, midLineY], 'Color', 'w', 'LineWidth', 1);
end

topAx = axes('Position', [axPos(1), axPos(2)+axPos(4)+0.02, axPos(3), 0.02], 'Units', 'normalized');
set(topAx, 'XLim', [0, size(ccS, 2)]+0.5, 'YLim', [0, 1], 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
axes(topAx);

for i = 1:length(ssgroups)
    startIdx = x_Sum(i);
    endIdx = x_Sum(i+1);
    
    patch([startIdx, endIdx, endIdx, startIdx], [0, 0, 1, 1], colors(i,:), 'EdgeColor', 'none');
    midLineX =(startIdx + endIdx) / 2;
    line([midLineX, midLineX], [0, 1], 'Color', 'w', 'LineWidth', 1);
end

fileName = ['subgroups_' file '_corr'];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% - fig.1e sample df/f
plot_sampledff = 1;
plot_location3D = 1;

MF_id = cell2mat(cellfun(@(x) x(1:2)', ssgroups, 'UniformOutput', false));
MF_id = MF_id(:);
MF_no = length(MF_id);

if plot_sampledff

    h = figure('Position', [300, 100, 850, 520], 'visible', 'on'); hold on;
    pre = 15;
    scale_factor = 6;
    y_shift = 0;
    dff_all = permute(ROI_dff_all,[2,1,3]);
    dff_all_z = zeros(size(dff_all));

    for i = 1:size(dff_all, 2)
        X = dff_all(:, i, :);
        X_m = mean(X, 'all', 'omitnan');
        X_s = std(X, 0, 'all', 'omitnan');
        dff_all_z(:, i, :) = (X - X_m) / X_s;
    end

    % assign color
    colorAssign = [];
    ci = 1;
    for i = 1:length(ssgroups)
        csize = min(2, numel(ssgroups{i}));
        colorAssign(ci:ci+csize-1, :) = ...
        repmat(colors(i, :), csize, 1);
        ci = ci + csize;
    end


    for i = 1:MF_no
        zzf = squeeze(dff_all_z(:,MF_id(i),pre:end)); % Z-scored data
        y_shift = y_shift + nanmax(zzf(:));
        for j = 1:Ntr
            if j == 1; pre = 30;
            else; pre = 1; end
            
            ext_shift = 0;
            if i == 1; ext_shift = -1; end
            plot(((j-1)*T(end)+T(pre:end))/1e3, (i-1)*scale_factor+ext_shift + squeeze(dff_all_z(j,MF_id(i),pre:end)), ...
                'linewidth', 1.5, 'color', colorAssign(i,:));
        end
    end
    
    for j = 1:Ntr
        tpuff = ((j-1)*T(end)+5e3)/1e3;
        yTriTop = MF_no * scale_factor; %
        plot(tpuff, yTriTop, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.7,0.7,0.7], 'MarkerFaceColor', [0.7,0.7,0.7]);
    end

    xlim([0-30, T(end)*Ntr/1e3+20]);
    ylim([0-2, MF_no*scale_factor+1]);
    
    x0 = T(end)*Ntr/1e3-10;
    y0 = -5;
    
    plot(x0+[0,20], [y0,y0],'k-', 'LineWidth',3)
    text(x0+4, y0-1.5,'20 s', 'FontSize',15)
    
    plot([x0+19.5,x0+19.5], [y0,y0+2],'k-', 'LineWidth',3)
    text(x0+22,y0+.5,'2 \DeltaF/F', 'FontSize',15)
    
    xlabel('Time (s)');
    ylabel('MF');
    
    set(gca, 'LineWidth', 1, 'FontSize', 12, 'ycolor','none', 'xcolor', 'none');
    axis tight;
    
    fileName = ['sample_traces_' file];
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

end


%% fig.1g

plot_1g=1;

if plot_1g==1
    folder_names = {
          '171212_16_19_37';
          '191209_14_18_13';
          '200130_13_36_14 FunctAcq';
          '191018_13_39_41';
        };
    
    all_group_ids={};
    NN_dist = [];
    
    for ijk =1:length(folder_names)
    
        file = folder_names{ijk};
        quickAnalysis;
        dt = 10;
        Nmf0 = size(dff0_r,1);
        Nmf1 = length(group_ids);
        all_group_ids = [all_group_ids, group_ids];
        
        dd_all = pdist2(xyz(:,:), xyz(:,:));
        dd_all(eye(size(dd_all))==1)=nan;
    
        % NN distance (local)
        for i = 1:Nmf1
            mids = group_ids{i};
            k_oths = mids;
            if length(mids)>1
                k = mids(1);
                for j = 1:length(mids)
                    k_oths = setdiff(k_oths,k);
                    [a,b] = min(dd_all(k,k_oths));
                    NN_dist = [NN_dist,a];
                    k = k_oths(b);
                end
            end
        end
        fprintf("NN_dist size = %d \n",size(NN_dist,2))
    end
    
    lengths = cellfun(@length, all_group_ids);
    figure(Position = [100 100 300 350]);
    histogram(lengths, 'BinEdges', 0.5:max(lengths) + 0.5, 'FaceColor', [1, 0, 0],'EdgeColor', 'black', 'LineWidth', 1);
    xlabel('MFB number per axon');
    ylabel('Count');
    ylim([0 1000])
    set(gca, 'YTick', [0 10 100 1000],'Yscale','log');
    xticks([1,6,11]);
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, ...
        'FontName'   , 'Helvetica', ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'ZMinorTick'  , 'off'  ...
        )
    
    % print
    fileName = ['MFB per axon'];
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    
    figure;
    histogram(NN_dist, 'BinWidth', 10, 'FaceColor', [1, 0, 0], 'EdgeColor', 'black', 'LineWidth', 1); % RGB color for red
    xlabel('Distance (μm)');
    ylabel('Number');
    xlim([0, ceil(max(NN_dist)/10)*10]);
    yticks([0, 50, 100])
    
    fig = gcf;
    fig.Units = 'pixels';
    fig.Position = [100 100 300 350];
    set(gca, 'LineWidth', 1, 'FontSize', 15, ...
        'FontName'   , 'Helvetica', ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'   , ...
        'ZMinorTick'  , 'off'  )
    
    hold on;
    
    mdist = nanmean(NN_dist);
    stdist = nanstd(NN_dist);
    maxdist = max(ylim);
    plot(mdist, maxdist*0.9, 'v', 'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
    
    hold on
    plot(66, maxdist*0.981, 'kv', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
    r1 = mdist - stdist;
    r2 = mdist + stdist;
    ax = gca;
    maxdist = ax.YLim(2);
    plot([r1, r2], [maxdist*0.992, maxdist*0.992], 'k', 'LineWidth', 4);
    
    fileName = ['histogram_MFA_distance_all'];
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

end


%% fig.1h plot heatmap
t_r = (1:length(spd_r)) * dt;
tb = 10000/dt; 
t_end = t_r(end);

xl1 = 0; xl2 = t_end/dt;
xl1_z = 80; xl2_z = 250;
x0_tx = -30;

sbplt1 = [1:9];
sbplt2 = [10];
sbplt3 = [11];
sbp_no = 11;    

figure('Position', [100, 20, 610, 780]);

subplot(sbp_no,1,sbplt1); hold on
bsh = randperm(Nmf);
zz = dff_rz(:,valid_t);
Loco = MI_wheel_r(valid_t);
Whisk = MI_whisker_r(valid_t);
sorted_zz = sort(zz(~isnan(zz)), 'ascend');
n = numel(sorted_zz);
v1 = round(sorted_zz(max(floor(n * 0.005), 1)), 1);
v2 = round(sorted_zz(max(floor(n * 0.995), 1)), 1);
imagesc(zz(bsh, tb:end), [v1, v2]);

for j = 1:Ntr
    tpuff = ((j-1)*T(end)+5e3)/dt;
    yTriTop = size(dff_r,1)+3;
    plot(tpuff, yTriTop, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.7,0.7,0.7], 'MarkerFaceColor', [0.7,0.7,0.7]);
end

cb = colorbar();
set(cb, 'Position', [0.06, 0.79, 0.025, 0.05],'Ticks', [v1, v2],'FontSize', 11);
titleStr = {'Z-scored', '\DeltaF/F'};
title(cb, titleStr, 'FontSize', 10, 'FontWeight', 'normal');
xlim([xl1, size(dff_rz(bsh,tb:end),2)])
xticks([])
yticks([1, Nmf])
ylim([1, Nmf + 3])
ylabel({'Mossy fibre axon number'})
set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out')
% loco traces
subplot(sbp_no,1,sbplt2); hold on;axis off
plot(fastsmooth(Loco(tb:end), 10, 3, 1), 'Color', [0,153,76]/255, 'LineWidth', 1);
plot([0,0], [0.28,0.48], 'k-', 'LineWidth',2)
annotation('textbox', [0.02 0.13 0.1 0.1], 'String', 'Loco', ...
           'FontSize', 15, 'EdgeColor', 'none');
xlim([xl1, size(dff_rz(bsh,tb:end),2)])
xticklabels([])
set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out')
% wsk traces
subplot(sbp_no,1,sbplt3);hold on;axis off
plot(fastsmooth(Whisk(tb:end),10,3,1), 'Color', [255,163,26]/255, 'LineWidth',1)
plot([0,0], [0,.1]+0.15, 'k-', 'LineWidth',2)
annotation('textbox', [0.001 0.06 0.1 0.1], 'String', 'WMI', ...
           'FontSize', 15, 'EdgeColor', 'none');
xlim([xl1, size(dff_rz(bsh,tb:end),2)])
ylim([0 0.4])
xlabel('Time (s)')
annotation('line', [0.8 0.85], [0.12 0.12], 'LineWidth', 2);
annotation('textbox', [0.79, 0.08, 0.09, 0.04], 'String', '20 s', ...
           'EdgeColor', 'none', 'FontSize', 12);

fileName = ['heatMap_' file];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');



