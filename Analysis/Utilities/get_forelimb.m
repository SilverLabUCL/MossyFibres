% load dlc csv file
clear file_list folder_path dir
folder_path = 'C:\Users\Yizhou Xie\Downloads\all results\forelimb';
file_list = dir(fullfile(folder_path, '*.csv'));
allfilesinfo = 'C:\Users\Yizhou Xie\Downloads\all time files\Exp_List.xlsx';
T1 = readtable(allfilesinfo, 'Range', 'A2');
vid_filename = T1{find(strcmp(file, T1{:, 1}), 1),2};
vid_filename = vid_filename{1};
all_names = {file_list.name};
match_idx = find(contains(all_names, vid_filename), 1);
dlc_file_path = fullfile(folder_path, file_list(match_idx).name);
disp(dlc_file_path)
zd2 = importdata(dlc_file_path);
xy_L = zd2.data(:,2:3);
conf_L0 = zd2.data(:,4);
dis_L0 = fastsmooth(sqrt(nansum(diff(xy_L).^2,2)), 5,3,1);
xy_R = zd2.data(:,5:6);
conf_R0 = zd2.data(:,7);
dis_R0 = fastsmooth(sqrt(nansum(diff(xy_R).^2,2)), 5,3,1);

% find video name and micro name
tokens = regexp(dlc_file_path, '(\d{6}_\d{2}_\d{2}_\d{2})_VidRec', 'tokens');
if ~isempty(tokens)
    vid_name = tokens{1}{1};
end
allfilesinfo = 'C:\Users\Yizhou Xie\Downloads\all time files\Exp_List.xlsx';
T1 = readtable(allfilesinfo, 'Range', 'A2');
idx = strcmp(T1.Var2, vid_name);
micro_name = T1.Var1{idx};
eye_rel = string(['C:\Users\Yizhou Xie\Downloads\all time files\' vid_name ' VidRec\EyeCam-relative times.txt']);
opts = detectImportOptions(eye_rel, 'Delimiter', '\t', 'VariableNamesLine', 0);
TFL = readtable(eye_rel, opts);
FL_time_vals = TFL{:,2};

% trim
FL_time_vals = FL_time_vals(1:length(dis_L0));
if strcmp(micro_name,'200130_13_36_14 FunctAcq')==1
    FL_time_vals = TFL{:,2};
    FL_time_vals = FL_time_vals(122:121+length(dis_L0));
end

% get tr and tt
tt = (FL_time_vals - FL_time_vals(1) - Tshift)/1e3;
dis_L = interp1(tt, dis_L0, tr, 'linear', 'extrap');
dis_R = interp1(tt, dis_R0, tr, 'linear', 'extrap');

% process xy with confidence threshold
xy_L = zd2.data(:,2:3);
conf_L = zd2.data(:,4);
xy_R = zd2.data(:,5:6);
conf_R = zd2.data(:,7);
xy_L(conf_L < 0.8, :) = NaN;
xy_R(conf_R < 0.8, :) = NaN;
xy_L = fillmissing(xy_L, 'linear', 'MaxGap', 5);
xy_R = fillmissing(xy_R, 'linear', 'MaxGap', 5);

% compute raw speed (before smoothing)
dis_L2_raw = sqrt(nansum(diff(xy_L).^2, 2));
dis_R2_raw = sqrt(nansum(diff(xy_R).^2, 2));

% compute acceleration from raw speed (before smoothing)
vid_fs = 100; % replace with your actual video frame rate
acc_L_raw = [0; diff(dis_L2_raw)] * vid_fs;
acc_R_raw = [0; diff(dis_R2_raw)] * vid_fs;

% smooth acceleration (lighter than speed)
acc_L2 = fastsmooth(acc_L_raw, 5, 3, 1);
acc_R2 = fastsmooth(acc_R_raw, 5, 3, 1);

% smooth speed (as before)
dis_L2 = fastsmooth(dis_L2_raw, 10, 3, 1);
dis_R2 = fastsmooth(dis_R2_raw, 10, 3, 1);

% interpolate everything to neural signal time axis
dis_L2  = interp1(tt, dis_L2,  tr, 'linear', 'extrap');
dis_R2  = interp1(tt, dis_R2,  tr, 'linear', 'extrap');
acc_L2  = interp1(tt, acc_L2,  tr, 'linear', 'extrap');
acc_R2  = interp1(tt, acc_R2,  tr, 'linear', 'extrap');