%% Extract data from the output of winclust (cluster mazes) and 
% tracking.dat (tracking data) and save it in data.mat file 

%% Neural data
exp_directory = '~/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02';
exp_directory = uigetdir(exp_directory,'Select Experiment Folder');
if exp_directory == 0
    return;
end

listing = dir(fullfile(exp_directory,'Neuralynx','TT*','cl-maze*.*'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing); % number of maze-clusters

tt_no = str2double(extractAfter(folders,'TT'));
maze_no = floor(str2double(extractAfter(names,'cl-maze')));
cluster_no = str2double(extractAfter(names,'.'));
cluster_no(isnan(cluster_no))=0;    % cluster 0

mat_filename = fullfile(exp_directory,'data.mat');

absolue_paths = mat2cell(fullfile(folders,names),ones(N,1));
A = cellfun(@(x) importdata(x,',',13), absolue_paths, 'UniformOutput', false);

%% Postion data
addpath('../tracking');
[t, x, y, p, frame] = readtrackingdata(exp_directory);
offset = 48.5827 - 0.2; %%%% ONLY for DAY 2 !! %%%%%% after modification
   
% velocity
% 2020-03 5 ft = 1524 mm = 480 pixels (each pixel = 3.175 mm)
% 2020-10 3 ft = 914 mm = 840 pixels = norm([296 372]-[1136 348],2) 
% each pixel ~ 1.1 mm
ppcm = norm([296 372]-[1136 348])/91.4; % pixels per cm

pos.t = t + offset; % (+offset so Nlx time is the origin)
pos.x = x / ppcm; % cm
pos.y = y / ppcm; % cm

% velocity
pos.vx = gradient(pos.x)./gradient(pos.t); % Vx in cm/sec
pos.vy = gradient(pos.y)./gradient(pos.t); % Vy in cm/sec
pos.s = vecnorm([pos.vx pos.vy]')'; % speed in cm/sec
%pos.hd = atan2d(pos.vy,pos.vx); % estimation of hd based

% IO port status p and frame number
pos.p = p;
pos.frame = frame;
        
%% spike data
cluster(N).name ='';
for index = 1:N
    cluster(index).name = ['tt' num2str(tt_no(index)) '_m' num2str(maze_no(index)) '_c' num2str(cluster_no(index))];
    cluster(index).tt = tt_no(index); % [spike.tt] == 7
    cluster(index).m = maze_no(index); % [cluster([spike.m] == 1).no]
    cluster(index).cl = cluster_no(index);
    cluster(index).no = index;
    cluster(index).ti = str2double(A{index}.textdata{12})*1e-6; % sec
    cluster(index).tf = str2double(A{index}.textdata{13})*1e-6; % esc
    cluster(index).t = (A{index}.data(:,18))*1e-6; % sec 
    % interpolation for position
    cluster(index).x = interp1(pos.t, pos.x, cluster(index).t);
    cluster(index).y = interp1(pos.t, pos.y, cluster(index).t);
    % interpolation for velocity
    cluster(index).vx = interp1(pos.t, pos.vx, cluster(index).t);
    cluster(index).vy = interp1(pos.t, pos.vy, cluster(index).t);
    cluster(index).s = vecnorm([cluster(index).vx cluster(index).vy]')';
end

spike = cluster;
save(mat_filename,'pos','cluster','spike','ppcm', 'offset');
disp(['File ' mat_filename ' has been created!'])