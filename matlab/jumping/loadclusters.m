%% Cluster data
%clear
close all

%% Neural data
exp_directory = '~/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02';
exp_directory = uigetdir(exp_directory,'Select Experiment Folder');
if exp_directory == 0
    return;
end
Nlx_directory = fullfile(exp_directory,'Neuralynx');

listing = dir(fullfile(Nlx_directory,'**','cl-maze*.*'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing); % number of maze-clusters

tt_no = str2double(extractAfter(folders,'TT'));
maze_no = floor(str2double(extractAfter(names,'cl-maze')));
cluster_no = str2double(extractAfter(names,'.'));
cluster_no(isnan(cluster_no))=0;    % cluster 0

mat_filename = fullfile(exp_directory,'data.mat');
load(mat_filename,'A')

%% Postion data
addpath('../tracking');
[t, x, y] = readtrackingdata(exp_directory);

% velocity
% 2020-03 5 ft = 1524 mm = 480 pixels (each pixel = 3.175 mm)
% 2020-10 3 ft = 914 mm = 840 pixels = norm([296 372]-[1136 348],2) 
% each pixel ~ 1.1 mm
ppcm = norm([296 372]-[1136 348])/91.4; % pixels per cm

pos.t = t;
pos.x = x / ppcm; % cm
pos.y = y / ppcm; % cm
% theta = rad2deg(atan2(pos.y-240, pos.x-320));
% pos.th = wrapTo360(theta);

% velocity
pos.vx = gradient(pos.x)./gradient(pos.t); % Vx in cm/sec
pos.vy = gradient(pos.y)./gradient(pos.t); % Vy in cm/sec
pos.s = vecnorm([pos.vx pos.vy]')'; % speed in cm/sec
%pos.hd = atan2d(pos.vy,pos.vx); % estimation of hd based 
        
%% spike data

spike(N).name ='';
tic
for index = 1:N
    spike(index).name = ['tt' num2str(tt_no(index)) '_m' num2str(maze_no(index)) '_c' num2str(cluster_no(index))];
    spike(index).tt = tt_no(index); % [spike.tt] == 7
    spike(index).m = maze_no(index); % [spike([spike.m] == 1).no]
    spike(index).cl = cluster_no(index);
    spike(index).no = index;
    spike(index).ti = str2double(A{index}.textdata{12})*1e-6;
    spike(index).tf = str2double(A{index}.textdata{13})*1e-6;
    spike(index).t = (A{index}.data(:,18))*1e-6;
    % interpolation for position
    spike(index).x = interp1(pos.t, pos.x, spike(index).t);
    spike(index).y = interp1(pos.t, pos.y, spike(index).t);
    % interpolation for velocity
    spike(index).vx = interp1(pos.t, pos.vx, spike(index).t);
    spike(index).vy = interp1(pos.t, pos.vy, spike(index).t);
    spike(index).s = vecnorm([spike(index).vx spike(index).vy]')';
end
toc

save(mat_filename,'pos','spike','-append');