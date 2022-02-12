% This program plots the data for Rat 980 (originally 913)
% The data include the occupancy and histogram of the spikes and the
% sequential firing of multiple spikes before jumping
%
%   See also PLOTRATEMAP.
%       
% SGL 2022-01-26 (originally 2021-01-31)
%
clc; clear; close all;
[datafile,exp_directory] = uigetfile(fullfile('D:\Analysis', 'data.mat'), 'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename,'pos','cluster','exp','ppcm', 'daq');

[img_file,img_directory] = uigetfile(fullfile(exp_directory,'frames','*.pgm'), 'Select Image File');
if isequal(img_file, 0)
    error('Image file was not selected!')
end
img_filename = fullfile(img_directory, img_file);

% [csc_filename,csc_directory] = uigetfile(fullfile(exp_directory,'*.ncs'), 'Select a CSC File');
% if isequal(csc_filename, 0)
%     error('CSC file was not selected!')
% end
% csc_filename = fullfile(csc_directory, csc_filename);

colors = ["#EDB120" "#7E2F8E" "yellow" "#A2142F" "red" "magenta" "green" "#D95319"];
colors = repmat(colors', ceil(length(cluster)/length(colors))); % repeat the colors to match the total number of clusters
%colors(35) = "#00CCCC";

img = imread(img_filename);
xmax = ceil(size(img,2)/ppcm);

start = tic;

%% Exclusions 
% MODIFY FOR EACH DAY!!!
% 
cluster_exlude = [];

% % Excluding some clusters (e.g. inter-neurons) -- not a good idea
% if ~isempty(cluster_exlude)
%     cluster(cluster_exlude)=[];
%     disp(['clusters ' num2str(cluster_exlude) ' excluded.']); 
% end

%% Plotting the position vs time and marking neural recording times,
% maze times and lap detection times

figure(1)
[lap, idx_analysis] = lapdetector(exp_directory);
set(gcf, 'Position', [100 100 1536 1175]);
xlabel('Time (sec)')
ylabel('Horizontal position (cm)')
title('Overall view')
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'overall.jpg'])

%% Enhance pos and lap
pos.lap = zeros(size(pos.t));
pos.dir = strings(size(pos.t));
pos.status = strings(size(pos.t));
for l = 1:length(lap)
    % Adding laps and direction to pos
    idx = (pos.t >=lap(l).t(1)) & (pos.t <lap(l).t(2)); % lap
    pos.lap(idx) = lap(l).no;
    pos.dir(idx)= lap(l).dir;
    pos.status(idx)= lap(l).status;
end

for l = 1:length(lap)
    % Adding gap size to each lap
    lap(l).gap = 888.651/ppcm + [0 pos.len(pos.t == lap(l).t_jump)];  % two edges of the gap (change for each setup!!!)
    lap(l).gap_length = diff(lap(l).gap);
    % time of crossing the gap edges:
    idx = pos.lap==l;
    gapidx = pos.x(idx)>=lap(l).gap(1) & pos.x(idx)<lap(l).gap(2);
    gapidx = [0; diff(gapidx)];
    post = pos.t(idx);
    lap(l).t_cross = post(gapidx~=0)';
end 

for l = 1:length(lap)
    lap(l).corr = max([lap.gap_length]) - lap(l).gap_length;
    lap(l).gap_corr = lap(l).gap + lap(l).corr;
end

%% Adding laps and direction to clusters
for c=1:length(cluster)
    cluster(c).lap = zeros(size(cluster(c).t));
    cluster(c).dir = strings(size(cluster(c).t));
    cluster(c).status = strings(size(cluster(c).t));
    cluster(c).x_corr = cluster(c).x;
    for l = 1:length(lap)
        idx = (cluster(c).t >=lap(l).t(1)) & (cluster(c).t <lap(l).t(2)); % lap
        cluster(c).lap(idx)= lap(l).no;
        cluster(c).dir(idx)=lap(l).dir;
        cluster(c).status(idx)=lap(l).status;
        cluster(c).x_corr(idx) = cluster(c).x(idx) + lap(l).corr; 
    end
end

%% Plotting time vs horizontal position for all laps
figure(2)
N = max(nnz([lap.dir]=="left"),nnz([lap.dir]=="right"));

% rightward laps
i = 2;
for l=[lap([lap.dir]=="right").no]
    subplot(N,2,i);
    plot(pos.x(pos.lap == l), pos.t(pos.lap == l), 'b')
    hold on
    xlim([0 xmax])
    ylim(lap(l).t)
    title(['rightward lap no ' num2str(lap(l).no)])
    i = i + 2;
end

% leftward laps
i = 1;
for l=[lap([lap.dir]=="left").no]
    subplot(N,2,i);
    plot(pos.x(pos.lap == l), pos.t(pos.lap == l), 'b')
    hold on
    xlim([0 xmax])
    ylim(lap(l).t)
    title(['leftward lap, lap no ' num2str(lap(l).no)])
    i = i + 2;
end

%% interpolation at each lap
posi.dt = 1/1000; % interpolation 1 kHz
posi.t = [];
posi.x = [];
posi.x_corr = []; % corrected for the gap size
posi.y = [];
posi.vx = [];
posi.vy = [];
posi.lap = [];
posi.dir = [];
posi.status = [];
posi.filt.vx = [];
posi.filt.vy = [];
posi.filt.s = [];
posi.filt.ax = [];
posi.filt.ay = [];
% angles
posi.roll = [];
posi.pitch = [];
posi.yaw = [];
posi.filt.avx = [];
posi.filt.avy = [];
posi.filt.avz = [];
posi.filt.av = [];

for l = 1:length(lap)
    % interpolation
    interp.t = lap(l).t(1):posi.dt:lap(l).t(2);
    % interpolation for position
    interp.x = interp1(pos.t, pos.x, interp.t);
    interp.x_corr = interp.x + lap(l).corr;
    interp.y = interp1(pos.t, pos.y, interp.t);
    % interpolation for velocity
    interp.vx = interp1(pos.t, pos.vx, interp.t);
    interp.vy = interp1(pos.t, pos.vy, interp.t); 
    % angles
    interp.roll = interp1(pos.t, pos.roll, interp.t);
    interp.pitch = interp1(pos.t, pos.pitch, interp.t);
    interp.yaw = interp1(pos.t, pos.yaw, interp.t);
    
    interp.lap = zeros(size(interp.t));
    interp.dir = strings(size(interp.t));
    interp.status = strings(size(interp.t));
    interp.lap(:)=lap(l).no;
    interp.dir(:)=lap(l).dir;
    interp.status(:)=lap(l).status;
    
    posi.t = [posi.t interp.t];
    posi.x = [posi.x interp.x];
    posi.x_corr = [posi.x_corr interp.x_corr];
    posi.y = [posi.y interp.y];
    posi.vx = [posi.vx interp.vx];
    posi.vy = [posi.vy interp.vy];
    posi.dir = [posi.dir interp.dir];
    posi.status = [posi.status interp.status];
    posi.lap = [posi.lap interp.lap];
    % angles
    posi.roll = [posi.roll interp.roll];
    posi.pitch = [posi.pitch interp.pitch];
    posi.yaw = [posi.yaw interp.yaw];    

    % filtered data:
    interp.filt.vx = filtertheta(interp.t, interp.vx, 0.01, 10); % cm/sec
    interp.filt.vy = filtertheta(interp.t, interp.vy, 0.01, 10); % cm/sec
    interp.filt.s = vecnorm([interp.filt.vx interp.filt.vy]')';  % speed in cm/sec
    interp.filt.ax = gradient(interp.filt.vx)./gradient(interp.t);   % ax in cm/sec
    interp.filt.ay = gradient(interp.filt.vy)./gradient(interp.t);   % ax in cm/sec
    % angular velocities
    interp.avx = gradient(interp.roll)./gradient(interp.t);   % avx in deg/sec
    interp.avy = gradient(interp.pitch)./gradient(interp.t);   % avy in deg/sec
    interp.avz = gradient(interp.yaw)./gradient(interp.t);   % avz in deg/sec
    interp.filt.avx = filtertheta(interp.t, interp.avx, 0.01, 10); % deg/sec
    interp.filt.avy = filtertheta(interp.t, interp.avy, 0.01, 10); % deg/sec
    interp.filt.avz = filtertheta(interp.t, interp.avz, 0.01, 10); % deg/sec
    
    posi.filt.vx = [posi.filt.vx interp.filt.vx];
    posi.filt.vy = [posi.filt.vy interp.filt.vy];
    posi.filt.s = [posi.filt.s interp.filt.s];
    posi.filt.ax = [posi.filt.ax interp.filt.ax];
    posi.filt.ay = [posi.filt.ay interp.filt.ay];
    % angle
    posi.filt.avx = [posi.filt.avx interp.filt.avx];
    posi.filt.avy = [posi.filt.avy interp.filt.avy];
    posi.filt.avz = [posi.filt.avz interp.filt.avz];
end
posi.filt.av = vecnorm([posi.filt.avx; posi.filt.avy; posi.filt.avz]);  % angular speed in deg/sec
%% plotting the occupancy
figure(4); clf;
img = imread(img_filename);
img = imlocalbrighten(img);
img = medfilt2(img);
imshow(img);
figure(4)
hold on
h = quiver(posi.x * ppcm,posi.y * ppcm,posi.vx * ppcm,posi.vy * ppcm,2);
plot(repmat(lap(1).gap(1),2,1)* ppcm, ylim* ppcm,'y','LineWidth',2);
plot(repmat(lap(1).gap(2),2,1)* ppcm, ylim* ppcm,'y','LineWidth',2);
title('Occupancy with Velocity')
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'occupancy_with_velocity2.jpg'])

%% Rate map binning
hist.bin_size = 3; % cm
hist.edges = 0:hist.bin_size:xmax; % size of image to cm
hist.posi = histcounts(posi.x,hist.edges) * posi.dt; % seconds in each bin

%% detect the more exact time of jump
for l =1:length(lap)
    idx = (posi.t>lap(l).t_jump-0.2) & (posi.t<lap(l).t_jump+0.2) & (abs(posi.filt.vx)>30);
    idx = find(idx,1);
    lap(l).t_jump_exact = posi.t(idx);
end

%% CSC
% theta phase
start2 = tic;
for c=1:length(cluster)
    cluster(c).phase = nan*cluster(c).t;
end

for sh=1:4
    csc_filename=fullfile(exp_directory, ['CA1-Shank' num2str(sh)],'B','LFP32.ncs'); % use plotcsc to optimize it
    [time,data] = read_bin_csc(csc_filename);
    for l=1:length(lap)
        idx = time >= lap(l).t(1) & time <= lap(l).t(2);
        timecsc = time(idx);
        lfp = data(idx);
        [theta, phase] = filtertheta(timecsc,lfp);
        
        % looking at all the clusters in the same shank
        for c=[cluster([cluster.sh]==sh).no]
            idx = [cluster(c).lap]==l;
            if nnz(idx) > 0 % if there a firing for cluster c in this lap
                cluster(c).phase(idx) = interp1(timecsc,phase, cluster(c).t(idx));
            end
        end
    end
end
toc(start2)
%% Saving data
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
save(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'colors','xmax','daq','hist');
fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);
disp(['File ' mat_filename ' has been created!'])