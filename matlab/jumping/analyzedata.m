%%% This program plots the data for Rat 980 (originally 913)
% The data include the occupancy and histogram of the spikes and the
% sequential firing of multiple spikes before jumping
%
%   See also LAPDETECTOR, PLOTRATEMAP, PLOTTIME, PLOTHYSTERESIS, PLOTCSC.
%
% 	Date 2023-01-01 (originally 2021-01-31, 2022-01-26)
%   
%   Author Shahin G Lashkari
%
clc; clear; close all;
answer = inputdlg({'Rat', 'Date'},'Enter the rat number and the date',[1 30],{'1068', '2022-12-20'});
rat_no = answer{1};
date_str = answer{2};
%% Selecting the appropriate files
[datafile,exp_directory] = uigetfile(['E:\Rat' rat_no '\Analysis\' date_str '\raw_data.mat'],'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename,'pos','cluster','exp', 'daq');

[~, framesPath] = uigetfile(['E:\Rat' rat_no '\TrackingData\' date_str '\frames\*.pgm'],'Select a pgm File for Frame');
frames = dir(fullfile(framesPath,'*.pgm'));
img_filename = fullfile(frames(1).folder, frames(1).name);

% for phase calculations
csc_filename = fullfile(exp_directory, 'LFP','LFP.ncs'); % use apbin2lfp and plotcsc to optimize it

% Colors: https://sashamaps.net/docs/resources/20-colors/
colors = ["#e6194B" "#3cb44b" "#ffe119"  "#f58231" "#42d4f4" "#f032e6" "#fabed4" "#469990" "#9A6324"];
colors = repmat(colors', ceil(length(cluster)/length(colors))); % repeat the colors to match the total number of clusters

img = imread(img_filename);
exp.xmax = ceil(size(img,2)/exp.ppcm);
start = tic;

%% Exclusions
% MODIFY FOR EACH DAY!!!
%

cluster_exlude = find([cluster.size] >  max(std([cluster.size]),1e4)); % exclude one std dev above mean
% cluster_exlude = find(isoutlier([cluster.size]));
if exp.name == "21-Dec-2021"
    cluster_exlude = [33 48];
    colors(18)="#f58231";
    colors(28)="#3cb44b";
end

figure(10); clf
plot([cluster.no],[cluster.size], '*b'); hold on
plot(xlim, std([cluster.size])*[1 1], '--k');
plot([cluster(cluster_exlude).no],[cluster(cluster_exlude).size], 'pr');

% Excluding some clusters (e.g. inter-neurons) -- not a good idea
if ~isempty(cluster_exlude)
    cluster(cluster_exlude)=[];
    disp(['clusters ' num2str(cluster_exlude) ' excluded.']);
end

for c=1:numel(cluster)
    cluster(c).no = c;
end

%% Plotting the position vs time and marking neural recording times,
% maze times and lap detection times
disp('Detecting laps ...')
figure(1)
[lap, idx_analysis] = lapdetector(exp_directory);
set(gcf, 'Position', [100 100 1700 1175]);
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
    lap(l).gap = exp.gap_edge_in_px/exp.ppcm + [0 pos.len(pos.t == lap(l).t_jump)];  % two edges of the gap
    lap(l).gap_length = diff(lap(l).gap);
    if exp.name == "21-Dec-2021"
        lap(l).gap_depth = 2.54 * 9.125; % cm, 9.125 in
    elseif exp.name == "20-Dec-2022"
        lap(l).gap_depth = 2.54 * 10.5; % cm, 10.5 in
    elseif exp.name == "09-Nov-2022"
        lap(l).gap_depth = 2.54 * 10; % cm, 10 inch
    elseif exp.name == "16-Dec-2022"
        lap(l).gap_depth = 2.54 * 10.5; % cm, 10.5 inch    
    else
        error('specify the gap depth');
    end
    
    % time of crossing the gap edges:
    idx = pos.lap==l;
    gapidx = pos.x(idx)>=lap(l).gap(1) & pos.x(idx)<lap(l).gap(2);
    gapidx = [0; diff(gapidx)];
    post = pos.t(idx);
    lap(l).t_cross = post(gapidx~=0)';
    lap(l).t_land = lap(l).t_cross(end);
end

for l = 1:length(lap)
    lap(l).corr = max([lap.gap_length]) - lap(l).gap_length;
    %lap(l).gap_corr = lap(l).gap + lap(l).corr;
end

%% Adding laps and direction to clusters
for c=1:length(cluster)
    cluster(c).lap = zeros(size(cluster(c).t));
    cluster(c).dir = strings(size(cluster(c).t));
    cluster(c).status = strings(size(cluster(c).t));
    %cluster(c).x_corr = cluster(c).x;
    for l = 1:length(lap)
        idx = (cluster(c).t >=lap(l).t(1)) & (cluster(c).t <lap(l).t(2)); % lap
        cluster(c).lap(idx)= lap(l).no;
        cluster(c).dir(idx)=lap(l).dir;
        cluster(c).status(idx)=lap(l).status;
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
    xlim([0 exp.xmax])
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
    xlim([0 exp.xmax])
    ylim(lap(l).t)
    title(['leftward lap, lap no ' num2str(lap(l).no)])
    i = i + 2;
end

%% interpolation at each lap
disp('Interpolation ...')
posi.dt = 1/1000; % interpolation 1 kHz
posi.t = [];
posi.x = [];
posi.y = [];
posi.p = [];
posi.dx = [];
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
    interp.t = lap(l).t(1):posi.dt:lap(l).t(2); % it's 1xN should be Nx1
    % interpolation for position
    interp.x = interp1(pos.t, pos.x, interp.t);
    interp.y = interp1(pos.t, pos.y, interp.t);
    interp.p = interp1(pos.t, pos.p, interp.t);
    
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
    posi.y = [posi.y interp.y];
    posi.p = [posi.p; interp.p];
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
    interp.filt.vx = filterlfp(interp.t, interp.vx, 0.01, 10); % cm/sec
    interp.filt.vy = filterlfp(interp.t, interp.vy, 0.01, 10); % cm/sec
    interp.filt.s = vecnorm([interp.filt.vx interp.filt.vy]')';  % speed in cm/sec
    interp.filt.ax = gradient(interp.filt.vx)./gradient(interp.t);   % ax in cm/sec
    interp.filt.ay = gradient(interp.filt.vy)./gradient(interp.t);   % ax in cm/sec
    % angular velocities
    interp.avx = gradient(interp.roll)./gradient(interp.t);   % avx in deg/sec
    interp.avy = gradient(interp.pitch)./gradient(interp.t);   % avy in deg/sec
    interp.avz = gradient(interp.yaw)./gradient(interp.t);   % avz in deg/sec
    interp.filt.avx = filterlfp(interp.t, interp.avx, 0.01, 10); % deg/sec
    interp.filt.avy = filterlfp(interp.t, interp.avy, 0.01, 10); % deg/sec
    interp.filt.avz = filterlfp(interp.t, interp.avz, 0.01, 10); % deg/sec
    
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
posi.s = vecnorm([posi.vx; posi.vy]);
posi.filt.av = vecnorm([posi.filt.avx; posi.filt.avy; posi.filt.avz]);  % angular speed in deg/sec

% x correction for moving frame of referene 
posi.dx = zeros(size(posi.x));
for l = 1:length(lap)
   idx = posi.lap == l;
   posi.dx(idx) = lap(l).corr;
   
end
%% plotting the occupancy
figure(4); clf;
l = 55; %l = randi([1 length(lap)]);
img_filename = [frames(1).folder filesep 'frame-' num2str(lap(l).frame) '.pgm'];
img = imread(img_filename);
img = imlocalbrighten(img);
img = medfilt2(img);
imshow(img);
figure(4)
hold on
h = quiver(posi.x * exp.ppcm,posi.y * exp.ppcm,posi.vx * exp.ppcm,posi.vy * exp.ppcm,2);
plot(repmat(lap(l).gap(1),2,1)* exp.ppcm, ylim* exp.ppcm,'y','LineWidth',2); % lap(l).cross_color
plot(repmat(lap(l).gap(2),2,1)* exp.ppcm, ylim* exp.ppcm,'y','LineWidth',2); % lap(l).cross_color
title('Occupancy with Velocity')
%title(['Rat ' convertStringsToChars(lap(l).status) 'ing ' convertStringsToChars(lap(l).dir) 'ward'])
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'occupancy_with_velocity.jpg'])

%% Directional Speed 
s_thresh = 0; % cm/s
states = sort(unique([lap.status]), 'descend');
no_states = length(states);
for j=0:1
    figure(5+j); clf;
    i = 0;
    a = [];
    for state = states
        for dir = ["left" "right"]
            i = i+1;
            idx = pos.status==state & pos.dir==dir & pos.s >= s_thresh & pos.s <= 450;
            a(i) = subplot(no_states,2,i);
            if j==0
                plot(pos.x(idx),pos.s(idx),'.')
            else
                plot(pos.x(idx)+[lap(pos.lap(idx)).corr]',pos.s(idx),'.');
            end
            ylabel('Speed (cm/s)');
            title([convertStringsToChars(dir) 'ward ' convertStringsToChars(state)])
            if i >= 2*no_states-1
                xlabel('Horizontal position (cm)')
            end
        end
    end
    
    %toc
    set(gcf, 'Position', [0 0 1800 300+200*no_states]);
    linkaxes(a,'xy')
    zoom xon
    xlim([0 exp.xmax])
    if j==0
        sgtitle('Directional speed (stationary frame of ref)');
        saveas(gcf,fullfile(exp_directory, 'Analysis','Directional_speed_stationary.png'))
    else
        sgtitle('Directional speed (moving frame of ref)');
        saveas(gcf,fullfile(exp_directory, 'Analysis','Directional_speed_moving.png'))
    end
end
%% Rate map binning
hist.bin_size = 3; % cm
hist.edges = 0 : hist.bin_size : exp.xmax + max([lap.corr]); % size of image to cm
hist.centers = movmean(hist.edges, 2, 'Endpoints', 'discard');
hist.s_thresh = 5;
% hist.posi = histcounts(posi.x,hist.edges) * posi.dt; % seconds in each bin

%% detect the (more) exact time of jump
figure(7); clf;
for l =1:length(lap)
    idx = (posi.t>lap(l).t_jump-0.2) & (abs(posi.filt.vx)>=50) & (posi.t<lap(l).t_jump + 0.2);
    idx = find(idx,1); % finding last zero
    lap(l).t_jump_exact = posi.t(idx);
end
i = 0;
for dir = ["left" "right"]
    i = i + 1;
    subplot(2,1,i); hold on
    for l=[lap([lap.dir]==dir).no]
        idx = posi.lap == l;
        plot(posi.t(idx)-lap(l).t_jump_exact,abs(posi.filt.vx(idx)),'m');
        xlim([-2 2])
        
        idx = pos.lap == l;
        plot(pos.t(idx)-lap(l).t_jump_exact,50*ones(size(pos.t(idx))),'k');
        title([convertStringsToChars(dir) 'ward laps aligned']);
        ylabel('Horizontal Speed (cm/s)')
        xlim([-2 2])
%         box on
    end
end
    
%% CSC
% theta phase
disp('Theta analysis ...')
start2 = tic;

[csc.t, csc.lfp] = read_bin_csc(csc_filename);
[csc.theta, csc.phase] = filterlfp(csc.t, csc.lfp, 'theta');

for c=1:length(cluster)
    cluster(c).phase = nan*cluster(c).t;
    % x correction for moving frame of referene 
    cluster(c).dx = interp1(posi.t, posi.dx, cluster(c).t);
end
posi.phase = nan * posi.t;

for l=1:length(lap)
    idx1 = csc.t >= lap(l).t(1) & csc.t <= lap(l).t(2);
    for c=1:length(cluster)
        idx2 = [cluster(c).lap]==l;
        if nnz(idx2) > 0 % if there a firing for cluster c in this lap
            cluster(c).phase(idx2) = interp1(csc.t(idx1), csc.phase(idx1), cluster(c).t(idx2));
        end
    end
    idx3 = [posi.lap] == l;
    posi.phase(idx3) = interp1(csc.t(idx1), csc.phase(idx1), posi.t(idx3));
    if mod(l,20)==0
        disp(['lap ' num2str(l) ' from ' num2str(numel(lap)) ' laps.'])
    end
end

% exclude cluster firings outside laps
for c=1:length(cluster)
    idx = isnan(cluster(c).phase); % only ones outside laps are nan
    cluster(c).t(idx) = [];
    cluster(c).x(idx) = [];
    cluster(c).y(idx) = [];
    cluster(c).p(idx,:) = [];
    cluster(c).vx(idx) = [];
    cluster(c).vy(idx) = [];
    cluster(c).s(idx) = [];
    cluster(c).lap(idx) = [];
    cluster(c).dir(idx) = [];
    cluster(c).status(idx) = [];
    cluster(c).phase(idx) = [];
    cluster(c).dx(idx) = [];
end

toc(start2)
%% plot rate map
cluster = plotratemap(exp_directory, posi, lap, cluster, colors, hist, exp);
%% Saving data
mat_filename = fullfile(exp_directory,'processed_data.mat');
if exist(mat_filename, 'file') == 2
    delete(mat_filename);
    warning([mat_filename ' removed']);
end
save(mat_filename,'pos','posi', 'lap', 'cluster','exp','colors','daq','hist','csc', '-v7.3');
fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);
disp(['File ' mat_filename ' has been created!'])
%% Extracting fields
extractfields(exp_directory);
plot_fields(exp_directory);