% This program plots the data for Rat 913
% The data include the occupancy and histogram of the spikes and the
% sequential firing of multiple spikes before jumping
% ## This program is the updated version of plothist ##
% SGL 2021-01-31
clc; close all;
exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
[datafile,exp_directory] = uigetfile(fullfile(exp_directory,'data.mat'), 'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename, 'pos', 'cluster', 'ppcm', 'offset')
clearvars -except datafile exp_directory pos cluster ppcm offset;

[img_file,img_directory] = uigetfile(fullfile(exp_directory,'Videos','*.pgm'), 'Select Image File');
if isequal(img_file, 0)
    error('Image file was not selected!')
end
img_filename = fullfile(img_directory, img_file);

colors = ["#EDB120" "#7E2F8E" "yellow" "#A2142F" "red" "magenta" "green" "#D95319"];
colors = repmat(colors', ceil(length(cluster)/length(colors))); % repeat the colors to match the total number of clusters
colors(35) = "00CCCC";
%colors = lines;

img = imread(img_filename);
xmax = ceil(size(img,2)/ppcm);

tic
%% Plotting the position vs time and marking neural recording times,
% maze times and lap detection times

figure(1)
[lap, idx_analysis] = lapdetector(exp_directory, xmax);

% filtered data:
pos.filt.vx = filtertheta(pos.t, pos.vx, 0.01, 10); % cm/sec
pos.filt.vy = filtertheta(pos.t, pos.vy, 0.01, 10); % cm/sec
pos.filt.s = vecnorm([pos.filt.vx pos.filt.vy]')';  % speed in cm/sec
pos.filt.ax = gradient(pos.filt.vx)./gradient(pos.t);   % ax in cm/sec
pos.filt.ay = gradient(pos.filt.vy)./gradient(pos.t);   % ax in cm/sec

% select only the ones that are in the laps
t = pos.t(idx_analysis);
x = pos.x(idx_analysis);
y = pos.y(idx_analysis);
vx = pos.vx(idx_analysis);
vy = pos.vy(idx_analysis);
s = pos.filt.s(idx_analysis);
frame = pos.frame(idx_analysis);

%% % MODIFY FOR EACH DAY!
% Excluding some clusters (e.g. inter-neurons):
% MODIFY FOR EACH DAY!
cluster_exlude = [1 3];
cluster(cluster_exlude)=[];
disp(['clusters ' num2str(cluster_exlude) ' excluded.'])
for c=1:length(cluster)
    cluster(c).no = c;  % modify this in future
end
% exclude outliers
idx = (t < 2705.005 | t > 2705.025);
t = t(idx);
x = x(idx);
y = y(idx);
vx = vx(idx);
vy = vy(idx);
s = s(idx);
frame = frame(idx);

% Adding gap size to each lap

% mean_gap_len = num2cell([22.0,21.9,21.8,21.7,24.4,24.3,24.6,24.5,24.5,24.6,24.5,24.6,...
%     27.1,27.1,26.8,27.0,26.9,26.9,27.0,26.9,27.1,26.8,29.4,29.2,0.0,29.4,29.4,29.4,...
%     29.7,29.5,29.3,29.4,29.5]);
mean_gap_len = num2cell(repmat(24.5,21,1));
[lap.gap] = mean_gap_len{:};

%% Adding laps and direction to clusters
for c=1:length(cluster)
    cluster(c).lap = zeros(size(cluster(c).t));
    cluster(c).dir = strings(size(cluster(c).t));
    for l = 1:length(lap)
        idx = (cluster(c).t >=lap(l).t(1)) & (cluster(c).t <lap(l).t(2)); % lap
        cluster(c).lap(idx)= lap(l).no;
        cluster(c).dir(idx)=lap(l).dir;
    end
end

%% Plotting time vs horizontal position for all laps
figure(2)
N = max(nnz([lap.dir]=="left"),nnz([lap.dir]=="right"));

% rightward laps
i = 2;
for l=[lap([lap.dir]=="right").no]
    subplot(N,2,i);
    plot(x, t, 'b')
    hold on
    xlim([0 xmax])
    ylim(lap(l).t)
    title(['rightward lap ' num2str(lap(l).no)])
    i = i + 2;
end

% leftward laps
i = 1;
for l=[lap([lap.dir]=="left").no]
    subplot(N,2,i);
    plot(x, t, 'b')
    hold on
    xlim([0 xmax])
    ylim(lap(l).t)
    title(['leftward lap, lap no ' num2str(lap(l).no)])
    i = i + 2;
end

%% plotting the occupancy with arrows showing the velocity
figure(3)
imshow(img+50);
figure(3)
hold on
quiver(x * ppcm,y * ppcm,vx * ppcm,vy * ppcm,20);
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'occupancy_with_velocity.jpg'])
title('Occupancy with Velocity')

%% interpolation at each lap
dt = 1/1000; % interpolation 1 kHz
posi.t = [];
posi.x = [];
posi.y = [];
posi.lap = [];
posi.dir = [];

% Adding laps and direction to pos
pos.lap = zeros(size(pos.t));
pos.dir = strings(size(pos.t));

for l = 1:length(lap)
    
    idx = (pos.t >=lap(l).t(1)) & (pos.t <lap(l).t(2)); % lap
    pos.lap(idx) = lap(l).no;
    pos.dir(idx)= lap(l).dir;
    
    % interpolation
    
    interp.t = lap(l).t(1):dt: lap(l).t(2);
    % interpolation for position
    interp.x = interp1(pos.t, pos.x, interp.t);
    interp.y = interp1(pos.t, pos.y, interp.t);
    
    interp.lap = zeros(size(interp.t));
    interp.dir = strings(size(interp.t));
    interp.lap(:)=lap(l).no;
    interp.dir(:)=lap(l).dir;
    
    posi.t = [posi.t interp.t];
    posi.x = [posi.x interp.x];
    posi.y = [posi.y interp.y];
    posi.dir = [posi.dir interp.dir];
    posi.lap = [posi.lap interp.lap];
end

% figure(50)
% plot(posi.t,posi.x,'.b')

%% Rate map calculations
% Rate map binning
hist.bin_size = 3; % cm
hist.edges = 0:hist.bin_size:xmax; % size of image to cm
hist.posi = histcounts(posi.x,hist.edges) * dt; % seconds in each bin

% Plotting the spikes on the path
for l = 1:length(lap) % each lap
    
    % Figures of the rat in the mid-jump
    figure(3+l)
    set(gcf, 'Position', [100 100 1536 1175]);
    ax1 = subplot(8,1,1:5);
    img = imread([exp_directory filesep 'Videos' filesep 'frame-' num2str(lap(l).frame) '.pgm']);
    imshow(img+50);
    figure(3+l)
    hold on
    idx = (t>=lap(l).t(1) & t<=lap(l).t(2));
    xlap = x(idx);
    ylap = y(idx);
    vxlap = vx(idx);
    vylap = vy(idx);
    quiver(xlap * ppcm,ylap * ppcm,vxlap * ppcm,vylap * ppcm);
    
    % occupancy histogram
    hist.posi = histcounts(posi.x(posi.lap==l), hist.edges) * dt; % seconds in each bin
    ax2 = subplot(8,1,6);
    histogram('BinCounts', hist.posi, 'BinEdges', hist.edges);
    ylabel('Occupancy (sec)')
    xlim([0 xmax])
    
    for c=1:length(cluster) % each cluster
        ax1 = subplot(8,1,1:5);
        h = plot(cluster(c).x([cluster(c).lap]==l) * ppcm, cluster(c).y([cluster(c).lap]==l) * ppcm,...
            'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        
        title(['lap no. ' num2str(lap(l).no) ', cluster no. ' num2str(c) ', ' ...
            num2str(cluster(c).name) '_l' num2str(num2str(l))], 'Interpreter', 'none');
        
        % spike histogram
        ax3 = subplot(8,1,7);
        hist.cluster = histcounts(cluster(c).x([cluster(c).lap]==l), hist.edges); % spikes in each bin
        histogram('BinCounts', hist.cluster, 'BinEdges', hist.edges);
        ylabel('Number of spikes')
        
        % Rate map
        ax4 = subplot(8,1,8);
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
        
        ylabel('Rate Map (Hz)')
        xlabel('Horizontal position (cm)')
        
        linkaxes([ax2 ax3 ax4],'x')
        
        
        saveas(gcf,[exp_directory filesep 'Analysis' filesep 'cluster' num2str(c) '_lap' num2str(l) '.jpg'])
        %set(h, 'Visible','off')
        delete(h);
    end
end

%% CSC
csc_filename= fullfile(exp_directory,'Neuralynx','CSC12.ncs');
addpath('../jumping');
% theta phase
for c=1:length(cluster)
    cluster(c).phase = nan*cluster(c).t;
end
for l=1:length(lap)
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase] = filtertheta(timecsc,lfp);
    
    % looking at all the clusters
    for c=1:length(cluster)
        idx = [cluster(c).lap]==l;
        if nnz(idx) > 0 % if there a firing for cluster c in this lap
            cluster(c).phase(idx) = interp1(timecsc,phase, cluster(c).t(idx));
        end
    end
end

%% Directional rate map for all the laps

for i = 1:3
    figure(200+i)
    for j=1:14
        c = (i-1) * 14 + j;
        if c > length(cluster)
            continue;
        end
        % leftward rate mapj
        hist.posi = histcounts(posi.x(posi.dir=="left"), hist.edges) * dt; % seconds in each bin
        hist.cluster = histcounts(cluster(c).x([cluster(c).dir]=="left"), hist.edges); % spikes in each bin
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        
        a(2*j-1) = subplot(14,2,14*2-(2*j-1));
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
        ylabel(['cluster ' num2str(c)]);
        if j == 14
            title('Rate Map (leftward)')
        elseif j ==1
            xlabel('Horizontal position (cm)')
        end
        
        % rightward rate map
        hist.posi = histcounts(posi.x(posi.dir=="right"), hist.edges) * dt; % seconds in each bin
        hist.cluster = histcounts(cluster(c).x([cluster(c).dir]=="right"), hist.edges); % spikes in each bin
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        
        a(2*j) = subplot(14,2,14*2-(2*j-1)+1);
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
        ylabel(['cluster ' num2str(c)]);
        linkaxes([a(2*j-1) a(2*j)],'y')
        if j == 14
            title('Rate Map (rightward)')
        elseif j ==1
            xlabel('Horizontal position (cm)')
        end
        
    end
    
    set(gcf, 'Position', [100 100 1600 1300]);
    linkaxes(a,'x')
    zoom xon
    xlim([0 xmax])
    saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_ratemap-' num2str(i) '.svg']))
end
toc
%% Saving data
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
save(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax','hist');
disp(['File ' mat_filename ' has been created!'])