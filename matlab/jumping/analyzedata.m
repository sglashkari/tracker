% This program plots the data for Rat 980 (originally 913)
% The data include the occupancy and histogram of the spikes and the
% sequential firing of multiple spikes before jumping
%
%
% SGL 2022-01-26 (originally 2021-01-31)
%
clc; clear; close all;
[datafile,exp_directory] = uigetfile(fullfile('D:\Analysis', 'data.mat'), 'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename, 'pos', 'exp', 'cluster', 'ppcm', 'offset')

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

tic

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

%% Enhance pos and lap
pos.lap = zeros(size(pos.t));
pos.dir = strings(size(pos.t));
for l = 1:length(lap)
    % Adding laps and direction to pos
    idx = (pos.t >=lap(l).t(1)) & (pos.t <lap(l).t(2)); % lap
    pos.lap(idx) = lap(l).no;
    pos.dir(idx)= lap(l).dir;
end

for l = 1:length(lap)
    % Adding gap size to each lap
    lap(l).gap = 888.651/ppcm + [0 pos.len(pos.t == lap(l).t_jump)];  % two edges of the gap (change for each day!!!) 11 inch = 29.845 cm
    % time of crossing the gap edges:
    idx = pos.lap==l;
    gapidx = pos.x(idx)>lap(l).gap(1) & pos.x(idx)<=lap(l).gap(2);
    gapidx = [0; diff(gapidx)];
    post = pos.t(idx);
    lap(l).t_cross = post(gapidx~=0)';
end 

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
    title(['rightward lap no ' num2str(lap(l).no)])
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
figure(3); clf;
img = imread(img_filename);
img = imlocalbrighten(img);
img = medfilt2(img);
imshow(img);
figure(3)
hold on
tic
h = quiver(x * ppcm,y * ppcm,vx * ppcm,vy * ppcm,2);
toc
plot(repmat(lap(1).gap(1),2,1)* ppcm, ylim* ppcm,'y','LineWidth',2);
plot(repmat(lap(1).gap(2),2,1)* ppcm, ylim* ppcm,'y','LineWidth',2);
title('Occupancy with Velocity')
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'occupancy_with_velocity.jpg'])

%% interpolation at each lap
dt = 1/1000; % interpolation 1 kHz
posi.t = [];
posi.x = [];
posi.y = [];
posi.vx = [];
posi.vy = [];
posi.lap = [];
posi.dir = [];

for l = 1:length(lap)
    
    % interpolation
    
    interp.t = lap(l).t(1):dt: lap(l).t(2);
    % interpolation for position
    interp.x = interp1(pos.t, pos.x, interp.t);
    interp.y = interp1(pos.t, pos.y, interp.t);
    % interpolation for velocity
    interp.vx = interp1(pos.t, pos.vx, interp.t);
    interp.vy = interp1(pos.t, pos.vy, interp.t);    
    
    interp.lap = zeros(size(interp.t));
    interp.dir = strings(size(interp.t));
    interp.lap(:)=lap(l).no;
    interp.dir(:)=lap(l).dir;
    
    posi.t = [posi.t interp.t];
    posi.x = [posi.x interp.x];
    posi.y = [posi.y interp.y];
    posi.vx = [posi.vx interp.vx];
    posi.vy = [posi.vy interp.vy];
    posi.dir = [posi.dir interp.dir];
    posi.lap = [posi.lap interp.lap];
end

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

%% Rate map calculations
% Rate map binning
hist.bin_size = 3; % cm
hist.edges = 0:hist.bin_size:xmax; % size of image to cm
hist.posi = histcounts(posi.x,hist.edges) * dt; % seconds in each bin

% Plotting the spikes on the path
for l = 1:length(lap) % each lap
    
    % Figures of the rat in the mid-jump
    figure(10+l)
    set(gcf, 'Position', [100 100 1770 1000]);
    ax1 = subplot(5,1,1:2);
    img = imread([exp_directory filesep 'frames' filesep 'frame-' num2str(lap(l).frame) '.pgm']);
    img = imlocalbrighten(img);
    img = medfilt2(img);
    imshow(img);
    figure(10+l)
    hold on
    idx = (t>=lap(l).t(1) & t<=lap(l).t(2));
    xlap = x(idx);
    ylap = y(idx);
    vxlap = vx(idx);
    vylap = vy(idx);
    quiver(xlap * ppcm,ylap * ppcm,vxlap * ppcm,vylap * ppcm);
    
    % occupancy histogram
    hist.posi = histcounts(posi.x(posi.lap==l), hist.edges) * dt; % seconds in each bin
    ax2 = subplot(5,1,3);
    histogram('BinCounts', hist.posi, 'BinEdges', hist.edges);
    ylabel('Occupancy (sec)')
    xlim([0 xmax])
    
    for c=1:length(cluster) % each cluster
        ax1 = subplot(5,1,1:2);
        h = plot(cluster(c).x([cluster(c).lap]==l) * ppcm, cluster(c).y([cluster(c).lap]==l) * ppcm,...
            'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        
        title([num2str(cluster(c).region) ': Cluster no. ' num2str(c) ', Lap no. ' num2str(lap(l).no) ', ' ... 
            cluster(c).name ], 'Interpreter', 'none');
        
        % spike histogram
        ax3 = subplot(5,1,4);
        hist.cluster = histcounts(cluster(c).x([cluster(c).lap]==l), hist.edges); % spikes in each bin
        histogram('BinCounts', hist.cluster, 'BinEdges', hist.edges);
        ylabel('Number of spikes')
        
        % Rate map
        ax4 = subplot(5,1,5);
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
        
        ylabel('Rate Map (Hz)')
        xlabel('Horizontal position (cm)')
        
        linkaxes([ax2 ax3 ax4],'x')
        xlim([0 xmax])
        
        saveas(gcf,[exp_directory filesep 'Analysis' filesep 'cluster' num2str(c) '_lap' num2str(l) '.jpg'])
        %set(h, 'Visible','off')
        delete(h);
    end
end

%% CSC
% theta phase
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

%% Directional rate map for all the laps
% groups of 14
for i = 1:ceil(length(cluster)/14)
    figure(200+i); clf;
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
        ylabel(['#' num2str(c) ':' cluster(c).region]);
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
        ylabel(['#' num2str(c) ':' cluster(c).region]);
        linkaxes([a(2*j-1) a(2*j)],'y')
        ylim(max([0 5],ylim))
        if j == 14
            title('Rate Map (rightward)')
        elseif j ==1
            xlabel('Horizontal position (cm)')
        end
        
        % edge of gap
        hold on
        plot(repmat(lap(l).gap(1),2,1), ylim,'b','LineWidth',2);
        plot(repmat(lap(l).gap(2),2,1), ylim,'b','LineWidth',2);
        a(2*j-1) = subplot(14,2,14*2-(2*j-1));
        hold on
        plot(repmat(lap(l).gap(1),2,1), ylim,'b','LineWidth',2);
        plot(repmat(lap(l).gap(2),2,1), ylim,'b','LineWidth',2);
    end
    
    set(gcf, 'Position', [100 100 1600 1100]);
    linkaxes(a,'x')
    zoom xon
    xlim([0 xmax])
    saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_ratemap-' num2str(i) '.jpg']))
end
%% Saving data
fprintf(['It totally took ' datestr(seconds(toc),'HH:MM:SS') ,'.\n']);
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
save(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax','hist');
disp(['File ' mat_filename ' has been created!'])