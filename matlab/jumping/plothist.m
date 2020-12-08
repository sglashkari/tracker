% This program plots the data for Rat 913
% The data include the occupancy and histogram of the spikes and the
% sequential firing of multiple spikes before jumping
% SGL 2020-11-28
clc; clear; close all;
exp_directory = '~/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02';
exp_directory = '~/Desktop/2020-11-11_Rat913-02';
[datafile,exp_directory] = uigetfile(fullfile(exp_directory,'data.mat'), 'Select Data File');
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename, 'pos', 'spike', 'ppcm')
colors = {'#D95319', '#EDB120', '#7E2F8E', 'yellow', '#A2142F', 'red', 'magenta','green'};

[img_file,img_directory] = uigetfile(fullfile(exp_directory,'Videos','frame-43000.pgm'), 'Select Image File');
img_filename = fullfile(img_directory, img_file);
I = imread(img_filename);

%% Plotting the position vs time and marking neural recording times, 
% maze times and lap detection times

figure(1)
% neural recording times
r1i = 105.221068;
r1f = 703.846944;
r2i = 945.026778;
r2f = 1341.248894;

xmax = round(size(I,2)/ppcm);

hold on
rectangle('Position',[r1i,0,r1f-r1i,xmax],'FaceColor','y')
rectangle('Position',[r2i,0,r2f-r2i,xmax],'FaceColor','y')

% maze times
m1i = 408.650168;
m1f = 677.921799;
m2i = 988.727761;
m2f = 1303.765711;

rectangle('Position',[m1i,0,m1f-m1i,xmax],'FaceColor',[0 .8 .8])
rectangle('Position',[m2i,0,m2f-m2i,xmax],'FaceColor',[0 .8 .8])

% analysis times
t1i = r1i;
t1f = 659;
t2i = r2i;
t2f = 1160; % 1160 will curtail the final lap and 1282 will be the whole version

t_crop_idx = ((pos.t > t1i) & (pos.t < t1f)) | ((pos.t > t2i) & (pos.t < t2f));
t = pos.t(t_crop_idx);
x = pos.x(t_crop_idx);
y = pos.y(t_crop_idx);
vx = pos.vx(t_crop_idx);
vy = pos.vy(t_crop_idx);

plot(pos.t, pos.x, '.b', t, x, '.k')
ylim([0 xmax])

% jump detection
x_thresh = 82;

% times that the rats jump (based on jump direction)
time_jump_leftward = t(diff(x > x_thresh) == -1);
time_jump_rightward = t(diff(x > x_thresh) == 1);
x_jump_leftward = x(diff(x > x_thresh) == -1);
x_jump_rightward = x(diff(x > x_thresh) == 1);

% lap detection
N = length(time_jump_leftward)-1;
t_left = ones(N,1); % lap at left
t_right = ones(N,1); % lap at right
x_left = ones(N,1);
x_right = ones(N,1);

for i=1:N
    % lap
    tl = t(t>=time_jump_leftward(i) & t<=time_jump_leftward(i+1));
    [x_left(i), idx] = min(x(t>time_jump_leftward(i) & t<time_jump_leftward(i+1)));
    t_left(i) = tl(idx);
    [x_right(i), idx] = max(x(t>time_jump_leftward(i) & t<time_jump_leftward(i+1)));
    t_right(i) = tl(idx);
end

% exception for beginning and end
tl = t(t<time_jump_leftward(1));
[x20, idx] = max(x(t<time_jump_leftward(1)));
l20 = tl(idx);
x_right = [x20;x_right];
t_right = [l20;t_right];
tl = t(t>time_jump_leftward(N+1));
[x_left(N+1), idx] = min(x(t>time_jump_leftward(N+1)));
t_left(N+1) = tl(idx);

% picking the frame where the jump has happened
addpath('../tracking');
[time, ~, ~, ~, frame] = readtrackingdata(exp_directory);

for i=1:N+1
    % leftward laps
    lap(2*i-1).dir = 'left';
    lap(2*i-1).no = 2*i-1;
    lap(2*i-1).t = [t_right(i) t_left(i)];
    lap(2*i-1).t_jump = time_jump_leftward(i);  % time of jump in the lap
    lap(2*i-1).frame = frame(time==lap(2*i-1).t_jump); % frame of jump in the lap
    
    % rightward laps
    if i==N+1
        break;
    end
    lap(2*i).dir = 'right';
    lap(2*i).no = 2*i;
    lap(2*i).t = [t_left(i) t_right(i+1)];
    lap(2*i).t_jump = time_jump_rightward(i); % time of jump in the lap
    lap(2*i).frame = frame(time==lap(2*i).t_jump); % frame of jump in the lap
end


plot(t_left, x_left, 'pk', 'MarkerSize',15)
plot(t_right, x_right, 'sk', 'MarkerSize',15)

plot(time_jump_leftward, x_jump_leftward, 'hc', 'MarkerSize',15)
plot(time_jump_rightward, x_jump_rightward, 'hr', 'MarkerSize',15)

set(gcf, 'Position', [100 100 2000 1500]);
saveas(gcf,fullfile(exp_directory, 'Analysis','overall_view.jpg'))
        
%% Plotting time vs horizontal position for all laps
figure(2)
for l = 1:length(lap)
    subplot(N+1,2,lap(l).no);
    plot(x, t, 'b')
    hold on
    xlim([0 xmax])
    ylim(lap(l).t)
end

%% plotting the occupancy with arrows showing the velocity
figure(3)
imshow(2*I+25);
figure(3)
hold on
quiver(x * ppcm,y * ppcm,vx * ppcm,vy * ppcm,2);
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'occupancy_with_velocity.jpg'])
title('Occupancy with Velocity')

%% Fgiures of the rat in the mid-jump
for l = 1:length(lap)
    if l==23
        continue;
    end
    figure(3+l)
    subplot(8,1,1:5);
    I = imread([exp_directory filesep 'Videos' filesep 'frame-' num2str(lap(l).frame) '.pgm']);
    imshow(2*I+25);
    figure(3+l)
    hold on
    idx = (t>=lap(l).t(1) & t<=lap(l).t(2));
    xlap = x(idx);
    ylap = y(idx);
    vxlap = vx(idx);
    vylap = vy(idx);
    quiver(xlap * ppcm,ylap * ppcm,vxlap * ppcm,vylap * ppcm);
    title(['lap no. ' num2str(lap(l).no) ', direction ' lap(l).dir])
end

%% interpolation
dt = 1/1000; % interpolation 1 kHz
posi.t = pos.t(1):dt: pos.t(end);
posi.t = posi.t(((posi.t > t1i) & (posi.t < t1f)) | ((posi.t > t2i) & (posi.t < t2f))); % pos when neural recording is on
% interpolation for position
posi.x = interp1(pos.t, pos.x, posi.t);
posi.y = interp1(pos.t, pos.y, posi.t);

%% Rate map binning
hist.bin_size = 3; % cm
hist.edges = 0:hist.bin_size:xmax; % size of image to cm
hist.posi = histcounts(posi.x,hist.edges) * dt; % seconds in each bin

%% Plotting the spikes on the path
% loooking at all the clusters in the whole experiment (maze 3 is the whole recorded experiment)
for j=1:length(spike([spike.m]==3))
    % only looking at spike.no = j
    cluster = spike([spike.no] == j);
    color = colors{mod(j,length(colors))+1};
    
    for l = 1:length(lap) % each lap
        if l==23
            continue;
        end
        
        figure(3+l)
        ax1 = subplot(8,1,1:5);
        hold on

        this_lap.spike.idx = (cluster.t >=lap(l).t(1) & cluster.t <=lap(l).t(2)); % lap
        this_lap.spike.t = cluster.t(this_lap.spike.idx);
        this_lap.spike.x = cluster.x(this_lap.spike.idx);
        this_lap.spike.y = cluster.y(this_lap.spike.idx);
        
        this_lap.pos.idx = posi.t >=lap(l).t(1) & posi.t <=lap(l).t(2); %lap
        this_lap.pos.t = posi.t(this_lap.pos.idx);
        this_lap.pos.x = posi.x(this_lap.pos.idx);
        this_lap.pos.y = posi.y(this_lap.pos.idx);
        
        h_img = plot(this_lap.spike.x * ppcm, this_lap.spike.y * ppcm,'o','MarkerEdgeColor','black', 'MarkerFaceColor', color);
        title([cluster.name '_l' num2str(l)], 'Interpreter', 'none');
        
        % histogram
        hist.posi = histcounts(this_lap.pos.x,hist.edges) * dt; % seconds in each bin
        ax2 = subplot(8,1,6);
        histogram('BinCounts', hist.posi, 'BinEdges', hist.edges);
        ylabel('Occupancy (sec)')
        
        ax3 = subplot(8,1,7);
        hist.cluster = histcounts(this_lap.spike.x,hist.edges); % spikes in each bin
        histogram('BinCounts', hist.cluster, 'BinEdges', hist.edges);
        ylabel('Number of spikes')
        
        ax4 = subplot(8,1,8);
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',color); %rate map histogram
        
        ylabel('Rate Map (Hz)')
        xlabel('Horizontal position (cm)')
        
        linkaxes([ax2 ax3 ax4],'x')
        xlim([0 xmax])
        set(gcf, 'Position', [100 100 1536 1175]);
        saveas(gcf,[exp_directory filesep 'Analysis' filesep cluster.name '-lap_' num2str(l) '.jpg'])
        set(h_img,'Visible','off')
        
        % storing histogram info for left and right laps
        if strcmp(lap(l).dir,'left')
            if l <= 2
                hist.left.posi = hist.posi;
                hist.left.cluster = hist.cluster;
            else
                hist.left.posi = hist.left.posi + hist.posi;
                hist.left.cluster = hist.left.cluster + hist.cluster;
            end
        else
            if l <= 2
                hist.right.posi = hist.posi;
                hist.right.cluster = hist.cluster;
            else
                hist.right.posi = hist.right.posi + hist.posi;
                hist.right.cluster = hist.right.cluster + hist.cluster;
            end
        end
        
        spike(j).hist = hist;
    end
    
end

%% CSC
csc_filename= fullfile(exp_directory,'Neuralynx','CSC4.ncs');
addpath('../jumping');
for l=1:length(lap)
    if strcmp(lap(l).dir,'left')
        continue;
    end
    figure(100+l)
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    theta = filtertheta(timecsc,lfp);
    
    ax1 = subplot(5,1,1);
    this_lap.pos.idx = pos.t >=lap(l).t(1) & pos.t <=lap(l).t(2); %lap
    this_lap.pos.t = pos.t(this_lap.pos.idx);
    this_lap.pos.x = pos.x(this_lap.pos.idx);
    plot(this_lap.pos.t,this_lap.pos.x);
    ylabel('Horizontal position (cm)')
    title(['lap ' num2str(l)]);
    
    ax2 = subplot(5,1,2);
    this_lap.pos.s = pos.s(this_lap.pos.idx);
    plot(this_lap.pos.t,this_lap.pos.s);
    ylabel('Speed (cm/s)')
    
    ax3 = subplot(5,1,3);
    plot(timecsc,lfp * 1e6);
    ylim([-1000 1000])
    ylabel('LFP (\muV)')  %% micro Volts ???
    
    ax4 = subplot(5,1,4);
    plot(timecsc,theta *1e6,'r');
    ylim([-500 500])
    ylabel('Theta (\muV)')
    
    ax5 = subplot(5,1,5);
    hold off; hold on
    % loooking at all the clusters in the whole experiment (maze 3 is the whole recorded experiment)
    for j=1:length(spike([spike.m]==3))
        
        cluster = spike([spike.no] == j);
        color = colors{mod(j,length(colors))+1};
        
        this_lap.spike.idx = cluster.t >=lap(l).t(1) & cluster.t <=lap(l).t(2); % lap
        this_lap.spike.t = cluster.t(this_lap.spike.idx);
        this_lap.spike.no = j * ones(size(this_lap.spike.t));
        
        plot(this_lap.spike.t, this_lap.spike.no,'o','MarkerEdgeColor','black', 'MarkerFaceColor', color);
    end
    xlabel('Time (sec)')
    ylim([0 length(spike([spike.m]==3))+1])
    ylabel('Cluster Number')
    
    linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
    xlim(lap(l).t);
    
    set(gcf, 'Position', [100 100 1536 1175]);
    saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-lap_' num2str(l) '.jpg']))
    
end

%% Directional rate map for all the laps
figure(200)

for j=1:length(spike([spike.m]==3))
    color = colors{mod(j,length(colors))+1};
    
    % fix all spike([spike.no==j]) to spike(j)
    
    hist.left.ratemap = spike(j).hist.left.cluster ./ (spike(j).hist.left.posi + eps); % adding eps to avoid division by zero
    hist.right.ratemap = spike(j).hist.right.cluster ./ (spike(j).hist.right.posi + eps); % adding eps to avoid division by zero
    
    % leftward rate map
    a(2*j-1) = subplot(length(spike([spike.m]==3)),2,2*length(spike([spike.m]==3))-(2*j-1));
    histogram('BinCounts', hist.left.ratemap, 'BinEdges', hist.edges, 'FaceColor',color); %rate map histogram
    ylabel(['cluster ' num2str(j)]);
    
    if j == length(spike([spike.m]==3))
        title('Rate Map (left)')
    end
    if j ==1 
        xlabel('Horizontal position (cm)')
    end
    
    % rightward rate map
    a(2*j) = subplot(length(spike([spike.m]==3)),2,2*length(spike([spike.m]==3))-(2*j-1)+1);
    histogram('BinCounts', hist.right.ratemap, 'BinEdges', hist.edges, 'FaceColor',color); %rate map histogram
    ylabel(['cluster ' num2str(j)]);
    
    linkaxes([a(2*j-1) a(2*j)],'y')
    
    if j == length(spike([spike.m]==3))
        title('Rate Map (Right)')
    end
    if j ==1 
        xlabel('Horizontal position (cm)')
    end
    
    ylabel(['cluster ' num2str(j)]);
    zoom xon
end

set(gcf, 'Position', [100 100 1600 1300]);
linkaxes(a,'x')
xlim([0 xmax]) % fix [hist.edges(1) hist.edges(end)] to [0 xmax]
%xlabel('Horizontal position (cm)')
saveas(gcf,fullfile(exp_directory, 'Analysis','Directional_ratemap.jpg'))

%% Cross correlation

%{
figure(201)
sampling_freq = 1e3; %30KHz
l = 12;
t_lap = lap(l).t(1):1/sampling_freq:lap(l).t(2);% lap time
tjl = time_jump_rightward(l/2);

i = 2;
timerange = tjl -i-6:sampling_freq:tjl -i-5;

% for j=1:length(spike([spike.m]==3))
%     sp(j).idx = ismember(rec_t,spike(1).t);
% end
sp1.t = spike(1).t(spike(1).t >= (tjl +i-6)  & spike(1).t < (tjl+i-5) );
sp3.t = spike(3).t(spike(3).t >= (tjl +i-6)  & spike(3).t < (tjl+i-5) );

ax1 = subplot(4,1,1);
[xxx,zzz] = histcounts(spike(1).t,t_lap);
histogram('BinCounts', xxx, 'BinEdges', zzz)

ax2 = subplot(4,1,2);
[yyy,zzz] = histcounts(spike(3).t,t_lap);
histogram('BinCounts', yyy, 'BinEdges', zzz)

ax3 = subplot(4,1,3);
zzz(end) = [];
%[c,lags] = xcorr(spike(1).t,spike(3).t);
[c,lags] = xcorr(xxx, yyy);
stem(lags,c)

ax4 = subplot(4,1,4);
[c,lags] = xcorr(xxx);
stem(lags,c)

%[c,lags] = xcorr(sp(1).idx);
%plot(lags,c,'.')

%linkaxes([ax1 ax2], 'x')
%}
%% Interspike interval
figure(202)
for l=2:2:30
    
    sp1.t = spike(1).t(spike(1).t >= lap(l).t(1)  & spike(1).t < lap(l).t(2) );
    sp3.t = spike(3).t(spike(3).t >= lap(l).t(1)  & spike(3).t < lap(l).t(2) );
    t_min = min(sp1.t(5),sp1.t(5));
    t_max = max(sp1.t(end-4),sp1.t(end-4));
    t_mid = mean([t_min t_max]);
    tjl = time_jump_rightward(l/2);
    for i=1:1
        %subplot(11,1,2*i-1);
        figure(202)
        subplot(8,2,round(l/2))
        sp1.t = spike(1).t(spike(1).t >= t_min & spike(1).t < t_mid );
        sp3.t = spike(3).t(spike(3).t >= t_min & spike(3).t < t_mid );
        %      sp1.t = spike(1).t(spike(1).t >= (tjl +i-1)  & spike(1).t < (tjl+i-0) );
        %      sp3.t = spike(3).t(spike(3).t >= (tjl +i-1)  & spike(3).t < (tjl+i-0) );
        [xx,yy] = meshgrid(sp1.t, sp3.t);
        zz = xx - yy;
        range = 0.1;
        zz((zz>range) | (zz<-range)) = nan;
        histogram(zz,(-range:0.004:range))
        
        figure(203)
        subplot(8,2,round(l/2))
        sp1.t = spike(1).t(spike(1).t >= t_mid & spike(1).t < t_max );
        sp3.t = spike(3).t(spike(3).t >= t_mid & spike(3).t < t_max );
        %      sp1.t = spike(1).t(spike(1).t >= (tjl +i-1)  & spike(1).t < (tjl+i-0) );
        %      sp3.t = spike(3).t(spike(3).t >= (tjl +i-1)  & spike(3).t < (tjl+i-0) );
        [xx,yy] = meshgrid(sp1.t, sp3.t);
        zz = xx - yy;
        range = 0.1;
        zz((zz>range) | (zz<-range)) = nan;
        histogram(zz,(-range:0.004:range))        
    end
end