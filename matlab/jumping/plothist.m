% This program plots the data for Rat 913
% The data include the occupancy and histogram of the spikes and the
% sequential firing of multiple spikes before jumping
% SGL 2020-11-28
clc; close all;
exp_directory = '~/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02';
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

% lap detection times
t1i = 110;
t1f = 659;
t2i = 945;
t2f = 1282;

t_crop_idx = ((pos.t > t1i) & (pos.t < t1f)) | ((pos.t > t2i) & (pos.t < t2f));
t = pos.t(t_crop_idx);
x = pos.x(t_crop_idx);
y = pos.y(t_crop_idx);
vx = pos.vx(t_crop_idx);
vy = pos.vy(t_crop_idx);

plot(t, x, '.')
ylim([0 xmax])

% jump detection
x_thresh = 82;

% times that the rats jump (left to right and right to left)
time_right_jump = t(diff(x > x_thresh) == -1);
time_left_jump = t(diff(x > x_thresh) == 1);
x_right_jump = x(diff(x > x_thresh) == -1);
x_left_jump = x(diff(x > x_thresh) == 1);

% lap detection
N = length(time_right_jump)-1;
left = ones(N,1); % lap at left
right = ones(N,1); % lap at right
x1 = ones(N,1);
x2 = ones(N,1);

for i=1:N
    % lap
    tl = t(t>=time_right_jump(i) & t<=time_right_jump(i+1));
    [x1(i), idx] = min(x(t>time_right_jump(i) & t<time_right_jump(i+1)));
    left(i) = tl(idx);
    [x2(i), idx] = max(x(t>time_right_jump(i) & t<time_right_jump(i+1)));
    right(i) = tl(idx);
end

% exception for beginning and end
tl = t(t<time_right_jump(1));
[x20, idx] = max(x(t<time_right_jump(1)));
l20 = tl(idx);
x2 = [x20;x2];
right = [l20;right];
tl = t(t>time_right_jump(N+1));
[x1(N+1), idx] = min(x(t>time_right_jump(N+1)));
left(N+1) = tl(idx);

plot(left, x1, 'pk', 'MarkerSize',15)
plot(right, x2, 'sk', 'MarkerSize',15)

plot(time_right_jump, x_right_jump, 'hg', 'MarkerSize',15)
plot(time_left_jump, x_left_jump, 'hr', 'MarkerSize',15)

set(gcf, 'Position', [100 100 2000 1200]);
saveas(gcf,[exp_directory filesep 'overall_view.jpg'])
        
%% Plotting time vs horizontal position for all laps
figure(2)
for i = 1:N+1
    subplot(N+1,2,2*i-1);
    plot(x, t)
    xlim([0 xmax])
    ylim([right(i) left(i)])
    if i==N+1
        break;
    end
    subplot(N+1,2,2*i);
    plot(x, t)
    xlim([0 xmax])
    ylim([left(i) right(i+1)])
end

%% plotting the occupancy with arrows showing the velocity
figure(3)

imshow(2*I+50);
figure(3)
hold on
quiver(x * ppcm,y * ppcm,vx * ppcm,vy * ppcm,2);

set(gcf, 'Position', [100 100 1536 1175]);
saveas(gcf,[exp_directory filesep 'occupancy_with_velocity.jpg'])

%% picking the frame where the jump has happened (left to right jump: jr)
addpath('../tracking');
[time, ~, ~, ~, frame] = readtrackingdata(exp_directory);
jump_frames = ones(N+1,2);
for i=1:N+1
    jump_frames(i,1) = frame(time==time_right_jump(i));
    if i==N+1
        break;
    end
    jump_frames(i,2) = frame(time==time_left_jump(i));
end

%% Fgiures of the rat in the mid-jump
lap_no = zeros(2*N+2,1);
figure(4)
for i = 1:N+1
    figure(2*(i+1));
    subplot(8,1,1:5);
    I = imread([exp_directory filesep 'Videos' filesep 'frame-' num2str(jump_frames(i,1)) '.pgm']);
    imshow(2*I+25);
    figure(2*(i+1));
    hold on
    lap.idx = (t>=right(i) & t<=left(i));
    xlap = x(lap.idx);
    ylap = y(lap.idx);
    vxlap = vx(lap.idx);
    vylap = vy(lap.idx);
    quiver(xlap * ppcm,ylap * ppcm,vxlap * ppcm,vylap * ppcm);
    %saveas(gcf,[exp_directory filesep 'Analysis' filesep num2str(2*i-1) '.ppm'])
    lap_no(2*i-1) = 1;
    if i==N+1
        continue;
    end
    figure(2*(i+1)+1);
    subplot(8,1,1:5);
    if i == 7
        I = imread([exp_directory filesep 'Videos' filesep 'frame-' num2str(jump_frames(i-1,2)) '.pgm']);
    else
        I = imread([exp_directory filesep 'Videos' filesep 'frame-' num2str(jump_frames(i,2)) '.pgm']);
    end
    imshow(2*I+25);
    figure(2*(i+1)+1);
    hold on
    lap.idx = (t>=left(i) & t<=right(i+1));
    xlap = x(lap.idx);
    ylap = y(lap.idx);
    vxlap = vx(lap.idx);
    vylap = vy(lap.idx);
    quiver(xlap * ppcm,ylap * ppcm,vxlap * ppcm,vylap * ppcm);
    %saveas(gcf,[exp_directory filesep 'Analysis' filesep num2str(2*i) '.ppm'])
    lap_no(2*i) = 1;
end

%% interpolation
dt = 1/1000; % interpolation 1 kHz
posi.t = pos.t(1):dt: pos.t(end);
posi.t = posi.t(((posi.t > r1i) & (posi.t < r1f)) | ((posi.t > r2i) & (posi.t < r2f))); % pos when neural recording is on
% interpolation for position
posi.x = interp1(pos.t, pos.x, posi.t);
posi.y = interp1(pos.t, pos.y, posi.t);

%% Rate map histogram for each cluster for the whole experiment
figure(100)
for j=length(spike([spike.m]==3)):-1:1
    a(j) = subplot(length(spike([spike.m]==3)),1,length(spike([spike.m]==3))-j+1);
    color = colors{mod(j,length(colors))+1};
    hist.bin_size = 2; % 2 cm
    hist.edges = 0:hist.bin_size:size(I,2)/ppcm; % size of image to cm
    hist.posi = histcounts(posi.x,hist.edges) * dt; % seconds in each bin (2 cm)
    cluster = spike([spike.no] == j);
    hist.cluster = histcounts(cluster.x,hist.edges); % spikes in each bin (2 cm)
    hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
    histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',color); %rate map histogram
    if j == length(spike([spike.m]==3))
        title('Rate Map (Hz)')
    end
    ylabel(['cluster ' num2str(cluster.no)]);
    zoom xon
    
end
set(gcf, 'Position', [100 100 1600 1300]);
linkaxes(a,'x')
xlim([hist.edges(1) hist.edges(end)])
xlabel('Horizontal position (cm)')
saveas(gcf,[exp_directory filesep 'Ratemap.jpg'])

%% Plotting the spikes on the path
% loooking at all the clusters in the whole experiment (maze 3 is the whole recorded experiment)
for j=1:length(spike([spike.m]==3))
    % only looking at spike.no = j
    cluster = spike([spike.no] == j);
    color = colors{mod(j,length(colors))+1};
    
    for i = 1:N+1 % lap (should be written in a better format)
        figure(2*(i+1));
        ax1 = subplot(8,1,1:5);
        hold on
        lap.spike.idx = (cluster.t >=right(i) & cluster.t <=left(i)); % lap
        lap.spike.t = cluster.t(lap.spike.idx);
        lap.spike.x = cluster.x(lap.spike.idx);
        lap.spike.y = cluster.y(lap.spike.idx);
        
        lap.pos.idx = posi.t >=right(i) & posi.t <=left(i); %lap
        lap.pos.t = posi.t(lap.pos.idx);
        lap.pos.x = posi.x(lap.pos.idx);
        lap.pos.y = posi.y(lap.pos.idx);
        
        h_img = plot(lap.spike.x * ppcm, lap.spike.y * ppcm,'o','MarkerEdgeColor','black', 'MarkerFaceColor', color);
        title([cluster.name '_l' num2str(2*i-1)], 'Interpreter', 'none');
        
        % histogram
        hist.posi = histcounts(lap.pos.x,hist.edges) * dt; % seconds in each bin
        ax2 = subplot(8,1,6);
        histogram('BinCounts', hist.posi, 'BinEdges', hist.edges);
        ylabel('Occupancy (sec)')
        
        ax3 = subplot(8,1,7);
        hist.cluster = histcounts(lap.spike.x,hist.edges); % spikes in each bin
        histogram('BinCounts', hist.cluster, 'BinEdges', hist.edges);
        ylabel('Number of spikes')
        
        ax4 = subplot(8,1,8);
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',color); %rate map histogram
        
        ylabel('Rate Map (Hz)')
        xlabel('Horizontal position (cm)')
        
        linkaxes([ax2 ax3 ax4],'x')
        xlim([hist.edges(1) hist.edges(end)])
        set(gcf, 'Position', [100 100 1536 1175]);
        saveas(gcf,[exp_directory filesep 'Analysis' filesep cluster.name '-lap_' num2str(2*i-1) '.jpg'])
        set(h_img,'Visible','off')
        
        if i==N+1
            continue;
        end
        figure(2*(i+1)+1);
        ax1 = subplot(8,1,1:5);
        hold on
        lap.spike.idx = (cluster.t >=left(i) & cluster.t <=right(i+1));
        lap.spike.t = cluster.t(lap.spike.idx);
        lap.spike.x = cluster.x(lap.spike.idx);
        lap.spike.y = cluster.y(lap.spike.idx);
        
        lap.pos.idx = ((posi.t > t1i) & (posi.t < t1f)) | ((posi.t > t2i) & (posi.t < t2f)); % pos when neural recording is on
        lap.pos.idx = posi.t >=left(i) & posi.t <=right(i+1); %lap
        lap.pos.t = posi.t(lap.pos.idx);
        lap.pos.x = posi.x(lap.pos.idx);
        lap.pos.y = posi.y(lap.pos.idx);
        
        h_img = plot(lap.spike.x * ppcm, lap.spike.y * ppcm,'o','MarkerEdgeColor','black', 'MarkerFaceColor', color);
        title([cluster.name '_l' num2str(2*i)], 'Interpreter', 'none');
        
        % histogram
        hist.posi = histcounts(lap.pos.x,hist.edges) * dt; % seconds in each bin
        ax2 = subplot(8,1,6);
        histogram('BinCounts', hist.posi, 'BinEdges', hist.edges);
        ylabel('Occupancy (sec)')
        xlim([hist.edges(1) hist.edges(end)])
        
        ax3 = subplot(8,1,7);
        hist.cluster = histcounts(lap.spike.x,hist.edges); % spikes in each bin
        histogram('BinCounts', hist.cluster, 'BinEdges', hist.edges);
        ylabel('Number of spikes')
        
        ax4 = subplot(8,1,8);
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',color); %rate map histogram
        
        ylabel('Rate Map (Hz)')
        xlabel('Horizontal position (cm)')
        
        set(gcf, 'Position', [100 100 1536 1175]);
        linkaxes([ax2 ax3 ax4],'x')
        xlim([hist.edges(1) hist.edges(end)])
        
        saveas(gcf,[exp_directory filesep 'Analysis' filesep cluster.name '-lap_' num2str(2*i) '.jpg'])
        set(h_img,'Visible','off')
    end
end

%% CSC
for i = 1:N
    figure(100+i)
    csc_filename= fullfile(exp_directory,'Neuralynx','CSC4.ncs');
    timerange = [left(i) right(i+1)];
    addpath('../jumping');
    [timecsc,lfp] = readcsc(csc_filename, timerange * 1e6); % microseconds
    theta = filtertheta(timecsc,lfp, 5, 12);
    
    ax1 = subplot(5,1,1);
    lap.pos.idx = pos.t >=left(i) & pos.t <=right(i+1); %lap
    lap.pos.t = pos.t(lap.pos.idx);
    lap.pos.x = pos.x(lap.pos.idx);
    plot(lap.pos.t,lap.pos.x);
    ylabel('Horizontal position (cm)')
    title(['lap ' num2str(2*i)]);
    
    ax2 = subplot(5,1,2);
    lap.pos.s = pos.s(lap.pos.idx);
    plot(lap.pos.t,lap.pos.s);
    ylabel('Speed (cm/s)')
    
    ax3 = subplot(5,1,3);
    plot(timecsc,lfp * 1e3);
    ylim([-1 1])
    ylabel('LFP (\muV)')
    
    ax4 = subplot(5,1,4);
    plot(timecsc,theta *1e3,'r');
    ylim([-0.5 0.5])
    ylabel('Theta (\muV)')
    
    ax5 = subplot(5,1,5);
    hold off; hold on
    % loooking at all the clusters in the whole experiment (maze 3 is the whole recorded experiment)
    for j=1:length(spike([spike.m]==3))
        
        cluster = spike([spike.no] == j);
        color = colors{mod(j,length(colors))+1};
        
        lap.spike.idx = cluster.t >=left(i) & cluster.t <=right(i+1); % lap
        lap.spike.t = cluster.t(lap.spike.idx);
        lap.spike.no = j * ones(size(lap.spike.t));
        
        plot(lap.spike.t, lap.spike.no,'o','MarkerEdgeColor','black', 'MarkerFaceColor', color);
    end
    xlabel('Time (sec)')
    ylim([0 length(spike([spike.m]==3))+1])
    ylabel('Cluster Number')
    
    linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
    xlim(timerange);
    
    set(gcf, 'Position', [100 100 1536 1175]);
    saveas(gcf,[exp_directory filesep 'Analysis' filesep 'CSC-lap_' num2str(2*i) '.jpg'])

end