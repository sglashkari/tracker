% This program plots the data for Rat 913
% The data include the occupancy and histogram of the spikes and the
% sequential firing of multiple spikes before jumping
% SGL 2020-11-28
clc; close all;
exp_directory = '/home/shahin/Desktop/20-12-09';
[datafile,exp_directory] = uigetfile(fullfile(exp_directory,'data.mat'), 'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename, 'pos', 'spike', 'ppcm', 'offset')
clearvars -except datafile exp_directory pos spike ppcm offset;

[img_file,img_directory] = uigetfile(fullfile(exp_directory,'Videos','*.pgm'), 'Select Image File');
if isequal(img_file, 0)
    error('Image file was not selected!')
end
img_filename = fullfile(img_directory, img_file);

colors = ["#EDB120" "#7E2F8E" "yellow" "#A2142F" "red" "magenta" "green" "#D95319"];
colors = repmat(colors', ceil(length(spike)/length(colors))); % repeat the colors to match the total number of spikes

I = imread(img_filename);
xmax = ceil(size(I,2)/ppcm);

%% Plotting the position vs time and marking neural recording times, 
% maze times and lap detection times

figure(1)
% neural recording times
r1i = 105.221068;
r1f = 703.846944;
r2i = 945.026778;
r2f = 1341.248894;

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
t1i = 120; %r1i;
t1f = r1f;
t2i = r2i;
t2f = 1200; % 1160 will curtail the final lap and 1282 will be the whole version

t_crop_idx = ((pos.t > t1i) & (pos.t < t1f)) | ((pos.t > t2i) & (pos.t < t2f));
t = pos.t(t_crop_idx);
x = pos.x(t_crop_idx);
y = pos.y(t_crop_idx);
vx = pos.vx(t_crop_idx);
vy = pos.vy(t_crop_idx);
s = pos.s(t_crop_idx);
frame = pos.frame(t_crop_idx);

plot(pos.t, pos.x, '.b', t, x, '.k')
ylim([0 xmax])

% jump detection
x_thresh = 82;

% times that the rats jump (based on jump direction)
jump_criteria_lefttward = (diff(x > x_thresh) == -1) & (vx(2:end) < -50);
jump_criteria_rightward = (diff(x > x_thresh) == 1) & (vx(2:end) > 50);
time_jump_leftward = t(jump_criteria_lefttward);
time_jump_rightward = t(jump_criteria_rightward);
x_jump_leftward = x(jump_criteria_lefttward);
x_jump_rightward = x(jump_criteria_rightward);

% lap detection
N = length(time_jump_leftward)-1;
t_left = ones(N,1); % lap at left
t_right = ones(N,1); % lap at right
x_left = ones(N,1);
x_right = ones(N,1);

for i=1:N+1
    % exception handling for break time
    if time_jump_leftward(i) < r1f && time_jump_leftward(i+1) > r2i
        time_jump_leftward = [time_jump_leftward(1:i); r2i; time_jump_leftward(i+1:end)];
        x_jump_leftward = [x_jump_leftward(1:i); x_thresh; x_jump_leftward(i+1:end)];
    end

    %lap
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
tl = t(t>time_jump_leftward(N+2));
[x_left(N+2), idx] = min(x(t>time_jump_leftward(N+2)));
t_left(N+2) = tl(idx);

% picking the frame where the jump has happened
for i=1:N+1
    % leftward laps
    lap(2*i-1).dir = "left";
    lap(2*i-1).no = 2*i-1;
    lap(2*i-1).t = [t_right(i) t_left(i)];
    lap(2*i-1).t_jump = time_jump_leftward(i);  % time of jump in the lap
    lap(2*i-1).frame = pos.frame(pos.t==lap(2*i-1).t_jump); % frame of jump in the lap
    
    % rightward laps
    if i==N+1
        break;
    end
    lap(2*i).dir = "right";
    lap(2*i).no = 2*i;
    lap(2*i).t = [t_left(i) t_right(i+1)];
    lap(2*i).t_jump = time_jump_rightward(i); % time of jump in the lap
    lap(2*i).frame = pos.frame(pos.t==lap(2*i).t_jump); % frame of jump in the lap
end


plot(t_left, x_left, 'pk', 'MarkerSize',15)
plot(t_right, x_right, 'sk', 'MarkerSize',15)

plot(time_jump_leftward, x_jump_leftward, 'hc', 'MarkerSize',15)
plot(time_jump_rightward, x_jump_rightward, 'hr', 'MarkerSize',15)

set(gcf, 'Position', [100 100 2000 1500]);
saveas(gcf,fullfile(exp_directory, 'Analysis','overall_view.jpg'))

%% Adding gap size to each lap

mean_gap_len = num2cell([22.0,21.9,21.8,21.7,24.4,24.3,24.6,24.5,24.5,24.6,24.5,24.6,...
    27.1,27.1,26.8,27.0,26.9,26.9,27.0,26.9,27.1,26.8,29.4,29.2,0.0,29.4,29.4,29.4,...
    29.7,29.5,29.3,29.4,29.5]);
[lap.gap] = mean_gap_len{:};

%% Adding laps and direction to spike for m3 maze
for j=1:length(spike([spike.m]==3))
    spike(j).lap = zeros(size(spike(j).t));
    spike(j).dir = strings(size(spike(j).t));
    for l = 1:length(lap)
        idx = (spike(j).t >=lap(l).t(1)) & (spike(j).t <lap(l).t(2)); % lap
        spike(j).lap(idx)= lap(l).no;
        spike(j).dir(idx)=lap(l).dir;
    end
end

%% Plotting time vs horizontal position for all laps
figure(2)
for l = 1:length(lap) % for future l=[lap.no]
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
    if l==25
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

%% Adding laps and direction to pos and posi for m3 maze
pos.lap = zeros(size(pos.t));
pos.dir = strings(size(pos.t));
posi.lap = zeros(size(posi.t));
posi.dir = strings(size(posi.t));
for l = 1:length(lap)
    idx = (pos.t >=lap(l).t(1)) & (pos.t <lap(l).t(2)); % lap
    pos.lap(idx) = lap(l).no;
    pos.dir(idx)= lap(l).dir;
    
    idx = (posi.t >=lap(l).t(1) & posi.t <lap(l).t(2)); % lap
    posi.lap(idx)=lap(l).no;
    posi.dir(idx)=lap(l).dir;
    
end
%% Rate map binning
hist.bin_size = 3; % cm
hist.edges = 0:hist.bin_size:xmax; % size of image to cm
hist.posi = histcounts(posi.x,hist.edges) * dt; % seconds in each bin
%% Plotting the spikes on the path
% looking at all the clusters in m3 (m3 is the whole recorded experiment)
for j=1:length(spike([spike.m]==3))
    for l = 1:length(lap) % each lap
        if l==25
            continue;
        end
        
        figure(3+l)
        ax1 = subplot(8,1,1:5);
        hold on

        h = plot(spike(j).x([spike(j).lap]==l) * ppcm, spike(j).y([spike(j).lap]==l) * ppcm,...
            'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
        title([spike(j).name '_l' num2str(l)], 'Interpreter', 'none');
        
        % histogram
        hist.posi = histcounts(posi.x(posi.lap==l), hist.edges) * dt; % seconds in each bin
        ax2 = subplot(8,1,6);
        histogram('BinCounts', hist.posi, 'BinEdges', hist.edges);
        ylabel('Occupancy (sec)')
        
        ax3 = subplot(8,1,7);
        hist.cluster = histcounts(spike(j).x([spike(j).lap]==l), hist.edges); % spikes in each bin
        histogram('BinCounts', hist.cluster, 'BinEdges', hist.edges);
        ylabel('Number of spikes')
        
        ax4 = subplot(8,1,8);
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(j)); %rate map histogram
        
        ylabel('Rate Map (Hz)')
        xlabel('Horizontal position (cm)')
        
        linkaxes([ax2 ax3 ax4],'x')
        %xlim([0 xmax])
        set(gcf, 'Position', [100 100 1536 1175]);
        saveas(gcf,[exp_directory filesep 'Analysis' filesep spike(j).name '-lap_' num2str(l) '.jpg'])
        set(h, 'Visible','off')
    end
end

%% CSC
csc_filename= fullfile(exp_directory,'Neuralynx','CSC4.ncs');
addpath('../jumping');
% theta phase
for j=1:length(spike([spike.m]==3))
    spike(j).phase = nan(size(spike(j).t));
end
for l=1:length(lap)
    % filtering out the leftward laps
    if lap(l).dir=="left" %strcmp(lap(l).dir,"left")
        continue;
    end
    figure(100+l)
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase] = filtertheta(timecsc,lfp);
    
    ax1 = subplot(5,1,1);
    plot(pos.t(pos.lap==l),pos.x(pos.lap==l));
    ylabel('Horizontal position (cm)')
    title(['lap ' num2str(l)]);
    
    ax2 = subplot(5,1,2);
    plot(pos.t(pos.lap==l),pos.s(pos.lap==l));
    ylabel('Speed (cm/s)')
    
    ax3 = subplot(5,1,3);
    plot(timecsc,phase)
    ylim([-180 180])
    ylabel('Phase (Degrees)')
    
    ax4 = subplot(5,1,4);
    plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
    hold on; 
    plot(timecsc,theta *1e6,'r');
    ylim([-300 300])
    ylabel('Theta (\muV)')
    
    % looking at all the clusters in m3 (m3 is the whole recorded experiment)
    for j=1:length(spike([spike.m]==3))
        if nnz([spike(j).lap]==l) > 0 % if there a firing for cluster j in this lap
            ax5 = subplot(5,1,5);
            hold off; hold on
            plot(spike(j).t([spike(j).lap]==l), spike(j).no,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
            % theta phase
            ax3 = subplot(5,1,3); hold on;
            spike(j).phase([spike(j).lap]==l) = interp1(timecsc,phase, spike(j).t([spike(j).lap]==l));
            if j == 1
                plot(spike(j).t([spike(j).lap]==l), spike(j).phase([spike(j).lap]==l),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
            end
        end
    end
    ax5 = subplot(5,1,5);
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
    
    % leftward rate map
    hist.posi = histcounts(posi.x(posi.dir=="left"), hist.edges) * dt; % seconds in each bin
    hist.cluster = histcounts(spike(j).x([spike(j).dir]=="left"), hist.edges); % spikes in each bin
    hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
    
    a(2*j-1) = subplot(length(spike([spike.m]==3)),2,2*length(spike([spike.m]==3))-(2*j-1));
    histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(j)); %rate map histogram
    ylabel(['cluster ' num2str(j)]);
    if j == length(spike([spike.m]==3))
        title('Rate Map (leftward)')
    elseif j ==1 
        xlabel('Horizontal position (cm)')
    end
    
    % rightward rate map
    hist.posi = histcounts(posi.x(posi.dir=="right"), hist.edges) * dt; % seconds in each bin
    hist.cluster = histcounts(spike(j).x([spike(j).dir]=="right"), hist.edges); % spikes in each bin
    hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
    
    a(2*j) = subplot(length(spike([spike.m]==3)),2,2*length(spike([spike.m]==3))-(2*j-1)+1);
    histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(j)); %rate map histogram
    ylabel(['cluster ' num2str(j)]);
    linkaxes([a(2*j-1) a(2*j)],'y')
    if j == length(spike([spike.m]==3))
        title('Rate Map (rightward)')
    elseif j ==1 
        xlabel('Horizontal position (cm)')
    end
    
    ylabel(['cluster ' num2str(j)]);
    zoom xon
end

set(gcf, 'Position', [100 100 1600 1300]);
linkaxes(a,'x')
xlim([0 xmax])
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
f = figure(202);
t1_count = 0;
t2_count = 0;
n = 10;
for l=2:2:30
    
    sp1.t = spike(1).t(spike(1).t >= lap(l).t(1)  & spike(1).t < lap(l).t(2) );
    sp3.t = spike(3).t(spike(3).t >= lap(l).t(1)  & spike(3).t < lap(l).t(2) );
    t_min = min(sp1.t(n),sp1.t(n)); % skip first n spikes
    t_max = max(sp1.t(end-n+1),sp1.t(end-n+1)); % skip last n spikes
    
    tjl = time_jump_rightward(l/2);
    t_max = tjl;
    t_mid = mean([t_min t_max]);
    for i=1:1
        
        % [t_min t_mid]
        ax1 = subplot(4,8,l-1);
        sp1.t = spike(1).t(spike(1).t >= t_min & spike(1).t < t_mid );
        sp3.t = spike(3).t(spike(3).t >= t_min & spike(3).t < t_mid );

        [xx,yy] = meshgrid(sp1.t, sp3.t);
        zz = xx - yy;
        range = 0.05;
        %zz((zz>range) | (zz<-range)) = nan;
        h = histogram(zz,(-range:0.003:range),'FaceColor','b');
        t1_count = t1_count + h.Values;
        
        title(['lap ' num2str(l) ' t_min to t_mid'], 'Interpreter', 'none');
        ylabel('No. Incidents')
        xlabel('Time (sec)')
        
        % [t_mid t_max]
        ax2 = subplot(4,8,l);
        sp1.t = spike(1).t(spike(1).t >= t_mid & spike(1).t < t_max );
        sp3.t = spike(3).t(spike(3).t >= t_mid & spike(3).t < t_max );

        [xx,yy] = meshgrid(sp1.t, sp3.t);
        zz = xx - yy;
        %zz((zz>range) | (zz<-range)) = nan;
        h = histogram(zz,(-range:0.003:range),'FaceColor','r');
        t2_count = t2_count + h.Values;
        
        xlim([-range range])
        title(['lap ' num2str(l) ' t_mid to t_max'], 'Interpreter', 'none');
        ylabel('No. Incidents')
        xlabel('Time (sec)')
        linkaxes([ax1 ax2],'y')
    end
end
f.WindowState = 'maximized';
sgtitle('Interspike interval between clusters 1 and 3');

f = figure(203);
ax1 = subplot(1,2,1);
histogram('BinEdges',h.BinEdges,'BinCounts',t1_count,'FaceColor','b');
xlim([-range range])
title('All laps compined from t_min to t_mid', 'Interpreter', 'none');
ylabel('No. Incidents')
xlabel('Time (sec)')

ax2 = subplot(1,2,2);
histogram('BinEdges',h.BinEdges,'BinCounts',t2_count,'FaceColor','r');
xlim([-range range])
title('All laps compined from t_mid to t_max', 'Interpreter', 'none');
ylabel('No. Incidents')
xlabel('Time (sec)')

linkaxes([ax1 ax2],'y')
f.WindowState = 'maximized';
sgtitle('Interspike interval between clusters 1 and 3');

%% Self-Interspike interval (test)
% f = figure(204);
% t1_count = 0;
% t2_count = 0;
% for l=2:2:30
%     
%     sp1.t = spike(1).t(spike(1).t >= lap(l).t(1)  & spike(1).t < lap(l).t(2) );
%     sp3.t = spike(3).t(spike(3).t >= lap(l).t(1)  & spike(3).t < lap(l).t(2) );
%     t_min = min(sp1.t(10),sp1.t(10));
%     t_max = max(sp1.t(end-9),sp1.t(end-9));
%     t_mid = mean([t_min t_max]);
%     %tjl = time_jump_rightward(l/2)
%     for i=1:1
%         
%         % [t_min t_mid]
%         subplot(4,8,l-1)
%         sp1.t = spike(1).t(spike(1).t >= t_min & spike(1).t < t_mid );
%         sp3.t = spike(3).t(spike(3).t >= t_min & spike(3).t < t_mid );
% 
%         [xx,yy] = meshgrid(sp1.t, sp1.t);
%         zz = xx - yy;
%         range = 0.05;
%         zz((zz>range) | (zz<-range)) = nan;
%         h = histogram(zz,(-range:0.005:range),'FaceColor','b');
%         t1_count = t1_count + h.Values;
%         
%         title(['lap ' num2str(l) ' t_min to t_mid'], 'Interpreter', 'none');
%         ylabel('No. Incidents')
%         xlabel('Time (sec)')
%         
%         % [t_mid t_max]
%         subplot(4,8,l)
%         sp1.t = spike(1).t(spike(1).t >= t_mid & spike(1).t < t_max );
%         sp3.t = spike(3).t(spike(3).t >= t_mid & spike(3).t < t_max );
% 
%         [xx,yy] = meshgrid(sp1.t, sp1.t);
%         zz = xx - yy;
%         zz((zz>range) | (zz<-range)) = nan;
%         h = histogram(zz,(-range:0.005:range),'FaceColor','r');
%         t2_count = t2_count + h.Values;
%         
%         xlim([-range range])
%         title(['lap ' num2str(l) ' t_mid to t_max'], 'Interpreter', 'none');
%         ylabel('No. Incidents')
%         xlabel('Time (sec)')
%     end
% end
% f.WindowState = 'maximized';
% sgtitle('Interspike interval between clusters 1 and 1');
% 
% f = figure(205);
% subplot(1,2,1)
% histogram('BinEdges',h.BinEdges,'BinCounts',t1_count,'FaceColor','b');
% xlim([-range range])
% title('All laps compined from t_min to t_mid', 'Interpreter', 'none');
% ylabel('No. Incidents')
% xlabel('Time (sec)')
% 
% subplot(1,2,2)
% histogram('BinEdges',h.BinEdges,'BinCounts',t2_count,'FaceColor','r');
% xlim([-range range])
% title('All laps compined from t_mid to t_max', 'Interpreter', 'none');
% ylabel('No. Incidents')
% xlabel('Time (sec)')
% 
% f.WindowState = 'maximized';
% sgtitle('Interspike interval between clusters 1 and 1');

%% Saving data
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
save(mat_filename,'pos','posi', 'lap', 'spike','ppcm', 'offset', 'colors','xmax','x_thresh');
disp(['File ' mat_filename ' has been created!'])