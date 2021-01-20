% This program plots the data for Rat 913
% The data include the occupancy and histogram of the spikes and the
% sequential firing of multiple spikes before jumping
% SGL 2020-11-28
clc; close all;
exp_directory = '~/Desktop/20-12-09';
[datafile,exp_directory] = uigetfile(fullfile(exp_directory,'data.mat'), 'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename, 'pos', 'spike', 'ppcm', 'offset')
clearvars -except datafile exp_directory pos spike ppcm offset;
colors = ["#EDB120" "#7E2F8E" "yellow" "#A2142F" "red" "magenta" "green" "#D95319"];
colors = repmat(colors', ceil(length(spike)/length(colors))); % repeat the colors to match the total number of spikes
[img_file,img_directory] = uigetfile(fullfile(exp_directory,'Videos','*.pgm'), 'Select Image File');
img_filename = fullfile(img_directory, img_file);
I = imread(img_filename);
xmax = ceil(size(I,2)/ppcm);
spike = spike([spike.m]==3); % only looking at the clusters in m3 (m3 is the whole recorded experiment)

%% Plotting the position vs time and marking neural recording times, 
% maze times and lap detection times
close all
event_filename = fullfile(exp_directory, 'Neuralynx', 'Events.nev');
[Time,Data,Header,EventIDs,TTLs] = readevent(event_filename)

idx = string(Data)=='Starting Recording';
ri = Time(idx);
idx = string(Data)=='Stopping Recording';
rf = Time(idx);

% neural recording times
r1i = ri(1); %105.221068;
r1f = rf(1); %703.846944;
r2i = ri(2); %945.026778;
r2f = rf(2); %1341.248894;

% maze times
m1i = 408.650168;
m1f = 677.921799;
m2i = 988.727761;
m2f = 1303.765711;

% analysis times
t1i = r1i;
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

% ignore lap 25
lap(25) = [];

hold on
rectangle('Position',[r1i,0,r1f-r1i,xmax],'FaceColor','y')
rectangle('Position',[r2i,0,r2f-r2i,xmax],'FaceColor','y')

plot(pos.t, pos.x, '.b', t, x, '.k')
ylim([0 xmax])

plot(t_left, x_left, 'pk', 'MarkerSize',15)
plot(t_right, x_right, 'sk', 'MarkerSize',15)

plot(time_jump_leftward, x_jump_leftward, 'hc', 'MarkerSize',15)
plot(time_jump_rightward, x_jump_rightward, 'hr', 'MarkerSize',15)

set(gcf, 'Position', [100 100 2000 1500]);

%%
%{
%% Adding laps and direction to spike for m3 maze
for j=1:length(spike)
    spike(j).lap = zeros(size(spike(j).t));
    spike(j).dir = strings(size(spike(j).t));
    for l = [lap.no]
        idx = (spike(j).t >=lap(l).t(1)) & (spike(j).t <lap(l).t(2)); % lap
        spike(j).lap(idx)= lap(l).no;
        spike(j).dir(idx)=lap(l).dir;
    end
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
for l = [lap.no]
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

%% CSC
csc_filename= fullfile(exp_directory,'Neuralynx','CSC4.ncs');
addpath('../jumping');
% theta phase
for j=1:length(spike)
    spike(j).phase = nan(size(spike(j).t));
end
for l=[lap([lap.dir]=="left").no] % all leftward laps

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
    for j=1:length(spike)
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
    ylim([0 length(spike)+1])
    ylabel('Cluster Number')
    
    linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
    xlim(lap(l).t);
    
    set(gcf, 'Position', [100 100 1536 1175]);
    saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-lap_' num2str(l) '.jpg']))
    
end

%% Directional rate map for all the laps
figure(200)

for j=1:length(spike)
    
    % leftward rate map
    hist.posi = histcounts(posi.x(posi.dir=="left"), hist.edges) * dt; % seconds in each bin
    hist.cluster = histcounts(spike(j).x([spike(j).dir]=="left"), hist.edges); % spikes in each bin
    hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
    
    a(2*j-1) = subplot(length(spike),2,2*length(spike)-(2*j-1));
    histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(j)); %rate map histogram
    ylabel(['cluster ' num2str(j)]);
    if j == length(spike)
        title('Rate Map (leftward)')
    elseif j ==1 
        xlabel('Horizontal position (cm)')
    end
    
    % rightward rate map
    hist.posi = histcounts(posi.x(posi.dir=="right"), hist.edges) * dt; % seconds in each bin
    hist.cluster = histcounts(spike(j).x([spike(j).dir]=="right"), hist.edges); % spikes in each bin
    hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
    
    a(2*j) = subplot(length(spike),2,2*length(spike)-(2*j-1)+1);
    histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(j)); %rate map histogram
    ylabel(['cluster ' num2str(j)]);
    linkaxes([a(2*j-1) a(2*j)],'y')
    if j == length(spike)
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

%% Saving data
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
save(mat_filename,'pos','posi','spike','ppcm', 'offset');
disp(['File ' mat_filename ' has been created!'])

%}