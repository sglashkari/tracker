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