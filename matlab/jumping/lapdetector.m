function [lap, idx_analysis] = lapdetector(exp_directory, varargin)
%%LAPDETECTOR Lap Detector gets a data file (data.mat) that is generated by winclust2mat
% and detects the lap info from it.
%
%  Inputs:
%       (1) Experiment Directory
%       (2) Length of image (in cm)
%       (3) Mode of Lap detection: 
%           (a) mode = [] empty: the whole lap; 
%           (b) mode = T : T sec before to T sec after
%           (c) mode = [T1 T2] : T1 se before up to T2 
%               sec after. sign doesn't matter: 
%               mode = [-abs(mode(1)) abs(mode(end))])
%       (4) Time Analysis 2xN [Initial Final]
%               e.g. [1 2; 3 4; 5 6; 7 8]
%       (5) 'exclude' Exclusion times 2xM [Initial Fianal]
%
%   Outputs
%       (1) lap: A structure that includes
%               dir: direction "left"/"right"
%               no: number
%               t: time of lap [time_initial time_final]
%               t_jump: time of jump
%               frame: frame number for the jump
%       (2) idx_analysis: all the indeces during the analysis time
%
%   Examples
%       To run the function with any pop-up window:
%           lap = lapdetector(exp_directory, xmax)
%
%       To specify the time of analysis:
%           lap = lapdetector(exp_directory, xmax, mode, time_of_analysis)
%
%       All inputs and outputs:
%           [lap, idx_analysis] = lapdetector(exp_directory, xmax, mode, time_of_analysis)
%
%   See also WINCLUST2MAT, ANALYZEDATA.
%
%   SGL 2022-02-09 (originally 2021-01-31, 2021-03-28)
%
clc; close all

%% Variable-length input argument parsing
if nargin < 1
    [~,exp_directory] = uigetfile(fullfile('D:\Analysis', 'data.mat'), 'Select Data File');
end

mat_filename = fullfile(exp_directory, 'data.mat');
if isfile(mat_filename)
    clearvars -except mat_filename exp_directory varargin;
    load(mat_filename, 'pos', 'exp', 'ppcm','daq')
else
    error('Data file is not valid!')
end

xmax = 2048/ppcm;
mode = [-inf inf]; % the whole lap
time_of_analysis = [pos.t(1) pos.t(end)];
time_exclusion = [];
%x_thresh = [133 143];
x_thresh = [133 139];

for argidx = 1:2:nargin-1
    switch varargin{argidx}
        case 'xmax'
            xmax = varargin{argidx+1};
        case 'mode'
            mode = varargin{argidx+1};
        case 'time'
            time_of_analysis = varargin{argidx+1};
        case 'exclude'
            time_exclusion = varargin{argidx+1};
        case 'thresh'
            x_thresh = varargin{argidx+1};
    end
end
mode = [-abs(mode(1)) abs(mode(end))]; % e.g. 2 would be [-2 -2] from 2 seconds before to 2 seconds after
tai = time_of_analysis(:,1);
taf = time_of_analysis(:,2);

%% Plotting the position vs time and marking neural recording times,

%% time detection: recording, maze, and analysis

% time of recording
tri = 0;                    % time of recording (initial)
trf = seconds(exp.finish-exp.start); % time of recording (final)

hold on
for i=1:length(tri)
    rectangle('Position',[tri(i),0,trf(i)-tri(i),xmax],'FaceColor','y')
end

% time of maze
tmi = pos.t(1); % tri % time of maze initial
tmf = pos.t(end); % trf % time of maze final

for i=1:length(tmi)
    rectangle('Position',[tmi(i),0,tmf(i)-tmi(i),xmax],'FaceColor',[0 .8 .8])
end

% extract data only during the time of analysis
idx_analysis = zeros(length(pos.t),1);
for i=1:length(tai)
    idx_analysis = idx_analysis | ((pos.t > tai(i)) & (pos.t < taf(i)));
end
% excluding some outliers
for i=1:size(time_exclusion,1)
    idx_analysis = idx_analysis & ~((pos.t > time_exclusion(i,1)) & (pos.t < time_exclusion(i,2)));
end

t = pos.t(idx_analysis);
x = pos.x(idx_analysis);
vx = pos.vx(idx_analysis);
vx = filterlfp(t, vx, 0.01, 2); % cm/sec
    
plot(pos.t, pos.x, '.b', t, x, '.k')
ylim([0 xmax])
plot(t,abs(vx))
%% jump detection

daq.ditch = double(daq.loadcell(2,:)>1);
ditch = interp1(daq.t,daq.ditch,t,'linear','extrap');
ditch = movmax(ditch, [100 0]); % extend the ditch
plot(t,ditch*xmax)
% times that the rats jump (based on jump direction)
jump_criteria_leftward = (x > x_thresh(1)) & (x < x_thresh(2)) & (vx < -40) & (~ditch); % multiple criterion
jump_criteria_rightward = (x > x_thresh(1)) & (x < x_thresh(2)) & (vx > 40) & (~ditch); % multiple criterion

% detecting the first of such chage
jump_criteria_leftward = diff(jump_criteria_leftward) == 1;
jump_criteria_rightward = diff(jump_criteria_rightward) == 1;

% removing double detection
jump_criteria_leftward_sum = movsum(jump_criteria_leftward, [0 500]);
jump_criteria_rightward_sum = movsum(jump_criteria_leftward, [0 500]);
jump_criteria_leftward(jump_criteria_leftward_sum > 1) = 0;
jump_criteria_rightward(jump_criteria_rightward_sum > 1) = 0;

time_jump_leftward = t(jump_criteria_leftward);
time_jump_rightward = t(jump_criteria_rightward);
x_jump_leftward = x(jump_criteria_leftward);
x_jump_rightward = x(jump_criteria_rightward);

plot(time_jump_leftward, x_jump_leftward, 'hk', 'MarkerSize',15)
plot(time_jump_rightward, x_jump_rightward, 'hr', 'MarkerSize',15)

%% detecting beginning and end of laps (extremes)
time_jump = sort([time_jump_leftward;time_jump_rightward; tai; taf;taf(end)]); % exception handling for start of stop of analysis
l = 1;
N = length(time_jump);

for i = 1:1:N-2
    %lap extremes (max and min)
    jump_idx_lap = t>=time_jump(i) & t<time_jump(i+1);
    tl = t(jump_idx_lap);
    [x_right(l), idx] = max(x(jump_idx_lap));
    t_right(l) = tl(idx);
    [x_left(l), idx] = min(x(jump_idx_lap));
    t_left(l) = tl(idx);
    l = l+1;
end

% ignore extremes in the vicinity of the gap
t_right(x_right < x_thresh(2) + 40) = [];
x_right(x_right < x_thresh(2) + 40) = [];
t_left(x_left > x_thresh(1) - 40) = [];
x_left(x_left > x_thresh(1) - 40) = [];

plot(t_right, x_right, 'pk', 'MarkerSize',15)
plot(t_left, x_left, 'sk', 'MarkerSize',15)

%% lap detection
time_jump = sort([time_jump_leftward;time_jump_rightward]); % no exception handling
time_extreme = sort([t_right t_left]);
N = length(time_extreme);
l = 1;

for i=1:N-1
    % leftward laps
    if ismember(time_extreme(i),t_right) && ismember(time_extreme(i+1),t_left) ...
            && (nnz(time_jump >= time_extreme(i) & time_jump <= time_extreme(i+1))>0)
        lap(l).dir = "left";
        lap(l).no = l;
        % time of jump in the lap
        lap(l).t_jump = time_jump(time_jump >= time_extreme(i) & time_jump <= time_extreme(i+1));
        lap(l).frame = pos.frame(pos.t==lap(l).t_jump); % frame of jump in the lap
        % time of lap
        % [time_extreme(i) time_extreme(i+1)]; % from min to max or vice versa
        % [lap(l).t_jump-n lap(l).t_jump+m]; % from n sec before to m sec after
        lap(l).t = [max(time_extreme(i),lap(l).t_jump + mode(1)) min(time_extreme(i+1),lap(l).t_jump + mode(2))];
        l = l+1;
    end
    % rightward laps
    if ismember(time_extreme(i),t_left) && ismember(time_extreme(i+1),t_right) ...
            && (nnz(time_jump >= time_extreme(i) & time_jump <= time_extreme(i+1))>0)
        lap(l).dir = "right";
        lap(l).no = l;
        % time of jump in the lap
        lap(l).t_jump = time_jump(time_jump >= time_extreme(i) & time_jump <= time_extreme(i+1));
        lap(l).frame = pos.frame(pos.t==lap(l).t_jump); % frame no of jump in the lap
        % time of lap
        lap(l).t = [max(time_extreme(i),lap(l).t_jump + mode(1)) min(time_extreme(i+1),lap(l).t_jump + mode(2))];
        l = l+1;
    end
end

%% plotting each lap seperately IN ALTERNATING COLORS (first one is red)
for l=1:length(lap)
    idx = t >= lap(l).t(1) & t <= lap(l).t(2);
    if lap(l).dir == "right"
        color = '.g';
    else
        color = '.m';
    end
    if l == 1
        color = '.r';
    end
    plot(t(idx),x(idx),color)
end
xlim([tri trf])
set(gcf, 'Position', [50 50 2300 1300]);

%% Marking laps as jump or ditch
for l=1:length(lap)
    idx = t >= lap(l).t_jump+0.4 & t <= lap(l).t_jump+0.6;
    if mean(ditch(idx))>0.9 % 90% probability ditch
        lap(l).status = "ditch";
    else
        lap(l).status = "jump";
    end 
end
    
%% save file if the functiona in called, otherwise display the info
if nargout ~= 0
    analysis_directory = fullfile(exp_directory, 'Analysis');
    if ~exist(analysis_directory, 'dir')
        mkdir(analysis_directory)
    end
    saveas(gcf,fullfile(analysis_directory,'overall_view.svg'));
else
    format longG
    disp(exp.name);
    disp('Frame numbers: ')
    for l=1:length(lap)
        fprintf('%d,',lap(l).frame);
    end
    fprintf('\n')
    fprintf('Maze Times: \n')
    disp([tmi tmf])
    fprintf('Analysis Times: \n')
    disp([tai taf])
end

end