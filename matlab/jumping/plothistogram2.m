clc; close all;
%% Cluster data
% file_name = fullfile('C:\Users\dome3neuralynx\OneDrive - Johns Hopkins\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx\TT11','cl-maze1.1');
% A = 
% time = A.data(:,18);
clear
Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
listing = dir(fullfile(Nlx_directory,'**','cl-maze*.*'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing); % number of maze-clusters

tt_no = str2double(extractAfter(folders,"\TT"));
maze_no = floor(str2double(extractAfter(names,"cl-maze")));
cluster_no = str2double(extractAfter(names,"."));
cluster_no(isnan(cluster_no))=0;    % cluster 0

% T = table(maze_no, tt_no, cluster_no, 'VariableNames',{'maze','TT','cluster'});
absolue_paths = mat2cell(folders+'\'+names,ones(N,1));

A = cellfun(@(x) importdata(x,',',13), absolue_paths); %, 'UniformOutput', false);

% Postion data
pos_p_file_name = fullfile(Nlx_directory,'Pos.p.ascii');
B = importdata(pos_p_file_name,',',24);
time = seconds(B.data(:,1)*1e-6);
x = B.data(:,2);
y = 480 - B.data(:,3);
theta = rad2deg(atan2(y-240, x-320));
theta = wrapTo360(theta);
position = timetable(time,x,y,theta);

%% Generate a timetable (named all_cluster_mazes) with all the timestamps
start_time = zeros(N,1);
end_time = zeros(N,1);

maze_tt_cl = cellstr(strcat('m',string(maze_no),'_tt',string(tt_no),'_c',string(cluster_no)));

for index = 1:N
    start_time(index) = str2double(A(index).textdata{12});
    end_time(index) = str2double(A(index).textdata{13});
    %ts = start_time:1/30000:end_time;
    time = seconds((A(index).data(:,18))*1e-6);
    all_ones = ones(length(time),1);
    
    cluster_maze = timetable(time,all_ones,'VariableNames',maze_tt_cl(index));
    index
    if index == 1
        cluster_mazes = cluster_maze;
    else
        cluster_mazes = synchronize(cluster_mazes,cluster_maze);
    end
end

%% combining position and spikes
pos_clust = synchronize(position,cluster_mazes);
% pos_clust_table = timetable2table(pos_clust);
% pos_clust_table.time = seconds(pos_clust_table.time);
% pos_clust_array = table2array(pos_clust_table);

%% selecting a portion of the cells and plotting
close all


%%
% ...
%
%
Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
listing = dir(fullfile(Nlx_directory,'**','TT*.ntt'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing);
absolue_paths = folders+'\'+names;

%% maze
maze = 4;

TimeEV = readevent;
StartTime = TimeEV(1:2:end);
EndTime = TimeEV(2:2:end);

%% Pos.p
[TimePos,x,y] = readposp;
y = 480 - y;
AngularPosition = rad2deg(atan2(y-240, x-320));
AngularPosition = wrapTo360(AngularPosition);

%% lap detector
lap = 5;

TimePosF = TimePos(~isnan(AngularPosition));
AngularPositionF = AngularPosition(~isnan(AngularPosition));
AngularPositionF = AngularPositionF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));
TimePosF =  TimePosF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));

LapTime = TimePosF(diff(AngularPositionF) < -340);
TimeRange = [LapTime(lap) LapTime(lap+1)]; % seconds


%% CSC
FilenameCSC = fullfile(Nlx_directory,'CSC14.ncs');
[TimeCSC,lfp] = readcsc(FilenameCSC, TimeRange * 1e6); % microseconds
theta = filtertheta(TimeCSC,lfp);
plot(TimeCSC,theta,'r')

%% selecting a portion of the cells and plotting
close all
start_t = seconds(sort(unique(start_time))*1e-6)';
end_t = seconds(sort(unique(end_time))*1e-6)';
end_t(3) = [];

S = timerange(start_t(maze),end_t(maze));
pos_col = [false(2,1); true];
clust_col = (maze_no == maze) & (tt_no < 18);
maze1 = pos_clust(S, [pos_col; clust_col]);
maze_table = timetable2table(maze1);
maze_table.time = seconds(maze1.time);
maze_array = table2array(maze_table);
no_cells = size(maze1,2)-1;

all_cells = repmat(1:no_cells,size(maze1,1),1);
f = maze_array(:,3:end);
g = f.*all_cells;
h=~isnan(g);
c1 = find(h);
[r,c] = find(h);

array_combined = c; %g(c2)

array_combined(array_combined==0) = NaN;
array_combined(1)=0;
array_combined(2)=no_cells+1;

t = seconds(maze_array(r,1));

%% selecting a portion

%% Plotting
f = figure(1);
tiledlayout(3,1)
ax1 = nexttile;
rectangle('Position',[TimeRange(1),160,diff(TimeRange),40],'FaceColor','g','EdgeColor','g')
hold on
rectangle('Position',[TimeRange(1),100,diff(TimeRange),10],'FaceColor','c','EdgeColor','c')
hold on
rectangle('Position',[TimeRange(1),260,diff(TimeRange),10],'FaceColor','c','EdgeColor','c')
hold on
plot(TimePosF,AngularPositionF,'*')
ylim([0 360])
pan xon
zoom xon
ylabel('Angular Position (\circ)')
title(['maze ' num2str(maze)])



%%


frameRate = 30; %frame rate in Hz
edges = linspace(0, 360, 360);
Edges = movmean(edges,2);
t_spike = seconds(t);
HistOcc = histcounts(AngularPositionF, edges);

f3 = figure(3);
bar(HistOcc);

f2 = figure(2);
E=zeros(17,1);
RateMap = zeros(17,length(edges)-1);
for i=1:17
    t_spike_lap = t_spike(array_combined==i);
    SpikePos = interp1(TimePosF,AngularPositionF,t_spike_lap);

    HistClust = histcounts(SpikePos, edges);
    
    RateMap(i,:) = (HistClust./HistOcc)*frameRate;
    
    [~, argmaxRateMap] = max(RateMap(i,:));
    E(i)=Edges(argmaxRateMap);
end
[~, SortedIndex] = sort(E);
k=0;
for i=1:17
    k = SortedIndex(i);
    subplot(17,1,i);
    bar(Edges(2:end),RateMap(k,:));
    ylabel(['cell no ' num2str(k)])
end

%%
figure(1);
ax2 = nexttile;

rasterplot(t_spike(3:end),SortedIndex(array_combined(3:end)));
pan xon
zoom on
ylabel('Spike raster plot (Sorted)')

ax3 = nexttile;
plot(TimeCSC, lfp*1e3, 'Color', uint8([230 230 230])); % milliVolts
hold on
plot(TimeCSC, theta*1e3,'r'); % milliVolts
ylim([-0.5 0.5])
xlabel('Time (sec)')
ylabel('LFP (\muV)')
pan xon
zoom xon
linkaxes([ax1 ax2 ax3],'x')

ax1.XLim = TimeRange + [35e-3 -15e-3];
%ax1.XLim = [5780 5786];
% ax1.XLim = [5820 5826]; % m4l4
ax1.Title.String = ['maze ' num2str(maze) ', lap ' num2str(lap)];
f.WindowState = 'maximized';