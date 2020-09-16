clc; close all;
if ~exist('pos_clust','var')
    loaddata;
end

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
maze = 2;

TimeEV = readevent;
StartTime = TimeEV(1:2:end);
EndTime = TimeEV(2:2:end);

%% Pos.p
[TimePos,x,y] = readposp(fullfile(Nlx_directory,'pos.p'), StartTime(maze), EndTime(maze));
y = 480 - y;
AngularPosition = rad2deg(atan2(y-240, x-320));
AngularPosition = wrapTo360(AngularPosition);

%% lap detector
lap = 3;

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
tiledlayout(4,1)
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



%% Histogram .. 


frameRate = 30; %frame rate in Hz
edges = linspace(0, 360, 360);
Edges = movmean(edges,2);
t_spike = seconds(t);
HistOcc = histcounts(AngularPositionF, edges);


%% interpolation
f = 3000; % interpolation 3 kHz
TimePosInterp = TimePos(1):1/f: TimePos(end);

xInterp = interp1(TimePos,x,TimePosInterp);
yInterp = interp1(TimePos,y,TimePosInterp);
AngularPositionInterp = rad2deg(atan2(yInterp-240, xInterp-320));
AngularPositionInterp = wrapTo360(AngularPositionInterp);

%% histogram
edges = linspace(0, 360, 360);
HistOcc = histcounts(AngularPosition, edges);
HistOccInterp = histcounts(AngularPositionInterp, edges);

f3 = figure(3);
rectangle('Position',[160,0,40,3.5e4],'FaceColor','g','EdgeColor','g')
hold on
rectangle('Position',[100,0,10,3.5e4],'FaceColor','c','EdgeColor','c')
hold on
rectangle('Position',[260,0,10,3.5e4],'FaceColor','c','EdgeColor','c')
hold on
bar(HistOccInterp);
title('Occupancy (Interpolated)')

%% velocity
% 5 ft = 1524 mm = 480 pixels
% each pixel = 3.175 mm

vx = gradient(x)./gradient(TimePos); % pixels/sec
vy = gradient(y)./gradient(TimePos); % pixels/sec

velocity = sqrt(vx.^2+vy.^2); % pixels/sec
velocity = velocity * 0.3175; % cm/sec

velocityF = gradient(AngularPositionF)./gradient(TimePosF); % deg/sec (1 deg ~ 1.3 cm)

%% rate map of cells (histogram of occupancy corrected, and velocity)
E=zeros(17,1);
RateMap = zeros(17,length(edges)-1);
for i=1:17
    t_spike_cluster = t_spike(array_combined==i);
    SpikePos = interp1(TimePosF,AngularPositionF,t_spike_cluster);
    SpikeVelocity = interp1(TimePos,velocity,t_spike_cluster);
    
    HistClust = histcounts(SpikePos(SpikeVelocity>5), edges);
    
    RateMap(i,:) = (HistClust./HistOccInterp)*frameRate;
    
    [~, argmaxRateMap] = max(RateMap(i,:));
    E(i)=Edges(argmaxRateMap);
end

f4 = figure(4);
rectangle('Position',[160,0,40,1],'FaceColor','g','EdgeColor','g')
hold on
rectangle('Position',[100,0,10,1],'FaceColor','c','EdgeColor','c')
hold on
rectangle('Position',[260,0,10,1],'FaceColor','c','EdgeColor','c')
hold on
bar(sum(RateMap));
title('Rate map for all clusters (Interpolated)')

[~, SortedIndex] = sort(E);
k=0;

f2 = figure(2);
for i=1:17
    k = SortedIndex(i);
    subplot(17,1,i);
    rectangle('Position',[160,0,40,1],'FaceColor','g','EdgeColor','g')
    hold on
    rectangle('Position',[100,0,10,1],'FaceColor','c','EdgeColor','c')
    hold on
    rectangle('Position',[260,0,10,1],'FaceColor','c','EdgeColor','c')
    hold on
    bar(Edges(2:end),RateMap(k,:));
    xlim([0 360])
    ylim([0 0.35])
    ylabel(['cell no ' num2str(k)])
    if i==1
        title('Rate Map (Interpolated)')
    end
end

%%
f1 = figure(1);
ax2 = nexttile;

rasterplot(t_spike(3:end),SortedIndex(array_combined(3:end)));
pan xon
zoom on
ylabel('Spike raster plot (Sorted)')

ax3 = nexttile;
plot(TimeCSC, lfp*1e3, 'b');%'Color', uint8([230 230 230])); % milliVolts
hold on
plot(TimeCSC, theta*1e3,'r'); % milliVolts
ylim([-1.2 1.2])
ylabel('LFP (\muV)')
pan xon
zoom xon
linkaxes([ax1 ax2 ax3],'x')

% velocity
CSCVelocity = interp1(TimePos,velocity,TimeCSC);
ax4 = nexttile;
rectangle('Position',[TimeRange(1),0,diff(TimeRange),5],'FaceColor','y','EdgeColor','y')
hold on
plot(TimeCSC, CSCVelocity,TimePosF(2:end-1),velocityF(2:end-1));
xlabel('Time (sec)')
ylabel('Velocity (cm/s)')
zoom on
ylim([0 150])
linkaxes([ax1 ax2 ax3 ax4],'x')

ax1.XLim = TimeRange + [35e-3 -15e-3];
%ax1.XLim = [5770 5990];
%ax1.XLim = [5780 5786];
% ax1.XLim = [5820 5826]; % m4l4
ax1.Title.String = ['maze ' num2str(maze) ', lap ' num2str(lap)];
% f1.WindowState = 'maximized';
% f2.WindowState = 'maximized';
% f3.WindowState = 'maximized';
% f4.WindowState = 'maximized';