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
pos_clust_table = timetable2table(pos_clust);
pos_clust_table.time = seconds(pos_clust_table.time);
pos_clust_array = table2array(pos_clust_table);

%selecting a portion of the cells and plotting
close all
start_t = seconds(sort(unique(start_time))*1e-6)';
end_t = seconds(sort(unique(end_time))*1e-6)';
end_t(3) = [];

m = 4;
S = timerange(start_t(m),end_t(m));
pos_col = [false(2,1); true];
clust_col = (maze_no == m) & (tt_no < 18);
maze = pos_clust(S, [pos_col; clust_col]);
maze_table = timetable2table(maze);
maze_table.time = seconds(maze.time);
maze_array = table2array(maze_table);
no_cells = size(maze,2)-1;

all_cells = repmat(1:no_cells,size(maze,1),1);
f = maze_array(:,3:end);
g = f.*all_cells;
h=~isnan(g);
c1 = find(h);
[r,c] = find(h);

array_combined = c; %g(c2)

array_combined(array_combined==0) = NaN;
array_combined(1)=0;
array_combined(2)=no_cells+1;
array_combined = array_combined -0.25 + 0.5*rand(length(array_combined),1); % adding random jitter

t = seconds(maze_array(r,1));
th = maze_array(r,2);

TT1 = timetable(maze.time,maze.theta,'VariableNames',{'Angular_Position'});
TT2 = timetable(t,array_combined,'VariableNames',{'Cell_No'});

[TimeSec,Data] = readcsc;
TT3 = timetable(seconds(TimeSec),Data,'VariableNames',{'CSC'});
TT3 = TT3(S,:);

TT = synchronize(TT1,TT2,TT3);

figure(1)
stackedplot(TT,'.');
title(['maze ' num2str(m)])
 
head(TT);
tail(TT);

figure(2)
yyaxis left
count = histogram(t,'BinWidth',milliseconds(50));

ylabel('Frequency')
yyaxis right
plot(maze.time,maze.theta,'.')
ylim([0 360])
ylabel('Angular Position')
title(['maze ' num2str(m)])

%% 

% figure(3)
% tiledlayout(3,1)
% 
% ax1 = nexttile;
% rectangle('Position',[seconds(start_t(m)),160,seconds(end_t(m)-start_t(m)),40],'FaceColor','g','EdgeColor','g')
% hold on
% plot(seconds(maze.time),maze.theta,'*')
% ylim([0 360])
% zoom xon
% pan xon
% ylabel('Angular Position (\circ)')
% title(['maze ' num2str(m)])
% 
% timef = maze.time(~isnan(maze.theta));
% thetaf = maze.theta(~isnan(maze.theta));
% lap_start = timef(diff(thetaf) < -350);
% 
% ax2 = nexttile;
% [counts5,edges] = histcounts(t,'BinWidth',milliseconds(5)); % 5 milliseconds
% edges50 = movmean(seconds(edges),2);
% counts50 = movsum(counts5,10);
% bar(edges50(2:end),counts50)
% zoom xon
% pan xon
% ylabel('Frequency')
% 
% ax3 = nexttile;
% plot(TimeSec,Data*1e6);
% ylim([-1200 1200])
% xlabel('Time (sec)')
% ylabel('CSC (\muV)')
% zoom xon
% pan xon
% linkaxes([ax1 ax2 ax3],'x')
% %%
% lap = 2;
% ax1.XLim = seconds([lap_start(lap) lap_start(lap+1)]);
% ax1.Title.String = ['maze ' num2str(m) ', lap ' num2str(lap)];