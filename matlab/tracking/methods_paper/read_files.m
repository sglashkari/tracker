clc; close all;
% Cluster data
% file_name = fullfile('C:\Users\dome3neuralynx\OneDrive - Johns Hopkins\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx\TT11','cl-maze1.1');
% A = 
% time = A.data(:,18);
Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
listing = dir(fullfile(Nlx_directory,'**','cl-maze*.*'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing); % number of maze-clusters

tt_no = str2double(extractAfter(folders,"\TT"));
maze_no = str2double(extractBetween(names,"cl-maze","."));
cluster_no = str2double(extractAfter(names,"."));

% T = table(maze_no, tt_no, cluster_no, 'VariableNames',{'maze','TT','cluster'});
absolue_paths = mat2cell(folders+'\'+names,ones(N,1));

A = cellfun(@(x) importdata(x,',',13), absolue_paths);

%% Postion data
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

for index = 1:N
    start_time(index) = str2double(A(index).textdata{12});
    end_time(index) = str2double(A(index).textdata{13});
    %ts = start_time:1/30000:end_time;
    time = seconds((A(index).data(:,18))*1e-6);
    all_ones = ones(length(time),1);
    
    maze_tt_cl = cellstr(strcat('m',string(maze_no),'_tt',string(tt_no),'_c',string(cluster_no)));
    cluster_maze = timetable(time,all_ones,'VariableNames',maze_tt_cl(index));
    
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

%% selecting a portion of the cells and plotting

start_t = seconds(sort(unique(start_time))*1e-6)';
end_t = seconds(sort(unique(end_time))*1e-6)';

%maze = cell(5,1);
for m = 1:5
S = timerange(start_t(m),end_t(m));
pos_col = [false(2,1); true];
clust_col = (maze_no == m) & (tt_no < 11);
maze = pos_clust(S, [pos_col; clust_col]);

head(maze)
figure(m);
subplot(3,1,1)
stackedplot(maze,{'theta'},'.')
subplot(3,1,[2 3])
stackedplot(maze,'.')
title(['maze ' num2str(m)])
end

% head(maze_2)
% tail(maze_2)
