%% Cluster data
% file_name = fullfile('C:\Users\dome3neuralynx\OneDrive - Johns Hopkins\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx\TT11','cl-maze1.1');
% A = 
% time = A.data(:,18);
clear
Nlx_directory = '~/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02/Neuralynx';
Nlx_directory = uigetdir(Nlx_directory,'Select Neuralynx Folder');

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

A = cellfun(@(x) importdata(x,',',13), absolue_paths, 'UniformOutput', false);

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
    start_time(index) = str2double(A{index}.textdata{12});
    end_time(index) = str2double(A{index}.textdata{13});
    %ts = start_time:1/30000:end_time;
    time = seconds((A{index}.data(:,18))*1e-6);
    all_ones = ones(length(time),1);
    
    cluster_maze = timetable(time,all_ones,'VariableNames',maze_tt_cl(index));
    disp index;
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