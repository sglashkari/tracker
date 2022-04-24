clc; clear; close all

exp_directory = 'D:\Analysis\2021-12-21';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'colors','xmax','hist','daq');


cluster_no = 29; %[17 18 23 28 31 39]; %[7 28 29 40]; %24 ; %[23 28]; %[17 18 23 28]; % [22 29]; [3 13 17 18 22 23 28 29 31]; [cluster.no];
colors(18)="#f58231";
colors(28)="#3cb44b";
close all

s_thresh = 10;
direction = 'right';
showCSC = false;
isZoom = true;
timerange = [-2 2]; % 2 sec before to 2 sec after
x_reference = 'moving';
t_reference = 'landing';
opacity = 0.5;

mydir = fullfile(exp_directory, 'Analysis',['cluster' num2str(cluster_no)]);
if exist(mydir, 'dir')
    warning('Directory removed')
    rmdir(mydir, 's');
end
mkdir(mydir);
mydoc = publish('plottime','outputDir',mydir,'showCode', false);
movefile(mydoc, fullfile(mydir, ['cluster' num2str(cluster_no), '.html']));