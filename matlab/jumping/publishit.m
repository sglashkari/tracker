%   Author Shahin G Lashkari
%
clc; clear; close all;
answer = inputdlg({'Rat', 'Date'},'Enter the rat number and the date',[1 30],{'1068', '2022-12-20'});
rat_no = answer{1};
date_str = answer{2};
%% Selecting the appropriate files
[datafile,exp_directory] = uigetfile(['E:\Rat' rat_no '\Analysis\' date_str '\processed_data.mat'],'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename);
xmax = exp.xmax;
ppcm = exp.ppcm;

% https://sashamaps.net/docs/resources/20-colors/
cluster_no = [4 41 25 17]; %[17 18]; % 23 28]; % [28 35]; % [17 18 23 28]; %[18 28]; % [15 22 28 35]; %29; %[17 18 23 28]; %[17 18 23 28 31 39]; %[7 28 29 40]; %24 ; %[23 28]; %[17 18 23 28]; % [22 29]; [3 13 17 18 22 23 28 29 31]; [cluster.no];
cluster_no = [3 15 18 32 31];
cluster_no = [15 32];
% cluster_no = [35 28];
% colors(3)="#42d4f4";  % cyan 
% colors(15)="#42d4f4";  % cyan 
% colors(18)="#42d4f4";  % cyan 
% colors(32)="#e6194B"; % red
% colors(31)="#e6194B"; % red

% colors(17)="cyan";
% colors(18)="#f58231";
% %colors(28)="#0096FF";
% colors(22) = "yellow";
% colors(31)="#3cb44b";
% colors(28)="#911eb4"; % purple
% 
% % colors(17)="#e6194B"; % red
% % colors(18)="#3cb44b"; % green
% % colors(23)="#f58231"; % orange
% % colors(28)="#42d4f4";  % cyan 

close all

s_thresh = 10;
direction = 'left';
showCSC = true;
CSCLaps = []; % [] for all
if isempty(CSCLaps)
    CSCLaps = [lap([lap.dir]==direction).no]; % only the desired direction
end
if ~showCSC
    CSCLaps = [];
end
isZoom = true;
timerange = [-4 4]; % [-2 2] 2 sec before to 2 sec after
x_reference = 'stationary';
t_reference = 'take-off';
opacity = 0.25;

if showCSC
    mydir = fullfile(exp_directory, 'Analysis',['cluster' num2str(cluster_no) '_CSC']);
else
    mydir = fullfile(exp_directory, 'Analysis',['cluster' num2str(cluster_no)]);
end
if exist(mydir, 'dir')
    warning('Directory removed')
    rmdir(mydir, 's');
end
mkdir(mydir);
mydoc = publish('plottime','outputDir',mydir,'showCode', false);
movefile(mydoc, fullfile(mydir, ['cluster' num2str(cluster_no), '.html']));