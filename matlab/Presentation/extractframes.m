
%% Write an Video File for a Lap
close all
clc; clear;

exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename, 'lap', 'pos', 'cluster');
% 
s = '';

for l=20:20
    timerange = lap(l).t_jump + [-2 2];             % Time range of recording (2 sec before to 2 sec after jump)
    videoframes = pos.frame(pos.t>=timerange(1) & pos.t<=timerange(end)); % frames in the timerange
    command = sprintf('rsync -tvrPhe "ssh -p 25" dome3tracking@dome3router:/media/dome3tracking/Videos/2020-11-22/frame-{%d..%d}.raw /home/shahin/Downloads/test-raw/%d/',videoframes(1),videoframes(end),l);
    system(command);
    fprintf('\n The files for lap %d has been copied!\n',l);
end