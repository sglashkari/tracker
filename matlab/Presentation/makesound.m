%% Write an Audio File for Spike Times
% Create a WAVE file from the spike times and read the file
% back into MATLAB.   
clc; clear; close all;

exp_directory = 'D:\Analysis\2021-12-10';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename, 'lap', 'cluster');

l = 17;
c = 43;

fps = [20 25 30];


for f = fps
%% Write a WAVE file in the analysis folder. 
speed = 100./f;
%timerange = lap(l).t_jump + [-2 2];             % Time range of recording (2 sec before to after jump)
timerange = [255.0133  258.3058]*speed;
Fs = 44100;                                     % Sampling Frequency
t = timerange(1):1/Fs:timerange(end);           % One Second Time Vector
spiketimes = cluster(c).t(cluster(c).lap==l)*speed;   % Time of spike in the lap

% Create Tone
d = [spiketimes ones(size(spiketimes))];
% check https://www.mathworks.com/help/signal/gs/the-pulstran-function.html
y = pulstran(t, d,'rectpuls',1e-3);

figure(1)
plot(t, y)
xlabel('Time (s)')
xlim(timerange)
ylabel('Waveform')

filename = fullfile(exp_directory, 'Analysis',['cluster' num2str(c) '_lap' num2str(l) '_speed' num2str(f) 'fps.wav']);
audiowrite(filename,y,Fs);

%% check if the recording has been successful
clear y Fs
figure(1)
[y,Fs] = audioread(filename);   % Read the data back
sound(y,Fs); % Listen to the audio

end