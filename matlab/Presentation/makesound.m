%% Write an Audio File for Spike Times
% Create a WAVE file from the spike times and read the file
% back into MATLAB.   
close all
clc; clear;

exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename, 'lap', 'cluster');

l = 2;
c = 6;

%% Write a WAVE file in the analysis folder. 
timerange = lap(l).t_jump + [-2 2];             % Time range of recording (2 sec before to after jump)
Fs = 44100;                                     % Sampling Frequency
t = timerange(1):1/Fs:timerange(end);           % One Second Time Vector
spiketimes = cluster(c).t(cluster(c).lap==l);   % Time of spike in the lap

% Create Tone
d = [spiketimes ones(size(spiketimes))];
% check https://www.mathworks.com/help/signal/gs/the-pulstran-function.html
y = pulstran(t, d,'rectpuls',1e-3);

figure(1)
plot(t, y)
xlabel('Time (s)')
xlim(timerange)
ylabel('Waveform')

filename = fullfile(exp_directory, 'Analysis',['cluster' num2str(c) '_lap' num2str(l) '.wav']);
audiowrite(filename,y,Fs);

%% check if the recording has been successful
clear y Fs
figure(1)
[y,Fs] = audioread(filename);   % Read the data back
sound(y,Fs);                    % Listen to the audio