
%% Write an Video File for a Lap
close all
clc; clear;

exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename, 'lap', 'pos', 'cluster', 'ppcm', 'colors');

l = 2;
c = 6;

ppm_directory = ['/home/shahin/Downloads/test-raw/' num2str(l)];
framerate = 300;

tic
%%
timerange = lap(l).t_jump + [-2 2];             % Time range of recording (2 sec before to 2 sec after jump)

videoframes = pos.frame(pos.t>=timerange(1) & pos.t<=timerange(end)); % frames in the timerange
videotimes_var = pos.t(pos.t>=timerange(1) & pos.t<=timerange(end)); % times of the frames in the timerange (variable framerate)
videotimes_cte = (videotimes_var(1):1/framerate:videotimes_var(end))'; % times of the frames in the timerange (constant framerate)

system('counter=000');
command=sprintf('for f in frame-{%d..%d}.raw; do ffmpeg -f rawvideo -r 1 -s 1600x580 -pix_fmt gray -i %s/"$f" %s/lap%d-frame-$((counter)).ppm; ((counter++)); done'...
    ,videoframes(1),videoframes(end),ppm_directory,ppm_directory,l);
system(command);   

%fprintf('The error for the constant frame rate assumption is %.2f milliseconds at most.\n',max(videotimes_var-videotimes_cte)*1000)

%% plot spikes on ppm files
idx = (cluster(c).t>=videotimes_cte(1)) & (cluster(c).t<=(videotimes_cte(end)));
nnz(idx)

for frameNo = 1:length(videoframes)
    idx = (cluster(c).t>=videotimes_cte(frameNo)-3/framerate) & (cluster(c).t<=(videotimes_cte(frameNo)+3/framerate)); % spikes at the timerange of the frame
    if nnz(idx)
        disp(frameNo)
        ppm_filename = [ppm_directory filesep 'lap' num2str(l) '-frame-' num2str(frameNo-1) '.ppm'];  % ppm starts from frame-0
        I = imread(ppm_filename);
        imshow(I); hold on
        h = plot(cluster(c).x(idx) * ppcm, cluster(c).y(idx) * ppcm, 'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c), 'MarkerSize', 20);
        pause(0.2);
        imwrite(getframe(gca).cdata, ppm_filename);
    end
end

%% saves a wav file for the spikes
Fs = 44100;                                     % Sampling Frequency
t = timerange(1):1/Fs:timerange(end);           % Time Vector
spiketimes = cluster(c).t(cluster(c).lap==l);   % Time of spike in the lap

% Create Tone
d = [spiketimes ones(size(spiketimes))];
y = pulstran(t, d,'rectpuls',1e-3);

figure(2)
plot(t, y)
xlabel('Time (s)')
xlim(timerange)
ylabel('Waveform')

filename = sprintf('%s/audio.wav',ppm_directory);
audiowrite(filename,y,Fs);

%% create video with sound
command = sprintf('ffmpeg -framerate %d -s 1600x580 -i %s/lap%d-frame-%%d.ppm -i %s/audio.wav -c:v libx264 -crf 18 -r %d %s/video.avi', ...
    framerate, ppm_directory, l, ppm_directory, framerate, ppm_directory);
system(command);

%% play video
toc
system(['xdg-open ' ppm_directory '/video.avi'])