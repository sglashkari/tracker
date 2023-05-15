%% Write an Video File for a Lap
% This is not the final version
% in this version the video is created in the matlab
% for final version check makevideowithaudio

close all
clc; clear;

exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename, 'lap', 'pos', 'cluster', 'ppcm', 'colors');

vid_filename = '/home/shahin/Downloads/test-raw/1/compressed.avi';

l = 1;
c = 6;
timerange = lap(l).t_jump + [-2 2];             % Time range of recording (2 sec before to 2 sec after jump)

videoframes = pos.frame(pos.t>=timerange(1) & pos.t<=timerange(end)); % frames in the timerange
videotimes_var = pos.t(pos.t>=timerange(1) & pos.t<=timerange(end)); % times of the frames in the timerange (variable framerate)


%% How much error constant frame rate has

figure(1); hold on
plot(videotimes_var)
videotimes_cte = linspace(videotimes_var(1),videotimes_var(end),length(videotimes_var))'; % times of the frames in the timerange (constant framerate)
videotimes_cte = videotimes_cte + mean(videotimes_var-videotimes_cte);
%plot(videotimes_cte)

framerate = 1/mean(diff(videotimes_var))



videotimes_cte = (videotimes_var(1):1/framerate:videotimes_var(end))'; % times of the frames in the timerange (constant framerate)
plot(videotimes_cte)

fprintf('The error for the constant frame rate assumption is %.2f milliseconds at most.\n',max(videotimes_var-videotimes_cte)*1000)

%% Write a MP4 file in the analysis folder.

videoReader  = VideoReader(vid_filename);
videoWriter = VideoWriter(insertBefore(vid_filename,".","_marked"));
numframes = videoReader.NumFrame;
framerate = videoReader.FrameRate;

if (numframes == length(videoframes))
    disp('sizes match')
else
    disp('sizes don''t match')
end

videoWriter.Quality = 100;
videoWriter.FrameRate = framerate;
open(videoWriter);

figure(2);
for frame_no=1:length(videoframes)
    I = read(videoReader,frame_no);
    
    imshow(I); hold on
    idx = (cluster(c).t>=videotimes_cte(frame_no) & cluster(c).t<=videotimes_cte(frame_no)); % spike at the timerange of the frame
    if nnz(idx)
        plot(cluster(c).x(idx) * ppcm, cluster(c).y(idx) * ppcm, 'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        pause(0.1);
    end
    frame = getframe(gcf);
    writeVideo(videoWriter,frame);
    pause(1/framerate);
    %delete(h);
end
close(videoWriter)
% vid_filename = fullfile(exp_directory, 'Analysis',['cluster' num2str(c) '_lap' num2str(l) '.avi']);
% videoWriter = VideoWriter(vid_filename);