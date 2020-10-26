function track_led1(vid_filename)
%% TRACK_LED1 tracks the LED on the FreeLynx
%
% Input:
%   video file of jumping experiment (Rat 913)
%
% Output:
%   mat file
%
% Author: Shahin G Lashkari
% Date: October 26, 2020
%
clc; close all;
if nargin==0
    vid_filename = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\913_Jumping_Recording\2020-10-22_12-38-17.mp4';
    vid_filename = 'D:\OneDrive - Johns Hopkins\HU\913_Jumping_Recording\2020-10-23-Day0\2020-10-23-133050.mp4'; % day 0
    vid_filename = 'D:\OneDrive - Johns Hopkins\JHU\913_Jumping_Recording\2020-10-25_Rat913-01\Video\2020-10-25_15-59-09.mp4'; %day 1 
end
f = waitbar(0,'Please wait...');
videoReader  = VideoReader(vid_filename)
%videoWriter = VideoWriter(insertBefore(vid_filename,".","_tracked"));
num_frames = videoReader.NumFrame % 1683
frame_rate = videoReader.FrameRate;
%videoWriter.FrameRate = frame_rate;
%videoWriter.Quality = 100;
position=zeros(num_frames,3);
%open(videoWriter);
figure(1);
tic
for frame_no = 1:num_frames
    % Read image (i.e. a frame) from a video
    I = read(videoReader,frame_no);
    
    hold off
    imshow(I);
    hold on
    
    %% Projected image
    % ref: https://www.mathworks.com/matlabcentral/answers/495720-how-to-fit-an-ellipse-to-an-image-in-matlab#answer_405863
    xmin = 0;
    ymin = 0;
    I2 = imcrop(I,[xmin ymin 1920 1080]);
    
    % Binarize
    Igray = rgb2gray(I2);
    threshold = 250; % custom threshold value
    BW = Igray > threshold;
    
    % Extract the maximum area
    BW = imclearborder(BW);
    BW = bwareafilt(BW,1);
    
    % Calculate centroid, orientation and major/minor axis length of the ellipse
    if nnz(BW) ~= 0
        s = regionprops('table', BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
        x_proj = s.Centroid(1)+xmin;
        y_proj = s.Centroid(2)+ymin;
        position(frame_no,:) = [frame_no, x_proj, y_proj];
    else
        x_proj = [];
        y_proj = [];
    end
    
    plot(x_proj, y_proj, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    
    if mod(frame_no,100)==0
        elapsed_time = toc;
        remaining_time = num2str(round(elapsed_time*(num_frames/frame_no-1)/60));
        waitbar(frame_no/num_frames,f,['Remaining time is ' remaining_time ' minutes.']);
        disp(frame_no);
    end
    
    frame = getframe(gcf);
    %writeVideo(videoWriter,frame);
    pause(1/frame_rate);
    
end

%close(videoWriter)

k=1:length(position);
t=(k-1)/frame_rate;

waitbar(1,f,'Finishing');
save(strrep(vid_filename,'mp4','mat'),'position','t','frame_rate');
close(f)
end