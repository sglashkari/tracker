clc; close all;
vid_filename = 'D:\200629_test_45fps_noRecording.avi';
videoReader  = VideoReader(vid_filename);
% videoReader.NumFrame
% 1683
c=zeros(1683,4);
tic
for frame_no = 1:1683%videoReader.NumFrames
% Read image (i.e. a frame) from a video
I = read(videoReader,frame_no);
%hold off
%ax = imshow(I);
%hold on

%% Projected image
% ref: https://www.mathworks.com/matlabcentral/answers/495720-how-to-fit-an-ellipse-to-an-image-in-matlab#answer_405863
xmin = 0;
ymin = 450;
I2 = imcrop(I,[xmin ymin 680 400]);

% Binarize
Igray = rgb2gray(I2);
BW = imbinarize(Igray);

% Extract the maximum area
BW = imclearborder(BW,8);
BW = bwareafilt(BW,1);

% Calculate centroid, orientation and major/minor axis length of the ellipse
s = regionprops('table', BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
x_proj = s.Centroid(1)+xmin;
y_proj = s.Centroid(2)+ymin;

%% crown
[centers, radii, metric] = imfindcircles(I,[15 30]);
x_crown = centers(1);
y_crown = centers(2);

%plot(x_proj, y_proj, 'r+', x_crown, y_crown, 'b+', 'MarkerSize', 15, 'LineWidth', 2);
if mod(frame_no,100)==0
    disp(frame_no);
    toc
end
%pause(0.03);

c(frame_no,:) = [x_proj, y_proj, x_crown, y_crown];
end

k=1:length(c);
t=(k-1)/videoReader.FrameRate;
fr = videoReader.FrameRate;
save('200629_test_45fps_noRecording.mat','c','t','fr')