%% Troubleshooting
% catch circles for test if the algorithm works
clc
close all
[file, path] = uigetfile('D:\OneDrive - Johns Hopkins\JHU\913_Jumping_Recording\2020-10-25_Rat913-01\Videos\2020-10-25_15-59-09.mp4');
vid_filename = fullfile(path,file);
videoReader  = VideoReader(vid_filename)

% Read image (i.e. a frame) from a video
I = read(videoReader,8099); %  3478 2045
figure(1);
imshow(I);

%% Projected image
% ref: https://www.mathworks.com/matlabcentral/answers/495720-how-to-fit-an-ellipse-to-an-image-in-matlab#answer_405863
xmin = 350;
ymin = 300;
I2 = imcrop(I,[xmin ymin 850 330]);

xmin = 0;
ymin = 0;
I2 = imcrop(I,[xmin ymin 1920 1080]);

% Binarize
Igray = rgb2gray(I2);
threshold = 250; % custom threshold value
BW = Igray > threshold;
figure(2);
imshow(BW);

% Extract the maximum area
BW = imclearborder(BW);
BW = bwareafilt(BW,1);

% Calculate centroid, orientation and major/minor axis length of the ellipse
if nnz(BW) ~= 0
    s = regionprops('table', BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
    x_proj = s.Centroid(1)+xmin;
    y_proj = s.Centroid(2)+ymin;
else
    x_proj = [];
    y_proj = [];
end

figure(1)
hold on
plot(x_proj, y_proj, 'r+', 'MarkerSize', 20, 'LineWidth', 2);