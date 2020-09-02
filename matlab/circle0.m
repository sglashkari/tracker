%% Troubleshooting
% catch circles for test if the algorithm works
clc
close all
vid_filename = 'D:\domeLatencyMeasurement\200731_latencyTest\200731_latencyTest_60fps.avi';
videoReader  = VideoReader(vid_filename)

% Read image (i.e. a frame) from a video
I = read(videoReader,3823); %  3478 2045
figure(1);
hold off
ax = imshow(I);
hold on

%% Projected image
% ref: https://www.mathworks.com/matlabcentral/answers/495720-how-to-fit-an-ellipse-to-an-image-in-matlab#answer_405863
xmin = 0;
ymin = 350;
I2 = imcrop(I,[xmin ymin 500 150]);

% Binarize
Igray = rgb2gray(I2);
BW = imbinarize(Igray);
BW(:,1) = 0;

% Extract the maximum area
BW = imclearborder(BW,8);
BW = bwareafilt(BW,1);

figure(2);
imshow(BW);

% Calculate centroid, orientation and major/minor axis length of the ellipse
s = regionprops('table', BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
x_proj = s.Centroid(1)+xmin;
y_proj = s.Centroid(2)+ymin;

%% crown
[centers, radii, metric] = imfindcircles(I,[8 20]);

figure(1)
viscircles(centers, radii,'EdgeColor','b');
hold on
% [a,b]=max(centers(:,1))
% if b~=1
%     warning('check it out!')
%     if centers(2,1)+100<centers(1,1)
%         disp(b)
%         %centers
%         centers = centers(b,:);
%     else
%         centers = centers(1,:);
%     end
% else
%     centers = centers(1,:);
% end
if size(centers,1)>1
    fprintf('check frame\n')
    if abs(centers(2,1)-centers(1,1))<100
        if centers(2,1)>centers(1,1)
            centers = centers(2,:);
        else
            centers = centers(1,:);
        end
    else
        if centers(2,1)>centers(1,1)
            centers = centers(1,:);
        else
            centers = centers(2,:);
        end
    end
end


if isempty(centers)
    centers = [nan nan];
end
x_crown = centers(1);
y_crown = centers(2);

plot(x_proj, y_proj, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
plot(x_crown, y_crown, 'b+', 'MarkerSize', 12, 'LineWidth', 2);