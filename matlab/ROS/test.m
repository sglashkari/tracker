clear; 
clc; close all

%ros_bag_file_name = '/home/shahin/Desktop/domeExperimentCamera_2019-09-07-13-40-02.bag';
ros_bag_file_name = '/home/shahin/Desktop/domeExperimentCamera_2019-09-07-13-40-02.bag';

rosbag('info',ros_bag_file_name)
bag = rosbag(ros_bag_file_name);

duration = 10; % seconds
bag = select(bag,'Time',[bag.StartTime, bag.StartTime + duration]);

msgs = readMessages(bag);
t = bag.MessageList.Time;

disp('Data Extracted!')

%%
tic
S = nan(length(t),1);
for j=1:length(t)
    image = readImage(msgs{j});
    %S(j) = mean(image,'all');
    S(j) = nnz(image==255);
    %imshow(image);
    %pause(t(j+1)-t(j));
    if mod(j,100)==0
        disp(j);
        toc
    end
end
imshow(image);
figure
plot(t,S);
toc