%% Extracts frames from a rosbag
% Shahin 2020-12-01
clear; clc; close all

ros_bag_file_name = '/home/shahin/Desktop/domeExperimentCamera_2019-09-07-13-40-02.bag';

rosbag('info',ros_bag_file_name)
bag = rosbag(ros_bag_file_name);

duration = 10; % seconds
bag = select(bag,'Time',[bag.StartTime, bag.StartTime + duration]);

msgs = readMessages(bag);
t = bag.MessageList.Time;

for j=1:length(t)
    frame = readImage(msgs{j});
    imshow(frame);
end