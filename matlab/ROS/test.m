%clear; 
clc; close all

%ros_bag_file_name = '/home/shahin/Desktop/domeExperimentCamera_2019-09-07-13-40-02.bag';
%ros_bag_file_name = '/home/shahin/Desktop/domeExperimentCamera_2019-09-07-13-40-02.bag';
ros_bag_file_name = '/home/shahin/Desktop/test_latency/openArenaCamera_2021-01-20-20-29-16.bag';

rosbag('info',ros_bag_file_name)
bag = rosbag(ros_bag_file_name);

%duration = 10; % seconds
bag = select(bag,'Time',[bag.StartTime, bag.EndTime]);

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
        disp(j/length(t));
        toc
    end
end
imshow(image);
figure(2)
plot(t,S);
toc

%%
threshold = 50;
t_on = t(S>=threshold);
t_off = t(S<threshold);

Ton1 = t_on(1);
Toff = t_off(t_off>Ton1);
Toff1 = Toff(1);

Ton2 = t_on(15);
Toff = t_off(t_off>Ton2);
Toff2 = Toff(1);

t_selected = [Ton1 Ton2 Toff1 Toff2];
S_selected = interp1(t,S,t_selected);
figure(2)
hold on
plot(t_selected,S_selected,'or');
t_selected
save('/home/shahin/Desktop/test_latency/rosbag.mat','t','S','t_selected')