%% This program measures the difference between software PTP and LED 
% Shahin 2021-01-20 modified 2021-02-11
% 
clc; close all; format longG
% clear;
rosbag_filename = '/home/shahin/Desktop/test_latency/2021-02-11/openArenaCamera_2021-02-11-20-01-11.bag';
event_filename = '/home/shahin/Desktop/test_latency/2021-02-11/Events.nev';

brightness_threshhold = 50;
utc_correction = -18000; % 5 hr = 18000 sec (EST UTC−05:00)

%% Camera Time (LED)
rosbag('info',rosbag_filename)
bag = rosbag(rosbag_filename);

%duration = 10; % seconds
bag = select(bag,'Time',[bag.StartTime, bag.EndTime]);

msgs = readMessages(bag);
t = bag.MessageList.Time;

disp('Data Extracted!')

tic
S = nan(length(t),1);
for frame_no=1:length(t)
    image = readImage(msgs{frame_no});
    S(frame_no) = nnz(image==255);
    if mod(frame_no,100)==0
        fprintf('%.1f %% of the video is extracted in %.f sec, the remaining time is %.f sec.\n',...
           100*frame_no/length(t),toc, toc * (1-frame_no/length(t))/(frame_no/length(t)));
    end
end

figure(1); hold on
plot(t,S);

idx = diff(S>=brightness_threshhold);
idx_on = [false; (idx == 1)];
idx_off = [false; (idx == -1)];
t_on = t(idx_on);
t_off = t(idx_off);

idx_select = (t_off - t_on) > 0.25 & (t_off - t_on) < 0.35; % duration between on and off is about 300 ms
t_led_on = t_on(idx_select)
t_led_off = t_off(idx_select)

S_led_on = interp1(t,S,t_led_on);
S_led_off = interp1(t,S,t_led_off);
plot(t_led_on,S_led_on,'or')
plot(t_led_off,S_led_off,'oc')

figure(2)
imshow(readImage(msgs{t == t_led_on(1)}));

%% Neuralynx Time (GPIO)

FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, Event IDs, TTLs, Extras, Event Strings
HeaderExtractionFlag = 0;
ExtractionMode = 1;
ExtractionModeVector = [];
addpath('../../pkgs/releaseDec2015/binaries'); % Neuralynx packages for Linux/Mac packages

[Timestamps, EventIDs, TTLs, Extras, EventStrings] = Nlx2MatEV_v3(event_filename,...
        FieldSelectionFlags, HeaderExtractionFlag, ExtractionMode, ExtractionModeVector );
    
idx_start = strcmp(EventStrings,'Sync Signals');
t_start = Timestamps(idx_start);
idx_on = strcmp(EventStrings,'TTL Input on AcqSystem1_Base board 0 port 0 value (0x0002).');
t_gpio_on = Timestamps(idx_on)' * 1e-6
idx_off = strcmp(EventStrings,'TTL Input on AcqSystem1_Base board 0 port 0 value (0x0000).');
t_gpio_off = Timestamps(idx_off & (Timestamps>t_start(1))')' * 1e-6


%% Neuralynx Time (software PTP)

T1 = str2double(string(EventStrings(TTLs==1 & EventIDs==11)));
T2 = Timestamps(TTLs==2 & EventIDs==11)';
T3 = Timestamps(TTLs==3 & EventIDs==11)';
T4 =str2double(string(EventStrings(TTLs==4 & EventIDs==11)));
T1 = T1(T1~=0); % discard the incorrect times

% If there is any incompatibility between the legth of timestamps because 
% of the interruptions during the data collection, it will be fixed here:

l = [length(T1),length(T2),length(T3),length(T4)];
T1(min(l)+1:end)=[];
T2(min(l)+1:end)=[];
T3(min(l)+1:end)=[];
T4(min(l)+1:end)=[];

figure(3)
T2134 = (T2 - T1 + T3 - T4)/2;
plot(T1,T2134,'*')

[offset, S_T2134] = polyfit(T1,T2134, 0);
[fittedY,std_dev] = polyval(offset, T1, S_T2134);
offset_mean = mean(fittedY);
std_dev = mean(std_dev);
SEM = std_dev/sqrt(length(T1)); 

hold on;
plot(T1, fittedY, 'g-', 'LineWidth', 5);
plot(T1,fittedY+std_dev,'m--',T1,fittedY-std_dev,'m--', 'LineWidth', 3)

fprintf('The software PTP offset is %.f ± %.2f milliseconds.\n',offset_mean * 1e-3 - utc_correction * 1e3,SEM*1e-3);

%% Offset between the GPIO and LED (most accurate measurement)

t_diff_on = t_led_on - t_gpio_on(1:4) + utc_correction
t_diff_off = t_led_off - t_gpio_off(1:4) + utc_correction

fprintf('The LED offset is %.f ± %.f milliseconds.\n',1000*mean([t_diff_on;t_diff_off]),1000*std([t_diff_on;t_diff_off]))
%}