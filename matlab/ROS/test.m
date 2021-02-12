%% This program measures the difference between software PTP and LED 
% Shahin 2021-01-20
% 
clc; close all
% clear;
rosbag_filename = '/home/shahin/Desktop/test_latency/2021-02-11/openArenaCamera_2021-02-11-20-01-11.bag';
event_filename = '/home/shahin/Desktop/test_latency/2021-02-11/Events.nev';

%% Camera Time
% rosbag('info',rosbag_filename)
% bag = rosbag(rosbag_filename);
% 
% %duration = 10; % seconds
% bag = select(bag,'Time',[bag.StartTime, bag.EndTime]);
% 
% msgs = readMessages(bag);
% t = bag.MessageList.Time;

disp('Data Extracted!')

tic
% S = nan(length(t),1);
% for frame_no=1:length(t)
%     image = readImage(msgs{frame_no});
%     S(frame_no) = nnz(image==255);
%     if mod(frame_no,100)==0
%         fprintf('%.1f %% of the video is extracted in %.f sec, the remaining time is.\n',100*frame_no/length(t),toc);
%     end
% end
% imshow(image);
figure(2); hold on
plot(t,S);

% time_start = input('What is the start time? ');
% S2 = S(t>time_start);
% t2 = t(t>time_start);
% 
% plot(t2,S2,'k');
% 
brightness_threshhold = 50;
% 
% t_on = t(S>=brightness_threshhold);
% t_off = t(S2<brightness_threshhold);
% 
% Ton1 = t_on(1);
% Toff = t_off(t_off>Ton1);
% Toff1 = Toff(1);
% 
% Ton2 = t_on(15);
% Toff = t_off(t_off>Ton2);
% Toff2 = Toff(1);
% 
% t_camera = [Ton1 Ton2 Toff1 Toff2]
% S_camera = interp1(t2,S2,t_camera);
% figure(2)
% hold on
% plot(t_camera,S_camera,'or');

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

%save('/Users/shahin/Desktop/test_latency/rosbag.mat','t','S','t_selected')

%% Neuralynx Time

FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, Event IDs, TTLs, Extras, Event Strings
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];
addpath('../../pkgs/releaseDec2015/binaries'); % Neuralynx packages for Linux/Mac packages

[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV_v3(event_filename,...
        FieldSelectionFlags, HeaderExtractionFlag, ExtractionMode, ExtractionModeVector );
    
% addpath('../jumping');
% [Timestamps,~,~,EventIDs,TTLs] = readevent(event_file);

T1 = Timestamps(TTLs==1 & EventIDs==11);
T2 = Timestamps(TTLs==2 & EventIDs==11);
T3 = Timestamps(TTLs==3 & EventIDs==11);
T4 = Timestamps(TTLs==4 & EventIDs==11);
T5 = Timestamps(TTLs==5 & EventIDs==11);
T6 = Timestamps(TTLs==6 & EventIDs==11);

idx_start = strcmp(EventStrings,'Sync Signals');
t_start = Timestamps(idx_start);
idx_on = strcmp(EventStrings,'TTL Input on AcqSystem1_Base board 0 port 0 value (0x0002).');
t_gpio_on = Timestamps(idx_on)' * 1e-6
idx_off = strcmp(EventStrings,'TTL Input on AcqSystem1_Base board 0 port 0 value (0x0000).');
t_gpio_off = Timestamps(idx_off & (Timestamps>t_start(1))')' * 1e-6

%t_nlx = [t_on(1) t_off(1) t_on(2) t_off(2)] 

%%
utl_correction = -18000; % 5 hr = 18000 sec (EST UTC−05:00)
t_diff_on = t_led_on - t_gpio_on(1:4) + utl_correction
t_diff_off = t_led_off - t_gpio_off(1:4) + utl_correction


%}