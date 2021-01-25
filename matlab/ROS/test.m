%clear; 
% Shahin 2021-01-20
% 
clc; close all

ros_bag_file_name = '/Users/shahin/Desktop/test_latency/openArenaCamera_2021-01-20-20-29-16.bag';
event_file = '/Users/shahin/Desktop/test_latency/Events.nev';

%% Camera Time
rosbag('info',ros_bag_file_name)
bag = rosbag(ros_bag_file_name);

%duration = 10; % seconds
bag = select(bag,'Time',[bag.StartTime, bag.EndTime]);

msgs = readMessages(bag);
t = bag.MessageList.Time;

disp('Data Extracted!')

tic
S = nan(length(t),1);
for j=1:length(t)
    image = readImage(msgs{j});
    S(j) = nnz(image==255);
    if mod(j,100)==0
        disp(j/length(t));
        toc
    end
end
imshow(image);
figure(2)
plot(t,S);

threshold = 50;
t_on = t(S>=threshold);
t_off = t(S<threshold);

Ton1 = t_on(1);
Toff = t_off(t_off>Ton1);
Toff1 = Toff(1);

Ton2 = t_on(15);
Toff = t_off(t_off>Ton2);
Toff2 = Toff(1);

t_camera = [Ton1 Ton2 Toff1 Toff2]
S_camera = interp1(t,S,t_camera);
figure(2)
hold on
plot(t_camera,S_camera,'or');

save('/Users/shahin/Desktop/test_latency/rosbag.mat','t','S','t_selected')

%% Neuralynx Time


FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, Event IDs, TTLs, Extras, Event Strings
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];
addpath('../../pkgs/releaseDec2015/binaries'); % Neuralynx packages for Linux/Mac packages

[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV_v3( event_file,...
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
t_on = Timestamps(idx_on & (Timestamps>t_start)');
idx_off = strcmp(EventStrings,'TTL Input on AcqSystem1_Base board 0 port 0 value (0x0000).');
t_off = Timestamps(idx_off & (Timestamps>t_start)');

t_nlx = [t_on(1) t_off(1) t_on(2) t_off(2)] * 1e-6

t_camera - t_nlx