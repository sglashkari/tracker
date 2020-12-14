close all; clear;
exp_directory = '/home/shahin/Desktop/sync-2020-12-11';
exp_directory = '/home/shahin/onedrive/JHU/913_Jumping_Recording/2020-11-22_Rat913-03';
new_method = false;

if new_method == true
    N = 4;
else
    N = 2;
end
N =1;

%% new sync method:
% SGL 2020-12-12
addpath('../tracking');
[t, x, ~, p, frame, ~] = readtrackingdata(exp_directory); 
p = cast(p(:,3),'int8');
plot(t, 2*p)
hold on
plot(t(2:end),logical(diff(1-p)))
t_gpio = t(2:end);
t_gpio = t_gpio(diff(p)~=0)
t_gpio = t_gpio(1:2)


figure(2)
plot(t, x)

port = [0 1]';


event_file_name = fullfile(exp_directory,'Neuralynx', 'Events.nev');
addpath('../jumping');
[Time,Data,Header,EventIDs,TTLs] = readevent(event_file_name);

Time(EventIDs == 21 | EventIDs == 11)
Data(EventIDs == 21 | EventIDs == 11)
TTLs(EventIDs == 21 | EventIDs == 11)
EventIDs(EventIDs == 21 | EventIDs == 11)

t_nlx = Time(EventIDs == 21);
t_nlx = t_nlx(end-1:end)

offset_new = t_nlx - t_gpio

disp(['The new offset is ' num2str(mean(offset_new)) ' sec +/-' num2str(std(offset_new)*1e3) ' msec.']);

figure(3)
plot(t + mean(offset_new), p, t_nlx, port, '*')

xlim([min(t_nlx) - 5  max(t_nlx) + 5])

%% old sync method
% 


T5 = Time(TTLs == 5);
T5 = T5(end)
T6 = Time(TTLs == 6);
T6 = T6(end)

% Day1 530 542
% Day2 9038 9050
% Day3 3976 3988
% Day4 1165 1177 = 231 243
% 913: 
% Day2 133 153

frames_tracking_led = [26765 27486]';
light_on_cam = t(frames_tracking_led(1)+1); % + 1 frame 0 is index 1 in MATLAB
light_off_cam = t(frames_tracking_led(2)+1);

offset_1 = T5 - light_on_cam;
offset_2 = T6 - light_off_cam;

DT1 = T6-T5;
DT2 = light_off_cam-light_on_cam;
offset_old = (offset_1 + offset_2)/2;

disp(['The old offset was ' num2str(offset_old) ' sec +/-' num2str(std([offset_1 offset_2])*1e3) ' msec.']);

Delta_offset = offset_new - offset_old

disp(['The difference between the two offset is ' num2str(mean(Delta_offset)) ' sec +/-' num2str(std(Delta_offset)*1e3) ' msec.']);
