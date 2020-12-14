close all; clear;
d = 3;
switch d
    case 3 % Day3
        syncit_Day3;
        %exp_directory = '/home/shahin/onedrive/JHU/913_Jumping_Recording/2020-11-22_Rat913-03';
        error('Break');
    case 11 % sync-2020-12-11
        exp_directory = '/home/shahin/Desktop/sync-2020-12-11';
        port = [1 0 1 0]';
        frame_led = [22340 22415 23015 23090]';
        frames_tracking_led = [22290 22965]';
    case 12 % sync-2020-12-12
        exp_directory = '/home/shahin/Desktop/sync-2020-12-12';
        % off on   off  on   off   on
        port = [0 1 0 1 0 1]';
        frame_led = [743 793 1144 1219 1820 1895]';
        frames_tracking_led = [1093 1769]';
end
disp(exp_directory)
%% new sync method:
% SGL 2020-12-12
addpath('../tracking');
[t, x, ~, p, frame, ~] = readtrackingdata(exp_directory); 
p = cast(p(:,3),'int8');
ax1 = subplot(4,1,1);
plot(t, p)
xlabel('Camera Time (Sec)')
ylabel('Port')
ylim([-0.5 1.5])

ax2 = subplot(4,1,2);

plot(t(2:end),logical(diff(1-p)),'+','MarkerSize',12,'LineWidth',2)
t_gpio = t(2:end);
t_gpio = t_gpio(diff(p)~=0);


xlabel('Camera Time (Sec)')
ylabel('Port Change')
ylim([0.5 1.5])
hold on

ax3 = subplot(4,1,3);
plot(t, x)
xlabel('Camera Time (Sec)')
ylabel('Posiotion (cm)')

ax4 = subplot(4,1,4);
plot(frame, x)
xlabel('Camera Frame No.')
ylabel('Posiotion (cm)')

%error('Break');
t_led = t(frame_led + 1);

ax2 = subplot(4,1,2);
hold on
plot(t_led,ones(size(t_led)),'o','MarkerSize',12,'LineWidth',2)
legend('port','LED')

ax1 = subplot(4,1,1);
hold on
plot(t_led,port,'or','MarkerSize',8,'LineWidth',2)
legend('port','LED')

linkaxes([ax1 ax2 ax3], 'x')

t_led = t_led(end-3:end);
t_gpio = t_gpio(end-3:end);
port = port(end-3:end);

if ~isequal(t_led,t_gpio)
    warning('Check the frame numbers for LED.');
else
    disp('LED is in accord with GPIO.');
end

event_file_name = fullfile(exp_directory,'Neuralynx', 'Events.nev');
addpath('../jumping');
[Time,Data,Header,EventIDs,TTLs] = readevent(event_file_name);

t_nlx = Time(EventIDs == 21);

t_nlx = t_nlx(end-3:end);


offset_new = t_nlx - t_gpio;

disp(['The offset (TTL) is ' num2str(mean(offset_new)) ' sec +/-' num2str(std(offset_new)*1e3) ' msec.']);

figure(2)

plot(t(2:end) + mean(offset_new),logical(diff(1-p)),'+','MarkerSize',16,'LineWidth',2)
hold on;
plot(t_nlx, ones(size(t_nlx)),'o','MarkerSize',16,'LineWidth',2)
xlim([min(t_nlx)-0.5  max(t_nlx)+0.5])
ylim([0.5 1.5])
xlabel('Neuralynx Time (Sec)')
legend('Camera TTL Received','Neuralynx TTL Sent')

%% old sync method

T5 = Time(TTLs == 5);
T5 = T5(end);
T6 = Time(TTLs == 6);
T6 = T6(end);

light_on_cam = t(frames_tracking_led(1)+1); % + 1 frame 0 is index 1 in MATLAB
light_off_cam = t(frames_tracking_led(2)+1);

offset_1 = T5 - light_on_cam;
offset_2 = T6 - light_off_cam;

DT1 = T6-T5;
DT2 = light_off_cam-light_on_cam;
offset_old = (offset_1 + offset_2)/2;

disp(['The old offset (tracking LED) was ' num2str(mean([offset_1 offset_2])) ' sec +/-' num2str(std([offset_1 offset_2])*1e3) ' msec.']);
Delta_offset = mean(offset_new) - offset_old;

disp(['The difference between the two offset is ' num2str(Delta_offset) ' sec.']);