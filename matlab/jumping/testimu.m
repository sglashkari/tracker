close all;
TimeRange = [2122802110	2190641110];
nlx_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03/Neuralynx';
imu = readimu(nlx_directory,TimeRange);
dt = 2156.416735;
imu.t = imu.t - dt;
figure(1); hold on

ax1 = subplot(2,1,1); hold on
plot(imu.t, [imu.wx imu.wy imu.wz]*1e6);

ax2 = subplot(2,1,2); hold on
plot(imu.t, [imu.ax imu.ay imu.az]*1e6);

linkaxes([ax1 ax2],'x')

accelReadings = [imu.ax imu.ay imu.az];
gyroReadings = [imu.wx imu.wy imu.wz];

decim = 1;
Fs = 3000;
fuse = imufilter('SampleRate',Fs,'DecimationFactor',decim);

orientation = fuse(accelReadings,gyroReadings);

time = (0:decim:size(accelReadings,1)-1)/Fs;

figure(2)
plot(time,eulerd(orientation,'ZYX','frame'))
title('Orientation Estimate')
legend('Z-axis', 'Y-axis', 'X-axis')
xlabel('Time (s)')
ylabel('Rotation (degrees)')

%xlim([-2 2])
