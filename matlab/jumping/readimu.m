function imu = readimu(nlx_directory,TimeRange)
%READIMU reads the imu values from the Neuralynx folder
%   Detailed explanation goes here
if nargin == 0
    nlx_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03/Neuralynx';
    nlx_directory = uigetdir(nlx_directory,'Select Neuralynx Directory');
end
if nlx_directory == 0
    imu = [];
    return;
end

if nargin == 2
    [imu.t,imu.wx, header1] = readcsc(fullfile(nlx_directory,'BASE_AVX.ncs'), TimeRange); % microseconds
    [~,imu.wy] = readcsc(fullfile(nlx_directory,'BASE_AVY.ncs'),TimeRange); % microseconds
    [~,imu.wz] = readcsc(fullfile(nlx_directory,'BASE_AVZ.ncs'), TimeRange); % microseconds
    [~,imu.ax] = readcsc(fullfile(nlx_directory,'BASE_LAX.ncs'), TimeRange); % microseconds
    [~,imu.ay] = readcsc(fullfile(nlx_directory,'BASE_LAY.ncs'), TimeRange); % microseconds
    [~,imu.az] = readcsc(fullfile(nlx_directory,'BASE_LAZ.ncs'), TimeRange); % microseconds
else
    [imu.t,imu.wx] = readcsc(fullfile(nlx_directory,'BASE_AVX.ncs')); % microseconds
    [~,imu.wy] = readcsc(fullfile(nlx_directory,'BASE_AVY.ncs')); % microseconds
    [~,imu.wz] = readcsc(fullfile(nlx_directory,'BASE_AVZ.ncs')); % microseconds
    [~,imu.ax] = readcsc(fullfile(nlx_directory,'BASE_LAX.ncs')); % microseconds
    [~,imu.ay] = readcsc(fullfile(nlx_directory,'BASE_LAY.ncs')); % microseconds
    [~,imu.az] = readcsc(fullfile(nlx_directory,'BASE_LAZ.ncs')); % microseconds
    
end
header1
% g = 981; % cm/s^2
% 
% imu.wx = imu.wx / 12e-4 * 250;  % deg/s
% imu.wy = imu.wy / 12e-4 * 250;  % deg/s
% imu.wz = imu.wz / 12e-4 * 250;  % deg/s
% imu.ax = imu.ax / 12e-4 * 2 * g;% cm/s^2
% imu.ay = imu.ay / 12e-4 * 2 * g;% cm/s^2
% imu.az = imu.az / 12e-4 * 2 * g;% cm/s^2

imu.w = vecnorm([imu.wx imu.wy imu.wz]')';
imu.a = vecnorm([imu.ax imu.ay imu.az]')';

if nargout == 0
    clc
    ax1 = subplot(2,1,1);
    plot(imu.t,imu.w);
    ax2 = subplot(2,1,2);
    plot(imu.t,imu.a);
    linkaxes([ax1 ax2],'x')
    clear imu;
end

end