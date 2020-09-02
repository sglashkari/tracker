% This is for testing occupany sgl 2020-08-30
clc; close all;

Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
Filename = fullfile(Nlx_directory,'pos.p');

TimeEV = readevent;
StartTime = TimeEV(1:2:end);
EndTime = TimeEV(2:2:end);

%% Pos.p
maze = 4;
[TimePos,x,y] = readposp(Filename, StartTime(maze), EndTime(maze));
y = 480 - y;
AngularPosition = rad2deg(atan2(y-240, x-320));
AngularPosition = wrapTo360(AngularPosition);

%% interpolation
f = 3000; % interpolation 3 kHz
TimePosInterp = TimePos(1):1/f: TimePos(end);

xInterp = interp1(TimePos,x,TimePosInterp);
yInterp = interp1(TimePos,y,TimePosInterp);
AngularPositionInterp = rad2deg(atan2(yInterp-240, xInterp-320));
AngularPositionInterp = wrapTo360(AngularPositionInterp);

%% histogram
edges = linspace(0, 360, 360);
HistOcc = histcounts(AngularPosition, edges);
HistOccInterp = histcounts(AngularPositionInterp, edges);

%% plotting
f1 = figure(1);
plot(TimePos, AngularPosition,'o');
hold on
plot(TimePosInterp, AngularPositionInterp,'.');
f2 = figure(2);
bar(HistOccInterp);

%% velocity
% 5 ft = 1524 mm = 480 pixels
% each pixel = 3.175 mm

vx = diff(x)./diff(TimePos); % pixels/sec
vy = diff(y)./diff(TimePos); % pixels/sec

velocity = sqrt(vx.^2+vy.^2); % pixels/sec
velocity = velocity * 0.3175; % cm/sec

f3 = figure(3);

TimePos = TimePos(2:end);
TimePosVelFil = TimePos(velocity<5);
velocityVelFil = velocity(velocity<5);
plot(TimePosVelFil, velocityVelFil,'o')


