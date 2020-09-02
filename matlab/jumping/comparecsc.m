%% Plot all CSCchannels in search for best TT
% Shahin G Lashkari August 10th, 2020
clc;
clear; close all;

Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
listing = dir(fullfile(Nlx_directory,'**','TT*.ntt'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing);
absolue_paths = folders+'\'+names;

%% maze
maze = 4;

TimeEV = readevent;
StartTime = TimeEV(1:2:end);
EndTime = TimeEV(2:2:end);

%% Pos.p
[TimePos,x,y] = readposp;
y = 480 - y;
AngularPosition = rad2deg(atan2(y-240, x-320));
AngularPosition = wrapTo360(AngularPosition);

%% lap detector
lap = 7;

TimePosF = TimePos(~isnan(AngularPosition));
AngularPositionF = AngularPosition(~isnan(AngularPosition));
AngularPositionF = AngularPositionF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));
TimePosF =  TimePosF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));

LapTime = TimePosF(diff(AngularPositionF) < -350);
TimeRange = [LapTime(lap) LapTime(lap+1)]; % seconds


%% CSC
TimeTT=[];
f = figure(2);
for i=1:N
    TTT = readtt(convertStringsToChars(absolue_paths(i)), TimeRange * 1e6); % microseconds
    FilenameCSC = fullfile(Nlx_directory,['CSC' num2str(i) '.ncs']);
    [TimeCSC,Data,Header] = readcsc(FilenameCSC, TimeRange * 1e6); % microseconds
    ax = plot(TimeCSC,-1200*i+Data*1e6);
    Header(21)
    hold on
end
% linkaxes([ax{1} ax{2}],'x')

xlim(TimeRange);
% ax1.XLim = [5821.5 5822.5];
title(['maze ' num2str(maze) ', lap ' num2str(lap)]);
f.WindowState = 'maximized';
zoom reset