%% Plot all data in search for spikes during jump
% Shahin G Lashkari August 10th, 2020
clc; 
%clear; close all;

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
lap = 2;

TimePosF = TimePos(~isnan(AngularPosition));
AngularPositionF = AngularPosition(~isnan(AngularPosition));
AngularPositionF = AngularPositionF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));
TimePosF =  TimePosF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));

LapTime = TimePosF(diff(AngularPositionF) < -340);
TimeRange = [LapTime(lap) LapTime(lap+1)]; % seconds


%% TT
TimeTT=[];
for i=1:N
	TTT = readtt(convertStringsToChars(absolue_paths(i)), TimeRange * 1e6); % microseconds
    TimeTT = [TimeTT; TTT];
end

%% CSC
FilenameCSC = fullfile(Nlx_directory,'CSC14.ncs');
[TimeCSC,Data] = readcsc(FilenameCSC, TimeRange * 1e6); % microseconds

%% Plotting
f = figure(2);
tiledlayout(2,1)
ax1 = nexttile;
rectangle('Position',[TimeRange(1),160,diff(TimeRange),40],'FaceColor','g','EdgeColor','g')
hold on
rectangle('Position',[TimeRange(1),100,diff(TimeRange),10],'FaceColor','c','EdgeColor','c')
hold on
rectangle('Position',[TimeRange(1),260,diff(TimeRange),10],'FaceColor','c','EdgeColor','c')
hold on
plot(TimePosF,AngularPositionF,'LineWidth',4)
ylim([0 360])
pan xon
zoom xon
ylabel('Angular Position (\circ)')
title(['maze ' num2str(maze)])

% ax2 = nexttile;
% [counts5,edges] = histcounts(TimeTT,'BinWidth',5e-3); % 5 milliseconds
% edges50 = movmean(edges,11); % time of the middle of the bin
% counts50 = movsum(counts5,10); % add 
% bar(edges50(6:end-5),counts50(5:end-5))
% pan xon
% zoom on
% ylabel('Frequency')

ax3 = nexttile;
plot(TimeCSC,Data*1e6);
ylim([-1200 1200])
xlabel('Time (sec)')
ylabel('LFP (\muV)')
pan xon
zoom xon
linkaxes([ax1 ax3],'x')

ax1.XLim = TimeRange + [35e-3 -15e-3];
%ax1.XLim = [5780 5786];
% ax1.XLim = [5820 5826]; % m4l4
ax1.Title.String = ['maze ' num2str(maze) ', lap ' num2str(lap)];
%f.WindowState = 'maximized';


set(f,'Units','Inches');
width = 10;
height = 5;
set(f,'position',[0,0,width,height])
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'lfp','-dpdf','-r0')