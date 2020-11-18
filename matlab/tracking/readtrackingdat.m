%function [Time,X,Y, Dir, Header] = readtrackingdat(Filename, StartTime, EndTime)
%%READBINARY reads pos.p file with the header
% StartTime is the vector of initial times (one value for each maze)
% EndTime is the vector of final times (one value for each maze)
% Times are in seconds
% sgl 2020-08-30

% if nargin == 0
%     Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
%     Filename = fullfile(Nlx_directory,'pos.p');
% end
% if nargin < 3
%     TimeEV = readevent(fullfile(Nlx_directory,'Events.nev'));
%     StartTime = TimeEV(1:2:end);
%     EndTime = TimeEV(2:2:end);
% end
clear;
close all;
clc
Filename = '~/Downloads/tracking.dat'; %'~/flycapture/bin/tracking.dat';
FileID = fopen(Filename,'r');
%Header = arrayfun(@(~) fgetl(FileID), 1:24, 'UniformOutput', false)';
%header_length = ftell(FileID);
format(:,1) = {'double'; 'int32'; 'int32'; 'int32'; 'int32'; 'single'; 'single'};
format(:,2) = mat2cell(ones(7,2),ones(1,7));
format(:,3) = split('time,p1,p2,p3,p4,x,y',",");

% format(:,1) = {'uint32';'uint32';'double'; 'uint32'; 'uint32'; 'uint32'; 'uint32'; 'single'; 'single'};
% format(:,2) = mat2cell(ones(9,2),ones(1,9));
% format(:,3) = split('frame,flag,time,p1,p2,p3,p4,x,y',",");

m = memmapfile(Filename,'Format',format); %, 'Offset', header_length);
fclose(FileID);

Samples = m.Data;
x = [Samples.x]';
y = [Samples.y]';
p3 = [Samples.p3]';
%flag = [Samples.flag]';
time = [Samples.time]'; %seconds
plot((time-64)/64*pi)

hold on
t = unwrap((time-64)/64*pi)/pi*64;
plot(t)
figure(2)
plot(t,x,'.')
% % Time Range
% idx = [];
% for i=1:length(StartTime)
%     idx = [idx; find(Time >= StartTime(i) & Time <= EndTime(i))];
% end
% 
% % If no output only plot the results:
% if nargout == 0
%     close all
%     figure(1)
%     plot(X, 480 - Y,'.');
%     axis equal
%     axis([0 640 0 480])
%     figure(2)
%     plot(Time, Y,'.',Time(idx), Y(idx),'*');
%     clear Time
% else
%     X = X(idx);
%     Y = Y(idx);
%     Dir = Dir(idx);
%     Time = Time(idx);
% end
% 
% end
figure(3)
plot(diff(t))
(t(end)-t(1))/60
figure(4)
plot(p3)