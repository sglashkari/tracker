clc; clear; close all;
% addpath('..'); %adding readbinary to the path
% file_name = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx\TT14\TT14.ntt.parms';
% [cluster, timetable, table, header] = readnttparam(file_name);
% size(data)

%% these are the format of the file for ntt.param
% format(:,1) = mat2cell(['double'; repmat('single',[26 1])],ones(1,27));
% format(:,2) = mat2cell(ones(27,2),ones(1,27));
% format(:,3) = split('Timestamp,MPeakX,MPeakY,MPeakA,MPeakB,PreVallEyX,PreValleyY,PreValleyA,PreValleyB,SpikeID,EnergyX,EnergyY,EnergyA,EnergyB,MaxHeight,MaxWidth,XPos,YPos,Time,PeakX,PeakY,PeakA,PeakB,ValleyX,ValleyY,ValleyA,ValleyB',",");

% test line

subplot(3,1,3)
x = [0  89  89  112 112 242 242 269 269 360];
y = [0  0   -3  -3  0   0   -3  -3  0   0];
line(x,y,'LineWidth',2)
ylim([-4 1])
xlim([0 360])

t=0:360;
data = sin(0.1*t);%rand(1,360); % sample data
subplot(3,1,1)
bar(data)
subplot(3,1,2)
map = interp1([0;1],[1 1 1;0 0.45 0.74],0:0.01:1); % color map
imagesc(data)
colormap('parula')
cb=colorbar;
cb.Position = cb.Position + [cb.Position(1)*0.12 0 -cb.Position(3)*0.3 0];

