%clc; clear; close all;
% % addpath('..'); %adding readbinary to the path
% % file_name = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx\TT14\TT14.ntt.parms';
% % [cluster, timetable, table, header] = readnttparam(file_name);
% % size(data)
% 
% %% these are the format of the file for ntt.param
% % format(:,1) = mat2cell(['double'; repmat('single',[26 1])],ones(1,27));
% % format(:,2) = mat2cell(ones(27,2),ones(1,27));
% % format(:,3) = split('Timestamp,MPeakX,MPeakY,MPeakA,MPeakB,PreVallEyX,PreValleyY,PreValleyA,PreValleyB,SpikeID,EnergyX,EnergyY,EnergyA,EnergyB,MaxHeight,MaxWidth,XPos,YPos,Time,PeakX,PeakY,PeakA,PeakB,ValleyX,ValleyY,ValleyA,ValleyB',",");
% 
% % test line
% 
% subplot(3,1,3)
% x = [0  89  89  112 112 242 242 269 269 360];
% y = [0  0   -3  -3  0   0   -3  -3  0   0];
% line(x,y,'LineWidth',2)
% ylim([-4 1])
% xlim([0 360])
% 
% t=0:360;
% data = sin(0.1*t);%rand(1,360); % sample data
% subplot(3,1,1)
% bar(data)
% subplot(3,1,2)
% map = interp1([0;1],[1 1 1;0 0.45 0.74],0:0.01:1); % color map
% imagesc(data)
% colormap('parula')
% cb=colorbar;
% cb.Position = cb.Position + [cb.Position(1)*0.12 0 -cb.Position(3)*0.3 0];

% f1 = figure(1);
% % ax1 = subplot(2,1,1);
% t=0:0.1:20;
% x=sin(t);
% 
% %
% % ax1.Units = 'pixels';
% % ax1.Position = [50 50 1550 50];
% % axes(ax1)
% ax1 = axes('units','pixels','position',[100 200 1536 100]);
% plot(t,x)
% 
% % set(gca, 'units','pixels','outerposition',[100 200 1550 100]);
% ax2 = axes('units','pixels','position',[100 350 1536 100]);
% plot(t,x.^2)
% % xlim([0 16])
% % set(gca, 'units','pixels','outerposition',[100 350 1550 100]);
% ax3 = axes('units','pixels','position',[100 500 1536 750]);
% imshow(2*I+50)
% % linkaxes([ax1 ax2 ax3],'x')
% % xlim([0 16])
% 
% set(f1, 'units','pixels','position',[50 50 1800 1500]);

% Nlx_directory = '/home/shahin/onedrive/JHU/913_Jumping_Recording/2020-10-25_Rat913-01/Neuralynx';
% Nlx_directory = '/home/shahin/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02/Neuralynx';
% Nlx_directory = '/home/shahin/onedrive/JHU/913_Jumping_Recording/2020-11-22_Rat913-03/Neuralynx';
% [Time,Data,Header,EventIDs,TTLs] = readevent(fullfile(Nlx_directory,'Events.nev'));
% Data(10:end)

exp_directory = '~/Desktop/20-12-09';
%exp_directory = '~/Desktop/20-12-21';
% [datafile,exp_directory] = uigetfile(fullfile(exp_directory,'analyzed_data.mat'), 'Select Data File');
% if isequal(datafile, 0)
%     error('Data file was not selected!')
% end
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'spike','ppcm', 'offset', 'colors','xmax','x_thresh');
%%
j = 1; % cluster_no;

f = figure(1000);
clf(f);
for l=2:2:30
ax(round(l/2)) = subplot(15,1,round(l/2));
idx = [spike(j).lap]==l;
plot(spike(j).x(idx),spike(j).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
hold on
plot(spike(j).x(idx), 360 + spike(j).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
%title(['cluster ' num2str(j) ' lap ' num2str(l)])
ylabel('Phase (angle)')
%xlabel('Horizontal position (cm)')
zoom xon
hold on
end
xlabel('Horizontal position (cm)')
zoom xon
linkaxes(ax,'x')
xlim([0 xmax])

%%
f = figure(1001);
clf(f);
for l=2:2:30
ax(round(l/2)) = subplot(15,1,round(l/2));
idx = [spike(j).lap]==l;
x_offset = (spike(j).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
X = spike(j).x(idx) + x_offset;
phase = spike(j).phase(idx);
plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
hold on
plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
%title(['cluster ' num2str(j) ' lap ' num2str(l)])
ylabel('Phase (angle)')
%xlabel('Horizontal position (cm)')
zoom xon
hold on
end
xlabel('Horizontal position (cm)')
zoom xon
linkaxes(ax,'x')
xlim([0 xmax])

%%
f = figure(1002);
clf(f);
for l=2:2:30
idx = [spike(j).lap]==l;
x_offset = (spike(j).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
X = spike(j).x(idx) + x_offset;
phase = spike(j).phase(idx);
plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
hold on
plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
ylabel('Phase (angle)')
end
title(['cluster ' num2str(j)])
ylabel('Phase (angle)')
xlabel('Horizontal position (cm)')
xlim([0 xmax])

%%
f = figure(1003);
clf(f);
for l=16:2:30
ax(round(l/2)-7) = subplot(8,1,round(l/2)-7);
idx = [spike(j).lap]==l;
x_offset = (spike(j).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
X = spike(j).x(idx) + x_offset;
phase = spike(j).phase(idx);
plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
hold on
plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
title(['cluster ' num2str(j) ' lap ' num2str(l)])
ylabel('Phase (angle)')
end
xlabel('Horizontal position (cm)')
linkaxes(ax,'x')
xlim([0 xmax])
