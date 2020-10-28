function track_led2(mat_filename)
%% TRACK_LED2 plots the position from the mat file
%
% Input:
%   mat file of jumping experiment (Rat 913)
%
% Output:
%   plots of position vs time and x-y
%
%   See also: track_led0, track_led1
%
% Author: Shahin G Lashkari
% Date: October 26, 2020
%
clc;
close all;
if nargin==0
    [file, path] = uigetfile('D:\OneDrive - Johns Hopkins\JHU\913_Jumping_Recording\2020-10-25_Rat913-01\Videos\2020-10-25_15-59-09.mat');
    mat_filename = fullfile(path,file);
end
load(mat_filename,'position');
k = (position(:,2) > 0);
t = position(k,1)/30;
x = position(k,2);
y = position(k,3);

%% plot raw data
figure(1)
tiledlayout(2,1)
ax1 = nexttile;
plot(t,x,'.')
ax2 = nexttile;
plot(t,y,'.')
linkaxes([ax1 ax2],'x')
xlim([t(1) t(end)])

figure(2)
plot(x,y,'.')
axis equal

