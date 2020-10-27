function track_led2(vid_filename)
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
    vid_filename = 'D:\OneDrive - Johns Hopkins\JHU\913_Jumping_Recording\2020-10-25_Rat913-01\Video\2020-10-25_15-59-09.mp4'; %day 1 
end
load(strrep(vid_filename,'mp4','mat'),'position');
varargin
k = (position(:,2) > 0);
t = position(k,1)/30;
x = position(k,2);
y = position(k,3);

%% plot raw data
figure(1)
plot(t,x,t,y)
figure(2)
plot(x,y,'.')
axis equal

