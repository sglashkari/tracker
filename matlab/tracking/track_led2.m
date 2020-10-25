function track_led2(vid_filename)
%% CIRCLE2  calculates latency of dome for methods paper
%
% Input:
%   mat file calculated by circle.
%
% Output:
%   latency
%
% Method:
%   cross correlation and interpolation
%   
%   See also circle1, circle0, xcorr.
%
% Author: Shahin G Lashkari
% Date: August 1, 2020
%
clc;
close all;
if nargin==0
    vid_filename = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\913_Jumping_Recording\2020-10-23-Day0\2020-10-23-133050.mp4';
end
load(strrep(vid_filename,'mp4','mat'));
varargin
k = (position(:,2) > 0);
t = position(k,1)/10;
x = position(k,2);
y = position(k,3);

%% plot raw data
figure(1)
plot(t,x,t,y)
figure(2)
plot(x,y)

