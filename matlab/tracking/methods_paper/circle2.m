function circle2(vid_filename)
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
    vid_filename = 'D:\domeLatencyMeasurement\200731_latencyTest\200731_latencyTest_90fps.avi';
end
load(strrep(vid_filename,'avi','mat'));
varargin 
x1 = c(:,1); % t<12
y1 = c(:,2);
x2 = c(:,3); %t<12
y2 = c(:,4);
t =t(:); % t<12

%% plot raw data
figure(1)
subplot(2,1,1)
plot(t,x1,'r.')
subplot(2,1,2)
plot(t,x2,'b.')

%% plot normalized data
figure(2)
x1n = normalize(x1);
x2n = -1*normalize(x2);
x1m = x1n - nanmean(x1n(t<4));
x2m = x2n - nanmean(x2n(t<4));
plot(t,x1m,'r',t,x2m,'b')

%% cross correlation
figure(3)
[cc,lags] = xcorr(x1m,x2m);
stem(lags,cc)
% figure(5)
% plot(lags,cc)
[~, argmax] = max(cc);
lags(argmax)/frame_rate
x_fit = round(linspace(argmax-4,argmax+4,9));
y_fit = cc(x_fit);

p = polyfit(x_fit,y_fit,2);
x_fit2 = linspace(argmax-4,argmax+4,900);
y_fit2 = polyval(p,x_fit2);

[~, argmax] = max(y_fit2);
lags2 = interp1(x_fit,lags(x_fit),x_fit2(argmax));
lags2/frame_rate

figure(4)
plot(x_fit,y_fit,'ob',x_fit2,y_fit2,'.r')

figure(2)
end
% 0.139752 sec (for 200629)
      
% 0.1320 200704_latencyTest_30fps_noRecording_video [C]         0.1357      
% 0.1320 200704_latencyTest_45fps_noRecording_video * ---       0.1310
% 0.1242 200704_latencyTest_60fps_noRecording_video *           0.1230
% 0.1320 200704_latencyTest_75fps_noRecording_video [C]         0.1289
% 0.1320 200704_latencyTest_90fps_noRecording_video [C]         0.1302
% 0.1398 200704_latencyTest_30fps_yesRecording_video [C]        0.1401
% 0.1242 200704_latencyTest_45fps_yesRecording_video [C]        0.1256
% 0.1320 200704_latencyTest_60fps_yesRecording_video [C]        0.1293
% 0.1242 200704_latencyTest_75fps_yesRecording_video [C][F]     0.1223
% 0.1320 200704_latencyTest_90fps_yesRecording_video [C][F]     0.1333

% 0<t<15
% === no  ===== yes (+-0.0078)
% 0.1320        0.1398
% 0.1242        0.1242
% 0.1242        0.1242
% 0.1242        0.1165
% 0.1320        0.1320

% 0<t<15 (interpolation)
% === no  ===== yes (+-0.0078)
% 0.1355        0.1400
% 0.1280        0.1205
% 0.1236        0.1242
% 0.1210        0.1165
% 0.1301        0.1333

%% 2020-07-31 after c++ code improvement
% 30fps         0.1103  0.1100
% 45fps         0.1000  0.1011
% 60fps         0.1026  0.1032
% 75fps         0.1026  0.1033
% 90fps         0.1076  0.1074




