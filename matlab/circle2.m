clc; 
close all;
load('D:\domeLatencyMeasurement\200704_latencyTest\200704_latencyTest_45fps_yesRecording_video.mat');

x1 = c(t<15,1);
y1 = c(:,2);
x2 = c(t<15,3);
y2 = c(:,4);
t =t(t<15);

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
x1m = x1n - mean(x1n(t<4));
x2m = x2n - mean(x2n(t<4));
plot(t,x1m,'r',t,x2m,'b')

%% cross correlation
figure(3)
[c,lags] = xcorr(x1m,x2m);
stem(lags,c)
% figure(4)
% plot(lags,c)
[argvalue, argmax] = max(c);
lags(argmax)/frame_rate
x_fit = round(linspace(argmax-4,argmax+4,9));
y_fit = c(x_fit);

p = polyfit(x_fit,y_fit,2);
x_fit2 = linspace(argmax-4,argmax+4,900);
y_fit2 = polyval(p,x_fit2);

[argvalue, argmax] = max(y_fit2);
lags2 = interp1(x_fit,lags(x_fit),x_fit2(argmax));
lags2/frame_rate

figure(4)
plot(x_fit,y_fit,'ob',x_fit2,y_fit2,'.r')

figure(2)
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


