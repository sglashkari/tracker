clc; close all;
load('200629_test_45fps_noRecording.mat');

x1 = c(:,1);
y1 = c(:,2);
x2 = c(:,3);
y2 = c(:,4);

subplot(2,1,1)
plot(t,x1,'r.')
subplot(2,1,2)
plot(t,x2,'b.')

figure(2)
x1n = normalize(x1);
x2n = -1*normalize(x2);
x1m = x1n - mean(x1n(t<4));
x2m = x2n - mean(x2n(t<4));
plot(t,x1m,'r',t,x2m,'b')


figure(3)
[c,lags] = xcorr(x1m,x2m);
plot(lags,c)
[argvalue, argmax] = max(c);
lags(argmax)/fr

% 0.139752 sec (for 200629)