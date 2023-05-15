% SGL 2021-02-26
clc; clear; 
close all
addpath('../jumping');
exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax');

l = 2;

%% CSC
figure(1); 
ax(1) = subplot(4,1,1); hold on
csc_filename= fullfile(exp_directory,'Neuralynx','CSC12.ncs');
[timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
timecsc = timecsc - lap(l).t_jump;
[theta.data, ~, theta.mag] = filterlfp(timecsc,lfp,'theta');
[delta.data, ~, delta.mag] = filterlfp(timecsc,lfp,'delta');
[beta.data, ~, beta.mag] = filterlfp(timecsc,lfp,'beta');

plot(timecsc, lfp * 1e6, 'Color','#D0D0D0');
plot(timecsc, theta.data * 1e6,'r');
ylim([-400 400])

idx = (timecsc >= (lap(l).t_jump -2) & timecsc < (lap(l).t_jump + 2));
ax(2) = subplot(4,1,2); hold on
plot(timecsc, lfp * 1e6, 'Color','#D0D0D0');
plot(timecsc, beta.data * 1e6,'r');
xlabel('Beta ()')
ylim([-400 400])

n = 5000;
%timecsc = downsample(timecsc,n);
%lfp = downsample(lfp,n);
% ax(2) = subplot(4,1,2); hold on
% idx = pos.lap==l;
% plot(pos.t(idx),pos.filt.vx(idx));
% [v_filt] = filtertheta(pos.t(idx),pos.filt.vx(idx));
% plot(pos.t(idx),v_filt,'r');
% 
% %% sliding window
dt = 0.5; % 500 msec
% % L = length(timecsc); % size of total time
% % window.w = round (dt * L/(timecsc(end)-timecsc(1))); % width of the window
% % 
% % for k1 = 1:L-w+1
% %     window.idx = k1:k1+w-1;
% %     window.t = timecsc(mean(window.idx));
% %     window.lfp(k1,:) = lfp(window.idx);    
% %     
% % end
% 
% ax(3) = subplot(4,1,3);
% plot(timecsc, lfp-theta,'k');
% 
% idx = mag > 1.5e-4;
% 
% ax(4) = subplot(4,1,4); hold on
% plot(timecsc, mag,'k');
% plot(timecsc(idx), mag(idx),'*g');
% 
% subplot(4,1,1);
% plot(timecsc(idx), theta(idx),'*g');

ax(3) = subplot(4,1,3:4); hold on
%lfp = sin(10000*timecsc);
fs = round(1/(mean(diff(timecsc))));

% spectrogram(lfp,100,80,50,fs,'yaxis')

theta.zscore = zscore(movmean(theta.mag,fs*dt));
delta.zscore = zscore(movmean(delta.mag,fs*dt));


beta.zscore = zscore(movmean(beta.mag,fs*dt));
thetadelta.zscore = zscore(movmean(theta.mag./delta.mag,fs*dt));

plot(timecsc, theta.zscore)
plot(timecsc, delta.zscore)
ylim([-3 3])

idx = zscore(movmean(theta.mag,fs*dt)) > 0; % z-score > 0

subplot(4,1,1);
plot(timecsc(idx), theta.data(idx) * 1e6,'*g');


linkaxes(ax,'x')
xlim(lap(l).t)
xlim([(lap(l).t_jump -5) (lap(l).t_jump + 2)])
xlim([-5 2])

figure(2)
thetadelta.zscore = downsample(thetadelta.zscore, n);
beta.zscore = downsample(beta.zscore, n);
plot(thetadelta.zscore, beta.zscore, 'o')
