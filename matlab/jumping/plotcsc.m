%%% Plotting CSCs
clc; close all
exp_directory = '/home/shahin/Desktop/20-12-09';

addpath('../jumping');
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'lap');
l = 2;
cscs = 2:16;
for c=cscs
    
    csc_filename= fullfile(exp_directory,'Neuralynx',['CSC' num2str(c) '.ncs']);
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase, mag] = filtertheta(timecsc,lfp);
    
    ax1 = subplot(3,1,1);
    %plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc,theta *1e6);%,'r');
    ylim([-300 300])
    ylabel('Theta (\muV)')
    title(['lap ' num2str(l)])
    
    ax2 = subplot(3,1,2);
    hold on
    plot(timecsc,phase);%,'b');
    ylim([-200 200])
    ylabel('Phase (Angle)')
    
    ax3 = subplot(3,1,3);
    hold on
    plot(timecsc,mag*1e6);%,'b');
    %ylim([-200 200])
    ylabel('Magnitude')
    
    
    disp(['mean mag for ' num2str(c) ' is ' num2str(mean(mag) *1e6,'%.2f')]);
end
legend(num2str(cscs'))
linkaxes([ax1 ax2 ax3],'x')

%%
figure(2)
cscs = 11:12;
for c=cscs
    
    csc_filename= fullfile(exp_directory,'Neuralynx',['CSC' num2str(c) '.ncs']);
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase, mag] = filtertheta(timecsc,lfp);
    
    ax1 = subplot(3,1,1);
    %plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc,theta *1e6);%,'r');
    ylim([-300 300])
    ylabel('Theta (\muV)')
    title(['lap ' num2str(l)])
    
    ax2 = subplot(3,1,2);
    hold on
    plot(timecsc,phase);%,'b');
    ylim([-200 200])
    ylabel('Phase (Angle)')
    
    ax3 = subplot(3,1,3);
    hold on
    plot(timecsc,mag*1e6);%,'b');
    %ylim([-200 200])
    ylabel('Magnitude')
    
end
legend(num2str(cscs'))
linkaxes([ax1 ax2 ax3],'x')