%%% Plotting CSCs
clc; close all
exp_directory = '/home/shahin/Desktop/20-12-09';

addpath('../jumping');
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'lap');
l = 2;
for c=[2 4]
    
    csc_filename= fullfile(exp_directory,'Neuralynx','CSC4.ncs');
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase] = filtertheta(timecsc,lfp);
    
    ax1 = subplot(2,1,1);
    plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc,theta *1e6,'r');
    ylim([-300 300])
    ylabel('Theta (\muV)')
    
    ax2 = subplot(2,1,2);
    hold on
    plot(timecsc,phase,'b');
    ylim([-10 370])
    ylabel('Phase (Angle)')
    
end

linkaxes([ax1 ax2],'x')