%%% Plotting CSCs
exp_directory = '~/Desktop/20-12-09';
csc_filename= fullfile(exp_directory,'Neuralynx','CSC4.ncs');
addpath('../jumping');
[timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
[theta, phase] = filtertheta(timecsc,lfp);
plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
hold on;
plot(timecsc,theta *1e6,'r');
ylim([-300 300])
ylabel('Theta (\muV)')
