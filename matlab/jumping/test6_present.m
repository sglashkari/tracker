% SGL 2021-02-15
clc; clear; close all
exp_directory = 'D:\Analysis\2021-12-21';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'colors','xmax','hist');
direction = 'right';
timerange = [-4 2]; % 2 sec before to 2 sec after

%% CSC
f=figure(100);
f.WindowState = 'maximized';

lfp_all = [];
phase_all = [];
i = 0;

sh=1;
csc_filename=fullfile(exp_directory, ['CA1-Shank' num2str(sh)],'B','LFP32.ncs'); % use plotcsc to optimize it
[time,data] = read_bin_csc(csc_filename);

for l=[lap([lap.dir]==direction).no]  % only the desired direction
    %lap(l).t_jump_exact = lap(l).t_jump_exact + rand; % test!!!!
    i = i + 1;
    idx = (time >= lap(l).t_jump_exact + timerange(1));% & (time < lap(l).t_jump_exact + 2);
    idx = find(idx):find(idx)+30000*(diff(timerange));
    timecsc = time(idx);
    lfp = data(idx);
    [theta, phase, mag] = filtertheta(timecsc,lfp);
    ax2 = subplot(4,1,2); hold on;
    
    idx = (posi.t >= lap(l).t_jump_exact + timerange(1)) & (posi.t < lap(l).t_jump_exact + timerange(2));
    plot(posi.t(idx)-lap(l).t_jump_exact,abs(posi.filt.vx(idx)),'m');
    
    idx = (pos.t >= lap(l).t_jump_exact + timerange(1)) & (pos.t < lap(l).t_jump_exact + timerange(2));
    plot(pos.t(idx)-lap(l).t_jump_exact,30*ones(size(pos.t(idx))),'k');
    ylabel('Horizontal Speed (cm/s)')
    ylim([0 300])
    box on
    
    ax3 = subplot(4,1,3); hold on;
    plot(timecsc - lap(l).t_jump_exact,lfp * 1e6,'Color','#D0D0D0')
    plot(timecsc - lap(l).t_jump_exact,theta *1e6,'r');
    ylim([-400 400])
    ylabel('Theta (\muV)')
    box on
    
    lfp_all(i,:) = lfp';
    
    
    ax4 = subplot(4,1,4); hold on
    [~,idx] = min(abs(timecsc-lap(l).t_jump_exact));
    plot(timecsc-lap(l).t_jump_exact,phase,'Color','#ADD8E6');
    plot(timecsc(idx)-lap(l).t_jump_exact,phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', 'red');

    ylim([-180 180])
    ylabel('Phase (Degrees)')
    
    phase_all(i) = phase(idx);
    
    xlabel('Time (sec)')

    box on
    
    linkaxes([ax2 ax3 ax4],'x')
    xlim(timerange)
end

% average thetas
lfp_ave = mean(lfp_all); % mean(lfp_all,1)
t_lfp_ave = timecsc - lap(l).t_jump_exact; % just want to have the resolution of time: dt [-2:dt:2]
ax1 = subplot(4,1,1); hold off; hold on;
[theta_ave, phase_ave, mag_ave] = filtertheta(t_lfp_ave,lfp_ave);
plot(t_lfp_ave,lfp_ave*1e6,'Color','#808080')
plot(t_lfp_ave,theta_ave *1e6,'r');
%legend('LFP','IMU Lin Acc','IMU Ang Vel')
ylabel('Average of all laps (\muV)')
xlim(timerange)

linkaxes([ax1 ax2 ax3 ax4],'x')
ylim([-200 200])
box on

sgtitle([direction 'ward laps aligned']);
saveas(gcf,fullfile(exp_directory, 'Analysis',['Alligned-' direction '.svg']))

%% Histogram of phases at time zero
% figure(200)
% bin_size = 30; % deg
% edges = -180:bin_size:180; % deg
% counts = histcounts(phase_all, edges);
% histogram('BinCounts', counts, 'BinEdges', edges, 'FaceColor','red');
% %title(['cluster ' mat2str(c) ' all laps combined - direction: ' direction 'ward ']);
% ylabel('Count');
% xlabel('Phase (deg)');
% xlim([-180 180]);