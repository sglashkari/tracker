% SGL 2021-02-15
clc; clear; close all
% exp_directory = 'D:\Analysis\2021-12-10';
% mat_filename = fullfile(exp_directory,'analyzed_data.mat');
% load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'colors','xmax','hist');
exp_directory = 'E:\Rat1068\Analysis\2022-12-20\';
% exp_directory = 'E:\Rat1055\Analysis\2022-11-09\';
% exp_directory = 'E:\Rat980\Analysis\2021-12-21\';
load(fullfile(exp_directory,'processed_data.mat'),'pos','posi', 'lap', 'cluster', 'colors','hist', 'csc' ,'exp');
% xmax = exp.xmax;
% ppcm = exp.ppcm;
%%
direction = 'left';
timerange = [-4 2]; % 4 sec before to 2 sec after

%% CSC
f=figure(100); clf;
f.WindowState = 'maximized';

lfp_all_jump = [];
lfp_all_ditch = [];

j = 0;
d = 0;

for l=[lap([lap.dir]==direction).no]  % only the desired direction
    %lap(l).t_jump_exact = lap(l).t_jump_exact + rand; % test!!!!
    if l == numel(lap)
        continue;
    end
    
    
    if lap(l).status == "jump"
        j = j + 1;
        idx = (csc.t >= lap(l).t_jump_exact + timerange(1));% & (time < lap(l).t_jump_exact + 2);
        idx = find(idx):find(idx)+30000*(diff(timerange));
        lfp_all_ditch(j,:) = csc.lfp(idx)';
        ax(2) = subplot(4,1,2); hold on;
        
    else
        d = d + 1;
        idx = (csc.t >= lap(l).t_jump_exact + timerange(1));% & (time < lap(l).t_jump_exact + 2);
        idx = find(idx):find(idx)+30000*(diff(timerange));
        lfp_all_jump(d,:) = csc.lfp(idx)';
        ax(4) = subplot(4,1,4); hold on;
    end
    
    idx = (posi.t >= lap(l).t_jump_exact + timerange(1)) & (posi.t < lap(l).t_jump_exact + timerange(2));
    plot(posi.t(idx)-lap(l).t_jump_exact,abs(posi.filt.vx(idx)),'m');
    
    idx = (pos.t >= lap(l).t_jump_exact + timerange(1)) & (pos.t < lap(l).t_jump_exact + timerange(2));
    plot(pos.t(idx)-lap(l).t_jump_exact,50*ones(size(pos.t(idx))),'k');
    ylabel('Horizontal Speed (cm/s)')
    ylim([0 300])
    box on
    
    
    xlabel('Time (sec)')

    box on
    xlim(timerange)
end

%%
% average thetas
lfp_ave_jump = median(lfp_all_jump);
lfp_ave_ditch = median(lfp_all_ditch);
idx = (csc.t >= lap(1).t_jump_exact + timerange(1));
idx = find(idx):find(idx)+30000*(diff(timerange));
timecsc = csc.t(idx);
t_lfp_ave = timecsc - lap(1).t_jump_exact; % just want to have the resolution of time: dt [-2:dt:2]
ax(1) = subplot(4,1,1); hold off; hold on;
[theta_ave_jump] = filterlfp(t_lfp_ave,lfp_ave_jump);
plot(t_lfp_ave,lfp_ave_jump*1e6,'Color','#808080')
plot(t_lfp_ave,theta_ave_jump *1e6,'r');
%legend('LFP','IMU Lin Acc','IMU Ang Vel')
ylabel('Average of all laps (\muV)')
title('jump')

ax(3) = subplot(4,1,3); hold off; hold on;
[theta_ave_ditch] = filterlfp(t_lfp_ave,lfp_ave_ditch);
plot(t_lfp_ave,lfp_ave_ditch*1e6,'Color','#808080')
plot(t_lfp_ave,theta_ave_ditch *1e6,'r');
%legend('LFP','IMU Lin Acc','IMU Ang Vel')
ylabel('Average of all laps (\muV)')
title('ditch')
linkaxes(ax,'x')
xlim(timerange)


ylim([-300 300])
box on

sgtitle([direction 'ward laps aligned']);
saveas(gcf,fullfile(exp_directory, 'Analysis',['Alligned-' direction '.jpg']))
saveas(gcf,fullfile(exp_directory, 'Analysis',['Alligned-' direction '.fig']))
% % Histogram of phases at time zero
% figure(200)
% bin_size = 30; % deg
% edges = -180:bin_size:180; % deg
% counts = histcounts(phase_all, edges);
% histogram('BinCounts', counts, 'BinEdges', edges, 'FaceColor','red');
% %title(['cluster ' mat2str(c) ' all laps combined - direction: ' direction 'ward ']);
% ylabel('Count');
% xlabel('Phase (deg)');
% xlim([-180 180]);