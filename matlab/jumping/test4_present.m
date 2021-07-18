% SGL 2021-02-08
clc; clear; 
close all
addpath('../jumping');
exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax');
direction = 'right';

%% CSC
csc_filename= fullfile(exp_directory,'Neuralynx','CSC12.ncs');

% filtered data:
pos.filt.vx = filtertheta(pos.t, pos.vx, 0.01, 10);
pos.filt.vy = filtertheta(pos.t, pos.vy, 0.01, 10);
pos.filt.s = vecnorm([pos.filt.vx pos.filt.vy]')'; % speed in cm/sec
pos.filt.ax = gradient(pos.filt.vx)./gradient(pos.t); % ax in cm/sec
pos.filt.ay = gradient(pos.filt.vy)./gradient(pos.t); % ax in cm/sec

for l =1:length(lap)
    % detect the more exact time of jump
    idx = (pos.t>lap(l).t_jump-0.6) & (pos.t<lap(l).t_jump+0.4) & (pos.filt.s>40);
    idx = find(idx,1);
    lap(l).t_jump_exact = pos.t(idx);
end

f=figure(100);
f.WindowState = 'maximized';

imu_a_all = [];
imu_w_all = [];
lfp_all = [];
i = 0;
for l=1:length(lap)
    if lap(l).dir~=direction % filtering out undesired direction
        continue;
    end
    i = i + 1;
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase, mag] = filtertheta(timecsc,lfp);
    ax1 = subplot(4,1,1); hold on;
    ax2 = subplot(4,1,2); hold on;
    idx = pos.lap==l;
    %plot(pos.t(idx) - lap(l).t_jump_exact,pos.s(idx),'.b');
    plot(pos.t(idx) - lap(l).t_jump_exact,pos.filt.vx(idx));%,'b');
    %plot(pos.t(idx) - lap(l).t_jump_exact, 40*ones(size(pos.t(idx))),'k');
    %plot([lap(l).t_jump lap(l).t_jump], [0 300])
    %plot([lap(l).t_jump_exact lap(l).t_jump_exact], [0 300]);
    ylabel('Speed (cm/s)')
    box on
    %ylim([0 300])
    
 
    ax3 = subplot(4,1,3); hold on;
    %plot(timecsc - lap(l).t_jump_exact,lfp * 1e6,'Color','#D0D0D0')
    plot(timecsc - lap(l).t_jump_exact,theta *1e6);%,'r');
    ylim([-700 700])
    ylabel('Theta (\muV)')
    t_idx = (timecsc > lap(l).t_jump_exact - 2) & (timecsc < lap(l).t_jump_exact + 2); % 2 sec before and 2 sec after

    size(lfp(t_idx))
    lfp_all(i,:) = lfp(t_idx)';
    t_lfp_ave = timecsc(t_idx) - lap(l).t_jump_exact;
    box on
    
    ax4 = subplot(4,1,4); hold on;
    imu = readimu(fullfile(exp_directory,'Neuralynx'),lap(l).t * 1e6); % microseconds
    
    t_idx = (imu.t > lap(l).t_jump_exact - 2) & (imu.t < lap(l).t_jump_exact + 2); % 2 sec before and 2 sec after
    imu.filt.a = filtertheta(imu.t, imu.a, 0.01, 10);
    imu_a_all(i,:) = imu.a(t_idx)';
    imu_w_all(i,:) = imu.w(t_idx)';
    t_imu_ave = imu.t(t_idx) - lap(l).t_jump_exact;
    plot(imu.t - lap(l).t_jump_exact,imu.ax*1e6)
    xlabel('Time (sec)')
    ylabel('IMU Lin Acc (\muV)')
    box on
    
    linkaxes([ax2 ax3 ax4],'x')
    xlim([-2 2])
end

% average thetas
imu_a_ave = mean(imu_a_all);
imu_w_ave = mean(imu_w_all);
lfp_ave = mean(lfp_all);
ax1 = subplot(4,1,1); hold off; hold on;
plot(t_lfp_ave,lfp_ave*1e6,'Color','#808080')
plot(t_imu_ave,imu_a_ave*1e6,'r');
plot(t_imu_ave,imu_w_ave*1e6,'g');
legend('LFP','IMU Lin Acc','IMU Ang Vel')
ylabel('Average of all laps (\muV)')

% lfp_ave_filt = filtertheta(t_lfp_ave,lfp_ave);
% plot(t_lfp_ave,lfp_ave_filt*1e6,'g');
% 
% imu_a_ave_filt = filtertheta(t_imu_ave,imu_a_ave);
% plot(t_imu_ave,imu_a_ave_filt*1e6,'m');

linkaxes([ax1 ax2 ax3 ax4],'x')
%ylim([-400 400])
box on

sgtitle([direction 'ward laps aligned']);

saveas(gcf,fullfile(exp_directory, 'Analysis',['Alligned-' direction '.svg']))