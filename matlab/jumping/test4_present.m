% SGL 2021-02-08
clc; clear; close all
addpath('../jumping');
exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax');
x_thresh = 84;

%lap_no = 21;
cluster_no = 31;

direction = 'left';
N = nnz([lap.dir]==direction);
c = cluster_no;
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

theta_all = [];
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
    plot(pos.t(idx) - lap(l).t_jump_exact,pos.s(idx),'.b');
    plot(pos.t(idx) - lap(l).t_jump_exact,pos.filt.s(idx),'r');
    plot(pos.t(idx) - lap(l).t_jump_exact, 40*ones(size(pos.t(idx))),'g')
    %plot([lap(l).t_jump lap(l).t_jump], [0 300])
    %plot([lap(l).t_jump_exact lap(l).t_jump_exact], [0 300]);
    ylabel('Speed (cm/s)')
    box on
    %ylim([0 300])
    
 
    ax3 = subplot(4,1,3); hold on;
    plot(timecsc - lap(l).t_jump_exact,lfp * 1e6,'Color','#D0D0D0')
    plot(timecsc - lap(l).t_jump_exact,theta *1e6,'r');
    ylim([-700 700])
    ylabel('Theta (\muV)')
    t_idx = (timecsc > lap(l).t_jump_exact - 2) & (timecsc < lap(l).t_jump_exact + 2);% 1 sec before and 1 sec after

    theta_all(i,:) = theta(t_idx)';
    lfp_all(i,:) = lfp(t_idx)';
    t_theta_ave = timecsc(t_idx) - lap(l).t_jump_exact;
    box on
    
    ax4 = subplot(4,1,4); hold on;
    plot(timecsc - lap(l).t_jump_exact,(theta-lfp)*1e6)
    xlabel('Time (sec)')
    ylabel('Error (\muV)')
    box on
    
    linkaxes([ax2 ax3 ax4],'x')
    xlim([-2 2])
end

% average thetas
theta_ave = mean(theta_all);
lfp_ave = mean(lfp_all);
ax1 = subplot(4,1,1); hold off; hold on;
plot(t_theta_ave,lfp_ave*1e6,'Color','#808080');
%plot(t_theta_ave,theta_ave*1e6,'r'); 
ylabel('Average of all lfps (\muV)')
theta_lft_ave = filtertheta(t_theta_ave,lfp_ave);
%plot(t_theta_ave,theta_lft_ave*1e6,'r')
ylim([-400 400])
box on

sgtitle([direction 'ward laps aligned']);

saveas(gcf,fullfile(exp_directory, 'Analysis',['Alligned-' direction 'ward.svg']))