% SGL 2021-02-08
clc; close all
addpath('../jumping');
exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax');
x_thresh = 84;

%lap_no = 21;
cluster_no = 31;

direction = 'right';
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
i = 0;
for l=1:length(lap)
%     if lap(l).dir~=direction % filtering out undesired direction
%         continue;
%     end
    i = i + 1;
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase] = filtertheta(timecsc,lfp);
    
    idx = pos.lap==l;
    ax1 = subplot(5,1,1); hold on;
    %plot(pos.t(idx),pos.x(idx));
    %ylabel('Horizontal position (cm)')
    %plot(pos.t(idx) - lap(l).t_jump_exact,pos.filt.ax(idx));
    %ylabel('Horizontal acceleration (cm/s^2)')
    title([direction 'ward laps aligned']);
    
    ax2 = subplot(5,1,2); hold on;
    plot(pos.t(idx) - lap(l).t_jump_exact,pos.s(idx),'.b');
    plot(pos.t(idx) - lap(l).t_jump_exact,pos.filt.s(idx),'r');
    plot(pos.t(idx) - lap(l).t_jump_exact, 40*ones(size(pos.t(idx))),'g')
    plot([lap(l).t_jump lap(l).t_jump], [0 300])
    plot([lap(l).t_jump_exact lap(l).t_jump_exact], [0 300]);
    ylabel('Speed (cm/s)')
    %ylim([0 300])
    
    ax3 = subplot(5,1,3); hold on;
    plot(timecsc - lap(l).t_jump_exact, [0; diff(phase)<350])
    %idx = [cluster(c).lap]==l;
    %plot(cluster(c).t(idx) - lap(l).t_jump_exact, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    %ylim([-180 180])
    ylabel('Phase (Degrees)')
    
    
    ax4 = subplot(5,1,4); hold on;
    plot(timecsc - lap(l).t_jump_exact,lfp * 1e6,'Color','#D0D0D0')
    plot(timecsc - lap(l).t_jump_exact,theta *1e6,'r');
    ylim([-600 600])
    ylabel('Theta (\muV)')
    t_idx = (timecsc > lap(l).t_jump_exact - 1) & (timecsc < lap(l).t_jump_exact + 1);% 1 sec before and 1 sec after
     
    dt = length(theta(t_idx));
    theta_all(i,:) = theta(t_idx)';
    
    ax5 = subplot(5,1,5); hold on;
    plot(timecsc - lap(l).t_jump_exact,(theta-lfp)*1e6)
    xlabel('Time (sec)')
    ylabel('Error (\muV)')
    
    linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
    %xlim(lap(l).t);
    %zoom xon
    xlim([-1 1])
    
end

% average thetas
theta_ave = mean(theta_all);
t_theta_ave = timecsc(t_idx) - lap(l).t_jump_exact;
ax1 = subplot(5,1,1); hold off; hold on;
plot(t_theta_ave,theta_ave*1e6,'r'); % divide by number of laps
ylabel('Average of all thetas  (\muV)')
[~, phase_ave] = filtertheta(t_theta_ave,theta_ave);
%plot(t_theta_ave,phase_ave)
for l=1:length(lap)
    [~, phase_all(l,:)] = filtertheta(t_theta_ave,theta_all(l,:));
    phase_diff(l,:) = unwrap(angdiff(phase_all(l,:)*pi/180,phase_ave*pi/180))*180/pi;
end
plot(t_theta_ave, mean(phase_diff))
saveas(gcf,fullfile(exp_directory, 'Analysis',['Alligned-' direction 'ward-cl_' num2str(cluster_no) '.jpg']))