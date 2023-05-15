% SGL 2021-01-31
clc; clear;
addpath('../jumping');
exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax');
colors(35) = "#00CCCC";
x_thresh = 84;
v_thresh = 20;

%lap_no = 10;
cluster_no = [25 35];
cluster_no = [6 11 23 16 35];
cluster_no = [6 23 16];
legendCell = cellstr(num2str(cluster_no', 'cluster #%-d'));

direction = 'right'; 
direction = 'left';
N = nnz([lap.dir]==direction);
isCSC = 1;

if isCSC
    close all;
end

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
    lap(l).t_jump = lap(l).t_jump_exact; % replacing with time of jump exact
end

tic
%% CSC
csc_filename= fullfile(exp_directory,'Neuralynx','CSC12.ncs');
for l=1:length(lap)
    if lap(l).dir~=direction || ~isCSC % filtering out undesired direction
        continue;
    end
    figure(100+l)
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase] = filtertheta(timecsc,lfp);
    
    idx = pos.lap==l;
    ax1 = subplot(4,1,1);
    plot(pos.t(idx)-lap(l).t_jump,pos.x(idx));
    ylabel('Horizontal position (cm)')
    title(['lap ' num2str(l) ' - ' direction 'ward']);
    
    ax2 = subplot(4,1,2); hold on;
    plot(pos.t(idx)-lap(l).t_jump,pos.filt.s(idx),'b');
    ylabel('Speed (cm/s)')
    ylim([0 300])
    
    ax3 = subplot(4,1,3); hold on
    plot(timecsc-lap(l).t_jump,phase)
    % theta phase
    for c=cluster_no
        idx = [cluster(c).lap]==l;
        if nnz(idx) > 0 % if there a firing for one of cluster_no in this lap
            plot(cluster(c).t(idx)-lap(l).t_jump, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        end
    end
    legend([{'Theta Phase'};legendCell])
    ylim([-180 180])
    ylabel('Phase (Degrees)')
    
    ax4 = subplot(4,1,4);
    plot(timecsc-lap(l).t_jump,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc-lap(l).t_jump,theta *1e6,'r');
    ylim([-400 400])
    ylabel('Theta (\muV)')
    xlabel('Time (sec)')
    
    linkaxes([ax1 ax2 ax3 ax4],'x')
    xlim([-1 1])
    
    set(gcf, 'Position', [100 100 1536 1175]);
    saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-' direction '-cl_' mat2str(cluster_no) '-lap_' num2str(l) '.svg']))
%     if length(cluster_no) == 1
%         
%     else
%         mat2str(cluster_no)
%         
%         saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-' direction 'ward-lap_' num2str(l) '.svg']))
%     end
    
end
toc
%% Phase Plots
flag1 = 0;
flag2 = 0;

for c = cluster_no
    
    % Phase-Postion (lap by lap)
    f = figure(1000); hold on
    ax = [];
    i = 0;
    for l=1:length(lap)
        if lap(l).dir~=direction % filtering out undesired direction
            continue;
        end
        i = i + 1;
        ax(i) = subplot(N,1,i);
        idx = [cluster(c).lap]==l & abs([cluster(c).vx]) > v_thresh;
        plot(cluster(c).x(idx),cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        hold on
        plot(cluster(c).x(idx), 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        title(['cluster(s) ' mat2str(cluster_no) ' lap ' num2str(l) ' ' direction 'ward'])
        ylabel('Phase (angle)')
        plot([640 640]/ppcm, [-180 540],'b','LineWidth',2);
        plot([892 892]/ppcm, [-180 540],'b','LineWidth',2);
        ylim([-180 540])
        zoom xon
    end
    xlabel('Horizontal position (cm)')
    f.WindowState = 'minimized';
    zoom xon
    linkaxes(ax,'x')
    xlim([0 xmax])
    
    % Phase-Postion (all laps)
    
    f = figure(1001); hold on
    flag = 0;
    if c == cluster_no(1) % only draw line for the fist time
        plot([640 640]/ppcm, [-180 540],'b','LineWidth',2);
        plot([892 892]/ppcm, [-180 540],'b','LineWidth',2);
    end
    for l=1:length(lap)
        if lap(l).dir~=direction % filtering out undesired direction
            continue;
        end
        idx = [cluster(c).lap]==l & abs([cluster(c).vx]) > v_thresh;
        plot(cluster(c).x(idx),cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        h = plot(cluster(c).x(idx), 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        plot_no = find(c==cluster_no); % for legend
        if nnz(idx) > 0 && flag == 0 % if there a firing for one of cluster_no in this lap
            h1(plot_no) = h;
            flag = 1; % just do this once
        end
    end
    title(['cluster(s) ' mat2str(cluster_no) ' all laps combined - direction: ' direction 'ward'])
    xlabel('Horizontal position (cm)')
    ylabel('Phase (angle)')
    ylim([-180 540])
    xlim([0 xmax])
    
    % Phase-Time (lap by lap)
    f = figure(1002); hold on
    i = 0;
    for l=1:length(lap)
        if lap(l).dir~=direction % filtering out undesired direction
            continue;
        end
        i = i + 1;
        ax(i) = subplot(N,1,i);
        idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_jump;
        plot(cluster(c).t(idx)-lap(l).t_jump, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        hold on
        plot(cluster(c).t(idx)-lap(l).t_jump, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        title(['cluster(s) ' mat2str(cluster_no) ' lap ' num2str(l)])
        ylabel('Phase (angle)')
        xlim([-1 1])
        ylim([-180 540])
        zoom xon
        hold on
        
    end
    f.WindowState = 'minimized';
    xlabel('Time (sec)')
    
    % Phase-Time (all laps)
    f = figure(1003); hold on;
    flag = 0;
    for l=1:length(lap)
        if lap(l).dir~=direction % filtering out undesired direction
            continue;
        end
        idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_jump;
        plot(cluster(c).t(idx)-lap(l).t_jump, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        h = plot(cluster(c).t(idx)-lap(l).t_jump, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        plot_no = find(c==cluster_no); % for legend
        if nnz(idx) > 0 && flag == 0 % if there a firing for one of cluster_no in this lap
            h2(plot_no) = h;
            flag = 1; % just do this once
        end
        title(['cluster ' num2str(c) ' lap ' num2str(l)])
    end
    title(['cluster(s) ' mat2str(cluster_no) ' all laps combined - direction: ' direction 'ward'])
    legend(legendCell)
    zoom xon
    xlim([-1 1])
    ylim([-180 540])
    xlabel('Time (sec)')
    ylabel('Phase (angle)')
    
end
f = figure(1001); hold on
set(gcf, 'Position', [100 100 1800 500]);
legend(h1, legendCell)
saveas(gcf,fullfile(exp_directory, 'Analysis',['Phase-Pos-cl' mat2str(cluster_no) '-' direction '.svg']))
f = figure(1003); hold on;
set(gcf, 'Position', [100 700 1800 500]);
xlim([-2 2])
legend(h2,legendCell)
saveas(gcf,fullfile(exp_directory, 'Analysis',['Phase-Time-cl' mat2str(cluster_no) '-' direction '.svg']))