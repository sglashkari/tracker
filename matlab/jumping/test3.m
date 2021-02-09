% SGL 2021-01-31
clc; 
addpath('../jumping');
exp_directory = '/home/shahin/Desktop/2020-11-22_Rat913-03';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax');
x_thresh = 84;
v_thresh = 20;

%lap_no = 10;
cluster_no = 6;
i = 1;
direction = 'left';
N = nnz([lap.dir]==direction);
isCSC = 1;

if isCSC
    close all;
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
    ax1 = subplot(5,1,1);
    plot(pos.t(idx),pos.x(idx));
    ylabel('Horizontal position (cm)')
    title(['lap ' num2str(l) ' - ' direction 'ward']);
    
    ax2 = subplot(5,1,2); hold on;
    plot(pos.t(idx),pos.filt.s(idx),'b');
    ylabel('Speed (cm/s)')
    ylim([0 300])
    
    ax3 = subplot(5,1,3);
    plot(timecsc,phase)
    ylim([-180 180])
    ylabel('Phase (Degrees)')
    
    ax4 = subplot(5,1,4);
    plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc,theta *1e6,'r');
    ylim([-350 350])
    ylabel('Theta (\muV)')
    
    % looking at all the clusters
    for c=1:length(cluster)
        if nnz([cluster(c).lap]==l) > 0 % if there a firing for cluster c in this lap
            ax5 = subplot(5,1,5);
            hold off; hold on
            idx = [cluster(c).lap]==l;
            plot(cluster(c).t(idx), c,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
            % theta phase
            ax3 = subplot(5,1,3); hold on;
            if nnz(cluster_no==c)
                plot(cluster(c).t(idx), cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
            end
        end
    end
    ax5 = subplot(5,1,5);
    xlabel('Time (sec)')
    ylim([0 length(cluster)+1])
    ylabel('Cluster Number')
    
    linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
    %xlim(lap(l).t);
    %zoom xon
    xlim([lap(l).t_jump-1 lap(l).t_jump+1])
    
    set(gcf, 'Position', [100 100 1536 1175]);
    saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-cl_ ' num2str(cluster_no) '-lap_' num2str(l) '.svg']))
    
end
toc
%% Phase-Postion (lap by lap)
f = figure(1000);
clf(f);

c = cluster_no;
for l=1:length(lap)
    if lap(l).dir~=direction % filtering out undesired direction
        continue;
    end
    ax(i) = subplot(N,1,i);
    idx = [cluster(c).lap]==l & abs([cluster(c).vx]) > v_thresh;
    plot(cluster(c).x(idx),cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    hold on
    plot(cluster(c).x(idx), 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    title(['cluster ' num2str(c) ' lap ' num2str(l) ' ' direction 'ward'])
    ylabel('Phase (angle)')
    zoom xon
    hold on
    i = i + 1;
end
xlabel('Horizontal position (cm)')
f.WindowState = 'minimized';

zoom xon
linkaxes(ax,'x')
xlim([0 xmax])

%% Phase-Postion (all laps)

f = figure(1001);
clf(f);
for l=1:length(lap)
    if lap(l).dir~=direction % filtering out undesired direction
        continue;
    end
    idx = [cluster(c).lap]==l & abs([cluster(c).vx]) > v_thresh;
    plot(cluster(c).x(idx),cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    hold on
    plot(cluster(c).x(idx), 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
end
title(['cluster ' num2str(c) ' all laps combined - direction: ' direction 'ward'])
xlabel('Horizontal position (cm)')
ylabel('Phase (angle)')
xlim([0 xmax])

%% Phase-Time (lap by lap)
f = figure(1002);
clf(f);
i = 1;
for l=1:length(lap)
    if lap(l).dir~=direction % filtering out undesired direction
        continue;
    end
    ax(i) = subplot(N,1,i);
    idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_jump;
    plot(cluster(c).t(idx)-lap(l).t_jump, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    hold on
    plot(cluster(c).t(idx)-lap(l).t_jump, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    title(['cluster ' num2str(c) ' lap ' num2str(l)])
    ylabel('Phase (angle)')
    xlim([-1 1])
    ylim([-180 540])
    zoom xon
    hold on
    i = i + 1;
end
f.WindowState = 'minimized';
xlabel('Time (sec)')

%% Phase-Time (all laps)
f = figure(1003);
clf(f);
for l=1:length(lap)
    if lap(l).dir~=direction % filtering out undesired direction
        continue;
    end
    idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_jump;
    plot(cluster(c).t(idx)-lap(l).t_jump, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    hold on
    plot(cluster(c).t(idx)-lap(l).t_jump, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    title(['cluster ' num2str(c) ' lap ' num2str(l)])
end
title(['cluster ' num2str(c) ' all laps combined'])
    zoom xon
    xlim([-1 1])
    ylim([-180 540])
xlabel('Time (sec)')
ylabel('Phase (angle)')

%{
v_thresh = 20;
f = figure(1001);
clf(f);
for l=2:2:30
    ax(round(l/2)) = subplot(15,1,round(l/2));
    idx = [cluster(c).lap]==l & [cluster(c).vx] <= v_thresh;
    x_offset = (cluster(c).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
    X = cluster(c).x(idx) + x_offset;
    phase = cluster(c).phase(idx);
    plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    hold on
    plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    %title(['cluster ' num2str(j) ' lap ' num2str(l)])
    ylabel('Phase (angle)')
    %xlabel('Horizontal position (cm)')
    zoom xon
    hold on
end
xlabel('Horizontal position (cm)')
zoom xon
linkaxes(ax,'x')
xlim([0 xmax])

%%
f = figure(1002);
clf(f);
for l=2:2:30
    idx = [cluster(c).lap]==l;
    x_offset = (cluster(c).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
    X = cluster(c).x(idx) + x_offset;
    phase = cluster(c).phase(idx);
    phase(cluster(c).vx(idx) <= v_thresh) = nan;
    plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    hold on
    plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    ylabel('Phase (angle)')
end
title(['cluster ' num2str(c)])
ylabel('Phase (angle)')
xlabel('Horizontal position (cm)')
xlim([0 xmax])

%%
f = figure(1003);
clf(f);
for l=16:2:30
    ax(round(l/2)-7) = subplot(8,1,round(l/2)-7);
    idx = [cluster(c).lap]==l;
    x_offset = (cluster(c).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
    X = cluster(c).x(idx) + x_offset;
    phase = cluster(c).phase(idx);
    phase(cluster(c).vx(idx) <= v_thresh) = nan;
    plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    hold on
    plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
    title(['cluster ' num2str(c) ' lap ' num2str(l)])
    ylabel('Phase (angle)')
end
xlabel('Horizontal position (cm)')
linkaxes(ax,'x')
xlim([0 xmax])
%}