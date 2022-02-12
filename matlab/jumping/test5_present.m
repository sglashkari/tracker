% SGL 2021-01-31 (test3_present + rate map)
clc; clear; close all
exp_directory = 'D:\Analysis\2021-12-21';

mat_filename = fullfile(exp_directory,'data.mat');
load(mat_filename,'pos','daq');
pos.av = angvel(quaternion(pos.q),1/300,'frame');

mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'colors','xmax','hist');
%x_thresh = 138;
v_thresh = 0;

cluster_no = 2;
direction = 'right';

legendCell = cellstr(num2str(cluster_no', 'cluster #%-d'));
N = nnz([lap.dir]==direction);
isCSC = true;

%% CSC
tic
sh=cluster(cluster_no(1)).sh;
csc_filename=fullfile(exp_directory, ['CA1-Shank' num2str(sh)],'B','LFP32.ncs'); % use plotcsc to optimize it
[time,data] = read_bin_csc(csc_filename);

for l=[lap([lap.dir]==direction).no]  % only the desired direction
    if ~isCSC
        break;
    end
    figure(100+l);
    
    idx = time >= lap(l).t(1) & time <= lap(l).t(2);
    timecsc = time(idx);
    lfp = data(idx);
    [theta, phase] = filtertheta(timecsc,lfp);
    
    idx = pos.lap == l;
    ax1 = subplot(6,1,1); hold on
    plot(pos.t(idx)-lap(l).t_jump_exact,pos.x(idx));
    ylabel('Horizontal position (cm)')
    title(['lap ' num2str(l) ' - ' direction 'ward ' convertStringsToChars(lap(l).status)]);

    
    ax2 = subplot(6,1,2); hold on;
    plot(posi.t(posi.lap==l)-lap(l).t_jump_exact,abs(posi.filt.vx(posi.lap==l)),'m');
    plot(pos.t(idx)-lap(l).t_jump_exact,30*ones(size(pos.t(idx))),'k');
    ylabel('Horizontal Speed (cm/s)')
    ylim([0 300])

    % plot a line for every time the rat cross the edges
    for i=1:length(lap(l).t_cross)
        plot(repmat(lap(l).t_cross(i),1,2)-lap(l).t_jump_exact, ylim,'b','LineWidth',2);
    end
    
    idx = daq.t >= lap(l).t(1) & daq.t <= lap(l).t(2);
    ax3 = subplot(6,1,3); hold on
    plot(daq.t(idx)-lap(l).t_jump_exact,daq.filt.loadcell(:,idx));
    ylabel('Force (N)')
    
    ax4 = subplot(6,1,4); hold on
    %plot(posi.t(posi.lap==l)-lap(l).t_jump_exact, posi.filt.avx(posi.lap==l));
    plot(posi.t(posi.lap==l)-lap(l).t_jump_exact, posi.filt.avy(posi.lap==l));
    %plot(pos.t(pos.lap==l)-lap(l).t_jump_exact, pos.pitch(pos.lap==l));
    ylim([-600 600])
    ylabel('Pitch: Anglar velocity (deg/s)')
    
    ax5 = subplot(6,1,5); hold on
    plot(timecsc-lap(l).t_jump_exact,phase,'Color','#ADD8E6')
    % theta phase
    for c=cluster_no
        idx = [cluster(c).lap]==l;
        if nnz(idx) > 0 % if there a firing for one of cluster_no in this lap
            plot(cluster(c).t(idx)-lap(l).t_jump_exact, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        end
    end
    legend([{'Theta Phase'};legendCell])
    ylim([-180 180])
    ylabel('Phase (Degrees)')
    
    ax6 = subplot(6,1,6);
    plot(timecsc-lap(l).t_jump_exact,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc-lap(l).t_jump_exact,theta *1e6,'r');
    ylim([-400 400])
    ylabel('Theta (\muV)')
    xlabel('Time (sec)')
    
    set(gcf, 'Position', [100 100 1536 1175]);
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')
    xlim([-2 2])

    saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-' direction '-cl_' mat2str(cluster_no) '-lap_' num2str(l) '.png']))
end
figure(100+l+1);
toc
%% Phase Plots
flag1 = 0;
flag2 = 0;

for c = cluster_no
    
    % Phase-Postion (lap by lap)
    f = figure(1000); hold on
    ax = [];
    i = 0;
    for l=[lap([lap.dir]==direction).no]  % only the desired direction
        i = i + 1;
        ax(i) = subplot(N,1,i);
        idx = [cluster(c).lap]==l & abs([cluster(c).vx]) > v_thresh;
        plot(cluster(c).x(idx),cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        %hold on
        %plot(cluster(c).x(idx), 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        title(['cluster(s) ' mat2str(cluster_no) ' lap ' num2str(l) ' ' direction 'ward'])
        ylabel('Phase (angle)')
        plot(repmat(lap(l).gap(1),2,1), [-180 540],'b','LineWidth',2);
        plot(repmat(lap(l).gap(2),2,1), [-180 540],'b','LineWidth',2);
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
        plot(repmat(lap(l).gap(1),2,1), [-180 540],'b','LineWidth',2);
        plot(repmat(lap(l).gap(2),2,1), [-180 540],'b','LineWidth',2);
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
    title(['cluster(s) ' mat2str(cluster_no) ' all laps combined - direction: ' direction 'ward  (velocity filtered: v > ' num2str(v_thresh) ' cm/s)'])
%     if strcmp(direction,'left')
%         set ( gca, 'xdir', 'reverse' )
%     end
    xlabel('Horizontal position (cm)')
    ylabel('Phase (angle)')
    ylim([-180 540])
    xlim([0 xmax])
    
    % Phase-Time (lap by lap)
    f = figure(1002); hold on
    i = 0;
    for l=[lap([lap.dir]==direction).no]  % only the desired direction
        i = i + 1;
        ax(i) = subplot(N,1,i);
        idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_jump_exact;
        plot(cluster(c).t(idx)-lap(l).t_jump_exact, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        hold on
        plot(cluster(c).t(idx)-lap(l).t_jump_exact, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
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
        idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_jump_exact;
        plot(cluster(c).t(idx)-lap(l).t_jump_exact, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        h = plot(cluster(c).t(idx)-lap(l).t_jump_exact, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        plot_no = find(c==cluster_no); % for legend
        if nnz(idx) > 0 && flag == 0 % if there a firing for one of cluster_no in this lap
            h2(plot_no) = h;
            flag = 1; % just do this once
        end
        title(['cluster ' num2str(c) ' lap ' num2str(l)])
        
        % range of gap:
%         idx = pos.lap==l;
%         gapidx = pos.x(idx)>lap(l).gap(1) & pos.x(idx)<=lap(l).gap(2);
%         gapidx = [0; diff(gapidx)];
%         post = pos.t(idx);
        for i=1:length(lap(l).t_cross)
            plot(repmat(lap(l).t_cross(i),1,2)-lap(l).t_jump_exact, [-180 540],'b','LineWidth',0.25);
        end
    
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
saveas(gcf,fullfile(exp_directory, 'Analysis',['Phase-Pos-cl' mat2str(cluster_no) '-' direction '.png']))
f = figure(1003); hold on;
set(gcf, 'Position', [100 700 1800 500]);
xlim([-2 2])
legend(h2,legendCell)
saveas(gcf,fullfile(exp_directory, 'Analysis',['Phase-Time-cl' mat2str(cluster_no) '-' direction '.png']))

%}
%% Rat-map with velocity filter
N = length(cluster_no);
i = 0;
for c = cluster_no
    
    figure(2000)
    i = i+1;
    % leftward rate map
    idx = posi.dir=="left" & abs(posi.vx) > v_thresh;
    hist.posi = histcounts(posi.x(idx), hist.edges) * posi.dt; % seconds in each bin
    idx = [cluster(c).dir]=="left" & abs([cluster(c).vx]) > v_thresh;
    hist.cluster = histcounts(cluster(c).x(idx), hist.edges); % spikes in each bin
    hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
    
    a(2*i-1) = subplot(N,2,2*i-1);
    %a(i) = subplot(N,1,i);
    histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
    ylabel(['cluster ' num2str(c)]);
    xlabel('Horizontal position (cm)');
    if i == 1
        title('Rate Map (leftward)')
    end
    if i == N
        xlabel('Horizontal position (cm)')
    end
    
%     set ( gca, 'xdir', 'reverse' )
    
    
    % rightward rate map
    idx = posi.dir=="right" & abs(posi.vx) > v_thresh;
    hist.posi = histcounts(posi.x(idx), hist.edges) * posi.dt; % seconds in each bin
    idx = [cluster(c).dir]=="right" & abs([cluster(c).vx]) > v_thresh;
    hist.cluster = histcounts(cluster(c).x(idx), hist.edges); % spikes in each bin
    hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
    
    a(2*i) = subplot(N,2,2*i);
    histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
    ylabel(['cluster ' num2str(c)]);
    linkaxes([a(2*i-1) a(2*i)],'y')
    ymax = ylim;
    ylim([0 ymax(2)]) % start from 0 when no cell rate map is null
    
    if i == 1
        title('Rate Map (rightward)')
    end
    if i == N
        xlabel('Horizontal position (cm)')
    end
    
    % edge of gap
    hold on
    plot(repmat(lap(l).gap(1),2,1), ylim,'b','LineWidth',2);
    plot(repmat(lap(l).gap(2),2,1), ylim,'b','LineWidth',2);
%     set ( gca, 'xdir', 'reverse' )
    a(2*i-1) = subplot(N,2,2*i-1);
    hold on
    plot(repmat(lap(l).gap(1),2,1), ylim,'b','LineWidth',2);
    plot(repmat(lap(l).gap(2),2,1), ylim,'b','LineWidth',2);
end

toc
set(gcf, 'Position', [100 100 1000 200+100*N]);
linkaxes(a,'x')
zoom xon
xlim([0 xmax])
sgtitle(['Directional Ratemap (velocity filtered: v > ' num2str(v_thresh) ' cm/s)']);
saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_ratemap-cl' mat2str(cluster_no) '.png']))