
% close all
% ddd = randi(20,5);
% ddd = sort(ddd(:)/25);
% figure(1); histogram(ddd,'BinWidth',0.0333);
% [counts5,edges] = histcounts(ddd,'BinWidth',0.0333);
% hold on; plot(edges(2:end),counts5,'-or')
% counts5'
% edges'
% n = 10;
% edges50 = movmean(edges,n+1)'
% counts50 = movsum(counts5,n)'
%
% fl = floor(n/2);
% figure(2); bar(edges50(6:end-5),counts50(5:end-5))
% %bar(edges50(4:end-2),counts50(1:end-2))
% whos counts5 counts50 edges edges50
% ddd

close all
csc_filename= fullfile(exp_directory,'Neuralynx','CSC4.ncs');
addpath('../jumping');
lap_no = 10;
cluster_no = 6;

%% Directional rate map for all the laps
%figure(200)
%openfig(fullfile(exp_directory, 'Analysis','Directional_ratemap.fig'));


%% CSC

% theta phase
for j=1:length(cluster([cluster.m]==3))
    cluster(j).phase = nan(size(cluster(j).t));
end
for l=lap_no %1:length(lap)
    % filtering out the leftward laps
    if lap(l).dir=="right" %strcmp(lap(l).dir,"left")
        continue;
    end
    figure(100+l)
    [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
    [theta, phase] = filtertheta(timecsc,lfp);
    
    ax1 = subplot(5,1,1);
    plot(pos.t(pos.lap==l),pos.x(pos.lap==l));
    ylabel('Horizontal position (cm)')
    title(['lap ' num2str(l)]);
    
    ax2 = subplot(5,1,2);
    plot(pos.t(pos.lap==l),pos.s(pos.lap==l));
    ylabel('Speed (cm/s)')
    
    ax3 = subplot(5,1,3);
    plot(timecsc,phase)
    ylim([-180 180])
    ylabel('Phase (Degrees)')
    
    ax4 = subplot(5,1,4);
    plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc,theta *1e6,'r');
    ylim([-300 300])
    ylabel('Theta (\muV)')
    
    % looking at all the clusters in m3 (m3 is the whole recorded experiment)
    for j=1:length(cluster([cluster.m]==3))
        if nnz([cluster(j).lap]==l) > 0 % if there a firing for cluster j in this lap
            ax5 = subplot(5,1,5);
            hold off; hold on
            plot(cluster(j).t([cluster(j).lap]==l), cluster(j).no,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
            % theta phase
            cluster(j).phase([cluster(j).lap]==l) = interp1(timecsc,phase, cluster(j).t([cluster(j).lap]==l));
            cluster(j).theta([cluster(j).lap]==l) = interp1(timecsc,theta, cluster(j).t([cluster(j).lap]==l));
            if j == cluster_no
                ax3 = subplot(5,1,3); hold on;
                plot(cluster(j).t([cluster(j).lap]==l), cluster(j).phase([cluster(j).lap]==l),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
                ax4 = subplot(5,1,4); hold on;
                plot(cluster(j).t([cluster(j).lap]==l), cluster(j).theta([cluster(j).lap]==l)*1e6,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
            end
        end
    end
    ax5 = subplot(5,1,5);
    xlabel('Time (sec)')
    ylim([0 length(cluster([cluster.m]==3))+1])
    ylabel('Cluster Number')
    
    linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
    xlim(lap(l).t);
    
    set(gcf, 'Position', [100 100 1536 1175]);
    saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-lap_' num2str(l) '.jpg']))
    
end

%%
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'offset', 'colors','xmax','x_thresh');
x_thresh = 84;

j = cluster_no;

f = figure(1000);
clf(f);
for l=2:2:30
    ax(round(l/2)) = subplot(15,1,round(l/2));
    idx = [cluster(j).lap]==l;
    plot(cluster(j).x(idx),cluster(j).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
    hold on
    plot(cluster(j).x(idx), 360 + cluster(j).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
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

v_thresh = 20;
%%
f = figure(1001);
clf(f);
for l=2:2:30
    ax(round(l/2)) = subplot(15,1,round(l/2));
    idx = [cluster(j).lap]==l & [cluster(j).vx] <= v_thresh;
    x_offset = (cluster(j).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
    X = cluster(j).x(idx) + x_offset;
    phase = cluster(j).phase(idx);
    plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
    hold on
    plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
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
    idx = [cluster(j).lap]==l;
    x_offset = (cluster(j).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
    X = cluster(j).x(idx) + x_offset;
    phase = cluster(j).phase(idx);
    phase(cluster(j).vx(idx) <= v_thresh) = nan;
    plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
    hold on
    plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
    ylabel('Phase (angle)')
end
title(['cluster ' num2str(j)])
ylabel('Phase (angle)')
xlabel('Horizontal position (cm)')
xlim([0 xmax])

%%
f = figure(1003);
clf(f);
for l=16:2:30
    ax(round(l/2)-7) = subplot(8,1,round(l/2)-7);
    idx = [cluster(j).lap]==l;
    x_offset = (cluster(j).x(idx) < x_thresh) * (lap(l).gap - lap(1).gap);
    X = cluster(j).x(idx) + x_offset;
    phase = cluster(j).phase(idx);
    phase(cluster(j).vx(idx) <= v_thresh) = nan;
    plot(X,phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
    hold on
    plot(X, 360 + phase,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(j));
    title(['cluster ' num2str(j) ' lap ' num2str(l)])
    ylabel('Phase (angle)')
end
xlabel('Horizontal position (cm)')
linkaxes(ax,'x')
xlim([0 xmax])