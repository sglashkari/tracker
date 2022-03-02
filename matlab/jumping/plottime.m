%% Time plots lap by lap (jump and ditch)
% Phase plots lap by lap and combined
% SGL 2022-02-11 (originally test5_present)
clc; clear; close all
% [datafile,exp_directory] = uigetfile(fullfile('D:\Analysis', 'analyzed_data.mat'), 'Select Data File');
% if isequal(datafile, 0)
%     error('Data file was not selected!')
% end
exp_directory = 'D:\Analysis\2021-12-21';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'colors','xmax','hist','daq');

v_thresh = 10;
cluster_no = 40; %[17 18 23 28]; % [22 29]; [3 13 17 18 22 23 28 29 31]; [cluster.no];
direction = 'left';
showCSC = true;
timerange = [-2 2]; % 2 sec before to 2 sec after
x_reference = 'stationary';
t_reference = 'landing';

%%
legendCell = strings(size(cluster))';
for c=[cluster.no]
    if ismember(c,cluster_no)
        legendCell(c) = string(num2str(c, 'cluster #%-d'));
    else
        legendCell(c) = "";
    end
end

% time of interest: landing or take-off
for l = 1:length(lap)
    lap(l).t_land = lap(l).t_cross(end);
    if strcmp(x_reference,'stationary')
        lap(l).dx_interest = 0;
    elseif strcmp(x_reference,'moving')
        lap(l).dx_interest = lap(l).corr;
    end
    if strcmp(t_reference,'take-off')
        lap(l).t_interest = lap(l).t_jump_exact;
    elseif strcmp(t_reference,'landing')
        lap(l).t_interest = lap(l).t_land;
    end
end
no_states = length(unique([lap.status]));

%% CSC
tic
sh=cluster(cluster_no(1)).sh;
csc_filename=fullfile(exp_directory, ['CA1-Shank' num2str(sh)],'B','LFP32.ncs'); % use plotcsc to optimize it
[time,data] = read_bin_csc(csc_filename);
for l=[lap([lap.dir]==direction).no]  % only the desired direction
    if ~showCSC
        break;
    end
%     if lap(l).status ~= "jump"
%         continue;
%     end
    figure(100+l);
    legendCellModified = legendCell;
    
    idx = time >= lap(l).t_interest + timerange(1) & time <= lap(l).t_interest + timerange(2);
    timecsc = time(idx);
    lfp = data(idx);
    [theta, phase] = filterlfp(timecsc,lfp);
    
    idx = pos.t >= lap(l).t_interest + timerange(1) & pos.t <= lap(l).t_interest + timerange(2);
    ax1 = subplot(6,1,1); hold on
    plot(pos.t(idx)-lap(l).t_interest,pos.x(idx));
    ylabel('Horizontal position (cm)')
    title(['lap ' num2str(l) ', at time of ' t_reference ', ' direction 'ward ' convertStringsToChars(lap(l).status)]);
    
    
    ax2 = subplot(6,1,2); hold on;
    plot(pos.t(idx)-lap(l).t_interest,30*ones(size(pos.t(idx))),'k');
    idx = posi.t >= lap(l).t_interest + timerange(1) & posi.t <= lap(l).t_interest + timerange(2);
    plot(posi.t(idx)-lap(l).t_interest,abs(posi.filt.vx(idx)),'m');
    
    ylabel('Horizontal Speed (cm/s)')
    ylim([0 300])
    
    % plot a line for every time the rat cross the edges
    for i=1:length(lap(l).t_cross)
        plot(repmat(lap(l).t_cross(i),1,2)-lap(l).t_interest, ylim,lap(l).cross_color,'LineWidth',2);
    end
    
    idx = daq.t >= lap(l).t_interest + timerange(1) & daq.t <= lap(l).t_interest + timerange(2);
    ax3 = subplot(6,1,3); hold on
    plot(daq.t(idx)-lap(l).t_interest,daq.filt.loadcell(:,idx));
    ylabel('Force (N)')
    
    ax4 = subplot(6,1,4); hold on
    idx = posi.t >= lap(l).t_interest + timerange(1) & posi.t <= lap(l).t_interest + timerange(2);
    plot(posi.t(idx)-lap(l).t_interest, posi.filt.avy(idx));
    ylim([-600 600])
    ylabel('Pitch: Anglar velocity (deg/s)')
    
    ax5 = subplot(6,1,5); hold on
    plot(timecsc-lap(l).t_interest,phase,'Color','#ADD8E6')
    % theta phase
    for c=cluster_no
        %idx = [cluster(c).lap]==l;
        idx = cluster(c).t >= lap(l).t_interest + timerange(1) & cluster(c).t <= lap(l).t_interest + timerange(2);
        if nnz(idx) > 0 % if there a firing for one of cluster_no in this lap
            plot(cluster(c).t(idx)-lap(l).t_interest, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        else
            legendCellModified(c)="";
        end
    end
    legend([{'Theta Phase'};legendCellModified(legendCellModified~="")])
    ylim([-180 180])
    ylabel('Phase (Degrees)')
    
    ax6 = subplot(6,1,6);
    plot(timecsc-lap(l).t_interest,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc-lap(l).t_interest,theta *1e6,'r');
    ylim([-400 400])
    ylabel('Theta (\muV)')
    xlabel('Time (sec)')
    
    set(gcf, 'Position', [100 100 1536 1175]);
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')
    xlim(timerange)
    
    saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-' direction '-cl_' mat2str(cluster_no) '-lap_' num2str(l) '.png']))
end
figure(100+l+1);
toc
%% Phase Plots
N = nnz([lap.dir]==direction);
legendCellModified11 = legendCell;
legendCellModified12 = legendCell;
legendCellModified21 = legendCell;
legendCellModified22 = legendCell;
h11 = [];
h12 = [];
h21 = [];
h22 = [];
for c = cluster_no
    
    % Phase-Postion (lap by lap)
    f = figure(1000); hold on
    ax = [];
    i = 0;
    for l=[lap([lap.dir]==direction).no]  % only the desired direction
        i = i + 1;
        ax(i) = subplot(N,1,i);
        idx = [cluster(c).lap]==l & abs([cluster(c).vx]) >= v_thresh;
        plot(cluster(c).x(idx) + lap(l).dx_interest ,cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        hold on
        ylim([-180 180])
        set(ax(i),'yticklabel',[])
        if i == 1
            title(['Cluster(s) ' mat2str(cluster_no) ' (lap by lap) at ' x_reference ' frame of ref, ' direction 'ward']);
        end
        if c == cluster_no(1) && mod(i,5)==0
            ylabel(['Lap #' num2str(l)]);
        end        
        % gap range
        plot(repmat(lap(l).gap(1) + lap(l).dx_interest,2,1), ylim,lap(l).cross_color,'LineWidth',2);
        plot(repmat(lap(l).gap(2) + lap(l).dx_interest,2,1), ylim,lap(l).cross_color,'LineWidth',2);
        xlim([0 xmax])
    end
    xlabel('Horizontal position (cm)')
    zoom xon
    linkaxes(ax,'x')
    set(ax,'xticklabel',[])
    set(ax(i),'xticklabel',0:50:xmax)
    f.WindowState = 'maximized';
    
    % Phase-Postion (all laps)
    
    f = figure(1001); hold on
    i = 0;
    ax = [];
    for state = unique([lap.status])
        flag = 0;
        i = i + 1;
        ax(i) = subplot(length(unique([lap.status])),1,i);
        for l=[lap([lap.dir]==direction & [lap.status]==state).no]  % only the desired direction and state
            
            idx = [cluster(c).lap]==l & abs([cluster(c).vx]) >= v_thresh;
            plot(cluster(c).x(idx) + lap(l).dx_interest,cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on;
            ylim([-180 540])
            h = plot(cluster(c).x(idx) + lap(l).dx_interest, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
            plot_no = find(c==cluster_no); % for legend
            if nnz(idx) > 0 && flag == 0 % if there a firing for one of cluster_no in this lap
                h1(plot_no) = h;
                flag = 1; % just do this once
            end
            % gap range
            if c == cluster_no(1) % only draw line for the fist time
                plot(repmat(lap(l).gap(1) + lap(l).dx_interest,2,1), ylim,lap(l).cross_color,'LineWidth',0.15);
                plot(repmat(lap(l).gap(2) + lap(l).dx_interest,2,1), ylim,lap(l).cross_color,'LineWidth',0.15);
            end
        end
        title(['cluster(s) ' mat2str(cluster_no) ' all laps at ' x_reference ' frame of ref, ' direction 'ward ' convertStringsToChars(state) ' (velocity filtered: v >= ' num2str(v_thresh) ' cm/s)'])
        ylabel('Phase (angle)')
    end
    
    xlabel('Horizontal position (cm)')
    linkaxes(ax,'x')
    zoom xon
    xlim([0 xmax])
    
    % Phase-Time (lap by lap)
    f = figure(1002); hold on
    i = 0;
    ax=[];
    for l=[lap([lap.dir]==direction).no]  % only the desired direction
        i = i + 1;
        ax(i) = subplot(N,1,i);
        idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_interest; % needs to be updated
        plot(cluster(c).t(idx)-lap(l).t_interest, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        hold on
        if i == 1
            title(['Spike train for cluster(s) ' mat2str(cluster_no) ' (lap by lap) aligned at time of ' t_reference ', ' direction 'ward']);
        end
        if c == cluster_no(1) && mod(i,5)==0
            ylabel(['Lap #' num2str(l)]);
        end
        xlim(timerange)
        ylim([-180 180])
        set(ax(i),'yticklabel',[])
        % cross time
        for j=1:length(lap(l).t_cross)
            plot(repmat(lap(l).t_cross(j),1,2)-lap(l).t_interest, ylim,lap(l).cross_color,'LineWidth',2);
        end
        zoom xon
        hold on
    end
    
    linkaxes(ax,'x')
    zoom xon
    xlim(timerange)
    xlabel('Time (sec)')
    %set(ax,'xticklabel',[])
    %set(ax(i),'xticklabel',0:50:xmax)
    f.WindowState = 'maximized';
    
    % Phase-Time (all laps)
    f = figure(1003); hold on;
    i = 0;
    ax = [];
    for state = unique([lap.status])
        flag = 0;
        i = i + 1;
        h = [];
        ax(i) = subplot(length(unique([lap.status])),1,i);
        for l=[lap([lap.dir]==direction & [lap.status]==state).no]  % only the desired direction and state
            
            % cross time
            if c == cluster_no(1) % only draw line for the fist time
                for j=1:length(lap(l).t_cross)
                    plot(repmat(lap(l).t_cross(j),1,2)-lap(l).t_interest, [-180 540],lap(l).cross_color,'LineWidth',0.15);
                    ylim([-180 540])
                end
            end
            
            %idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_interest;
            idx = cluster(c).t >= lap(l).t_interest + timerange(1) & cluster(c).t <= lap(l).t_interest + timerange(2);
            if nnz(idx) > 0 % if there a firing for one of cluster_no in this lap
                plot(cluster(c).t(idx)-lap(l).t_interest, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
                h = plot(cluster(c).t(idx)-lap(l).t_interest, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
                flag = 1;
            end
        end
        % legend (skip if no cluster exists in the timerange)
        if flag == 0
            if i == 1
                legendCellModified21(c)="";
            else
                legendCellModified22(c)="";
            end
        end
        if i == 1
            h21 = [h21 h];
        else
            h22 = [h22 h];
        end
        ylabel('Phase (angle)')
        title(['cluster(s) ' mat2str(cluster_no) ' all laps aligned at time of ' t_reference ' - ' direction 'ward ' convertStringsToChars(state)]);
        %legend(legendCell)
    end
    zoom xon
    linkaxes(ax,'x')
    xlim(timerange)
    xlabel('Time (sec)')
end

figure(1000); hold on
saveas(gcf,fullfile(exp_directory, 'Analysis',['LapByLap-Pos-cl' mat2str(cluster_no) '-' direction '.png']))

figure(1001); hold on
set(gcf, 'Position', [100 100 1800 500]);
legend(h1, legendCell(legendCell~=""))
saveas(gcf,fullfile(exp_directory, 'Analysis',['Phase-Pos-cl' mat2str(cluster_no) '-' direction '.png']))

figure(1002); hold on
saveas(gcf,fullfile(exp_directory, 'Analysis',['LapByLap-Time-cl' mat2str(cluster_no) '-' direction '.png']))

figure(1003); hold on;
set(gcf, 'Position', [100 100 1800 500]);
xlim(timerange)
% legend
subplot(2,1,1);
legend(h21,legendCellModified21(legendCellModified21~=""))
subplot(2,1,2);
legend(h22,legendCellModified22(legendCellModified22~=""))

saveas(gcf,fullfile(exp_directory, 'Analysis',['Phase-Time-cl' mat2str(cluster_no) '-' direction '.png']))

%% Rat-map with velocity filter (for stationary or moving frame)
N = length(cluster_no);
figure(2000); clf;
i = 0;
for c = cluster_no
    
    for state = unique([lap.status])
        for dir = ["left" "right"]
            i = i+1;
            
            idx = posi.status==state & posi.dir==dir & abs(posi.vx) >= v_thresh;
            hist.posi = histcounts(posi.x(idx) + lap(l).dx_interest, hist.edges) * posi.dt; % seconds in each bin
            idx = [cluster(c).status]==state & [cluster(c).dir]==dir & abs([cluster(c).vx]) >= v_thresh;
            hist.cluster = histcounts(cluster(c).x(idx) + lap(l).dx_interest, hist.edges); % spikes in each bin
            hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
            
            a(i) = subplot(N*no_states,2,i);
            histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c),'FaceAlpha',0.8); hold on; %rate map histogram (20% transparent)
            ylabel(['#' num2str(c) ':' cluster(c).region]);
            if i <= no_states
                title(['Rate Map (' convertStringsToChars(dir) 'ward)'])
            end
            if i >= 2*N*no_states-1
                xlabel('Horizontal position (cm)')
            end
            xlim(max([lap.dx_interest])+[0 xmax])
            if mod(i,2*no_states)==0
                linkaxes([a(i) a(i-1) a(i-2) a(i-3)],'y')
                ylim(max([0 5],ylim)); % start from 0 when no cell rate map is null
                states = repmat(unique([lap.status]),2,1); % alternating between ditch and jump for gaps
                states = states(end:-1:1); % correcting the order
                directions = ["right" "left" "right" "left"]; % alternating between left and right
                for j=i:-1:i-3
                    subplot(N*no_states,2,j);
                    for l=[lap([lap.dir]==directions(i-j+1) & [lap.status]==states(i-j+1)).no]  % only the desired direction
                        % gap range
                        plot(repmat(lap(l).gap(1) + lap(l).dx_interest,1,2), ylim,lap(l).cross_color,'LineWidth',0.25);
                        plot(repmat(lap(l).gap(2) + lap(l).dx_interest,1,2), ylim,lap(l).cross_color,'LineWidth',0.25);
                    end
                    %Reverse the stacking order so that the patch overlays the line
                    chi=get(gca, 'Children');
                    set(gca, 'Children',flipud(chi));
                end
            end
            
        end
    end
end

toc
set(gcf, 'Position', [100 0 1000 200+100*N*no_states]);
linkaxes(a,'x')
zoom xon
xlim(max([lap.dx_interest])+[0 xmax])
sgtitle(['Directional ratemap, ' x_reference ' frame of ref (velocity filtered: v >= ' num2str(v_thresh) ' cm/s)']);
saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_ratemap-cl' mat2str(cluster_no) '_' x_reference '.png']))