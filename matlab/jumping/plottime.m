%% Time plots lap by lap (jump and ditch)
% Phase plots lap by lap and combined
% SGL 2022-02-11 (originally test5_present)

cluster_no
s_thresh
direction
showCSC
timerange
x_reference
t_reference

%posi.s = vecnorm([posi.vx; posi.vy]);
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
    else
        error("t_reference should be either take-off or landing");
    end
end

depth = lap(l).gap_depth;
no_states = length(unique([lap.status]));

%% CSC
%tic
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
    title(['lap ' num2str(l) ', at the time of ' t_reference ', ' direction 'ward ' convertStringsToChars(lap(l).status)]);
    
    
    ax2 = subplot(6,1,2); hold on;
    %plot(pos.t(idx)-lap(l).t_interest,30*ones(size(pos.t(idx))),'k');
    plot(timerange, [30 30],'k','LineWidth',0.1);
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
    idx = pos.t >= lap(l).t_interest + timerange(1) & pos.t <= lap(l).t_interest + timerange(2);
    plot(pos.t(idx)-lap(l).t_interest, pos.p(idx,3),'.');
    ylim([-30 15])
    ylabel('Elevation (cm)')
    
    ax5 = subplot(6,1,5); hold on
    plot(timecsc-lap(l).t_interest,phase,'Color','#ADD8E6')
    % theta phase
    for c=cluster_no
        %idx = [cluster(c).lap]==l;
        idx = cluster(c).t >= lap(l).t_interest + timerange(1) & cluster(c).t <= lap(l).t_interest + timerange(2);
        if nnz(idx) > 0 % if there a firing for one of cluster_no in this lap
            plot(cluster(c).t(idx)-lap(l).t_interest, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        else
            plot(0, -360,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
            %             legendCellModified(c)="";
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
    
    set(gcf, 'Position', [0 0 1600 1300]);
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')
    xlim(timerange)
    
    saveas(gcf,fullfile(exp_directory, 'Analysis',['CSC-' direction '-cl_' mat2str(cluster_no) '-lap_' num2str(l) '.png']))
end
%figure(100+l+1);
%toc
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
    for dir = ["left" "right"]
        N = nnz([lap.dir]==dir);
        i = 0;
        for l=[lap([lap.dir]==dir).no]  % only the desired direction
            i = i + 1;
            if dir == "left"
                ax(i) = subplot(N,2,2*i-1);
            else
                ax(i) = subplot(N,2,2*i);
            end
            idx = [cluster(c).lap]==l & abs([cluster(c).s]) >= s_thresh;
            plot(cluster(c).x(idx) + lap(l).dx_interest ,cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
            hold on
            ylim([-180 180])
            set(ax(i),'yticklabel',[])
            if i == 1
                title(['Cluster(s) ' mat2str(cluster_no) ' (lap by lap) at ' x_reference ' frame of ref, ' convertStringsToChars(dir) 'ward']);
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
    end
    f.WindowState = 'maximized';
    
    % Phase-Postion (all laps)
    
    f = figure(1001); hold on
    i = 0;
    ax = [];
    for state = unique([lap.status])
        for dir = ["left" "right"]
            flag = 0;
            i = i + 1;
            ax(i) = subplot(length(unique([lap.status])),2,i);
            for l=[lap([lap.dir]==dir & [lap.status]==state).no]  % only the desired direction and state
                
                idx = [cluster(c).lap]==l & abs([cluster(c).s]) >= s_thresh;
                %plot(cluster(c).x(idx) + lap(l).dx_interest,cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on;
                h = scatter(cluster(c).x(idx) + lap(l).dx_interest, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on;
                alpha(h,opacity)
                ylim([-180 540])
                %h = plot(cluster(c).x(idx) + lap(l).dx_interest, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
                h = scatter(cluster(c).x(idx) + lap(l).dx_interest, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
                alpha(h,opacity)
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
            title(['cluster(s) ' mat2str(cluster_no) ' all laps at ' x_reference ' frame of ref, ' convertStringsToChars(dir) 'ward ' convertStringsToChars(state) ' (speed filtered: s >= ' num2str(s_thresh) ' cm/s)'])
            ylabel('Phase (angle)')
            
            if i>length(unique([lap.status]))*2-2
                xlabel('Horizontal position (cm)')
            end
        end
    end
    
    
    linkaxes(ax,'x')
    zoom xon
    if isZoom
       xlim([-50+min([lap.gap]) 50+max([lap.gap])]);
    else
       xlim([0 xmax]);
    end
    
    % Phase-Time (lap by lap)
    f = figure(1002); hold on
    for dir = ["left" "right"]
        N = nnz([lap.dir]==dir)+1;
        i = 0;
        ax=[];
        
        for l=[lap([lap.dir]==dir).no]  % only the desired direction
            i = i + 1;
            if dir == "left"
                ax(i) = subplot(N,2,2*i-1);
            else
                ax(i) = subplot(N,2,2*i);
            end
            idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_interest; % needs to be updated
            plot(cluster(c).t(idx)-lap(l).t_interest, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
            hold on
            if i == 1
                title(['Spike train for cluster(s) ' mat2str(cluster_no) ' (lap by lap) aligned at the time of ' t_reference ', ' convertStringsToChars(dir) 'ward']);
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
    end
    f.WindowState = 'maximized';
    
    % Phase-Time (all laps)
    f = figure(1003); hold on;
    i = 0;
    ax = [];
    for state = unique([lap.status])
        for dir = ["left" "right"]
            flag = 0;
            i = i + 1;
            h = [];
            ax(i) = subplot(length(unique([lap.status])),2,i);
            for l=[lap([lap.dir]==dir & [lap.status]==state).no]  % only the desired direction and state
                
                % cross time
                if c == cluster_no(1) % only draw line for the fist time
                    for j=1:length(lap(l).t_cross)
                        plot(repmat(lap(l).t_cross(j),1,2)-lap(l).t_interest, [-180 540],lap(l).cross_color,'LineWidth',0.05);
                        ylim([-180 540])
                    end
                end
                
                %idx = [cluster(c).lap]==l;% & [cluster(c).t] < lap(l).t_interest;
                idx = cluster(c).t >= lap(l).t_interest + timerange(1) & cluster(c).t <= lap(l).t_interest + timerange(2);
                if nnz(idx) > 0 % if there a firing for one of cluster_no in this lap
                    %plot(cluster(c).t(idx)-lap(l).t_interest, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
                    h = scatter(cluster(c).t(idx)-lap(l).t_interest, cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
                    alpha(h,opacity)
                    %h = plot(cluster(c).t(idx)-lap(l).t_interest, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
                    h = scatter(cluster(c).t(idx)-lap(l).t_interest, 360 + cluster(c).phase(idx),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
                    alpha(h,opacity)
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
            title(['cluster(s) ' mat2str(cluster_no) ' all laps aligned at the time of ' t_reference ' - ' convertStringsToChars(dir) 'ward ' convertStringsToChars(state)]);
            %legend(legendCell)
            if i>length(unique([lap.status]))*2-2
                xlabel('Time (sec)')
            end
        end
    end
    zoom xon
    linkaxes(ax,'x')
    xlim(timerange)
end

figure(1000); hold on
saveas(gcf,fullfile(exp_directory, 'Analysis',['LapByLap-Pos-cl' mat2str(cluster_no) '.png']))

figure(1001); hold on
set(gcf, 'Position', [0 100 1800 700]);
legend(h1, legendCell(legendCell~=""))
saveas(gcf,fullfile(exp_directory, 'Analysis',['Phase-Pos-cl' mat2str(cluster_no) '.png']))

figure(1002); hold on
saveas(gcf,fullfile(exp_directory, 'Analysis',['LapByLap-Time-cl' mat2str(cluster_no) '.png']))

figure(1003); hold on;
set(gcf, 'Position', [0 0 1800 700]);
xlim(timerange)
% legend
subplot(2,2,1);
legend(h21,legendCellModified21(legendCellModified21~=""))
subplot(2,2,2);
legend(h22,legendCellModified22(legendCellModified22~=""))

saveas(gcf,fullfile(exp_directory, 'Analysis',['Phase-Time-cl' mat2str(cluster_no) '.png']))

%% Rat-map with velocity filter (for stationary or moving frame)
N = length(cluster_no);
M = min([10 N]); % groups of 10 max
figure(2001); clf;
clear a
i = 0;
j = 1;
for c = cluster_no
    
    for state = unique([lap.status])
        for dir = ["left" "right"]
            i = i+1;
            idx = posi.status==state & posi.dir==dir & abs(posi.s) >= s_thresh & posi.lap > 0;
            dx = [lap(posi.lap(idx)).dx_interest];
            dx = reshape(dx,1,[]);
            hist.posi = histcounts(posi.x(idx) + dx, hist.edges) * posi.dt; % seconds in each bin
            idx = [cluster(c).status]==state & [cluster(c).dir]==dir & abs([cluster(c).s]) >= s_thresh & cluster(c).lap > 0;
            dx = [lap(cluster(c).lap(idx)).dx_interest];
            dx = reshape(dx,[],1);
            hist.cluster = histcounts(cluster(c).x(idx) + dx, hist.edges); % spikes in each bin
            hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
            a(i) = subplot(M*no_states,2,i);
            histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c),'FaceAlpha',0.8); hold on; %rate map histogram (20% transparent)
            if N == M % for small sets
                ylabel('Spike count')
                title(['#' num2str(c) ':' cluster(c).region ', ' convertStringsToChars(dir) 'ward ' convertStringsToChars(state)])
            else % for large sets
                set(gca,'XTick', []);
                ylabel(['#' num2str(c)])
            end
            
            if isZoom
                xlim([-50+min([lap.gap]) 50+max([lap.gap])]);
            else
                %xlim(max([lap.dx_interest])+[0 xmax]);
                xlim([0 xmax]);
            end
            
            if mod(i,2*no_states)==0
                linkaxes([a(i) a(i-1) a(i-2) a(i-3)],'y')
                ylim(max([0 2],ylim)); % start from 0 when no cell rate map is null
                states = repmat(unique([lap.status]),2,1); % alternating between ditch and jump for gaps
                states = states(end:-1:1); % correcting the order
                directions = ["right" "left" "right" "left"]; % alternating between left and right
                for k=i:-1:i-3
                    subplot(M*no_states,2,k);
                    for l=[lap([lap.dir]==directions(i-k+1) & [lap.status]==states(i-k+1)).no]  % only the desired direction and state
                        % gap range
                        plot(repmat(lap(l).gap(1) + lap(l).dx_interest,1,2), ylim,lap(l).cross_color,'LineWidth',0.25);
                        plot(repmat(lap(l).gap(2) + lap(l).dx_interest,1,2), ylim,lap(l).cross_color,'LineWidth',0.25);
                    end
                    %Reverse the stacking order so that the patch overlays the line
                    chi=get(gca, 'Children');
                    set(gca, 'Children',flipud(chi));
                end
            end
            
            if i == 2*M*no_states || (c == cluster_no(end) && i+(j-1)*2*M*no_states == 2*N*no_states)
                subplot(M*no_states,2,i);
                xlabel('Horizontal position (cm)')
                subplot(M*no_states,2,i-1);
                xlabel('Horizontal position (cm)')
                set(gcf, 'Position', [0 0 1800 300+200*N*no_states]);
                sgtitle(['Directional ratemap ' num2str(j) ', ' x_reference ' frame of ref (speed filtered: s >= ' num2str(s_thresh) ' cm/s)']);
                saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_ratemap_' num2str(j) '-cl' mat2str(cluster_no) '_' x_reference '.png']));
                i = 0;
                j = j + 1;
                if c < cluster_no(end)
                    figure(2000+j); clf;
                end
                a = [];
            end
            
        end
    end
end

%% time histogram (for take-off or landing)

% binning
hist.time_bin_size = 0.05; % sec
hist.time_edges = timerange(1):hist.time_bin_size:timerange(2); % size of image to cm

N = length(cluster_no);
figure(3000); clf;
i = 0;
for c = cluster_no
    
    for state = unique([lap.status])
        for dir = ["left" "right"]
            i = i+1;
            
            cluster(c).t_rel = []; % time relative to t_reference
            for l=[lap([lap.dir]==dir & [lap.status]==state).no]  % only the desired direction and state
                idx = cluster(c).t >= lap(l).t_interest + timerange(1) & cluster(c).t <= lap(l).t_interest + timerange(2);
                cluster(c).t_rel = [cluster(c).t_rel; cluster(c).t(idx) - lap(l).t_interest];
            end
            n = length([lap([lap.dir]==dir & [lap.status]==state).no]); % number of laps
            
            hist.cluster = histcounts(cluster(c).t_rel, hist.time_edges); % spikes in each bin
            hist.time_hist = hist.cluster / (hist.time_bin_size*n) ; % spikes per bin per lap
            
            a(i) = subplot(N*no_states,2,i);
            h = histogram('BinCounts', hist.time_hist, 'BinEdges', hist.time_edges, 'FaceColor',colors(c),'FaceAlpha',0.8); hold on; %rate map histogram (20% transparent)
            hold on
            y = findobj(h,'type','histogram');
            [maxVal,argmaxVal] = max(y.Values);
            t_max = (y.BinEdges(argmaxVal)+y.BinEdges(argmaxVal+1))/2;
            %pd = fitdist(cluster(c).t_rel,'Poisson');
            %hold on; h = histfit(cluster(c).t_rel,y.NumBins,'normal');
            
            %plot(hist.time_edges,pdf(pd,hist.time_edges),'k')
            
            ylabel('Spike count');
            title(['#' num2str(c) ':' cluster(c).region ', ' convertStringsToChars(dir) 'ward ' convertStringsToChars(state)])
            if i >= 2*N*no_states-1
                xlabel('Time (sec)')
            end
            xlim(timerange)
            if mod(i,2*no_states)==0
                linkaxes([a(i) a(i-1) a(i-2) a(i-3)],'y')
                ylim(max([0 5],ylim)); % start from 0 when no cell rate map is null
                states = repmat(unique([lap.status]),2,1); % alternating between ditch and jump for gaps
                states = states(end:-1:1); % correcting the order
                directions = ["right" "left" "right" "left"]; % alternating between left and right
                %                 for j=i:-1:i-3
                %                     subplot(N*no_states,2,j);
                %
                %                     for l=[lap([lap.dir]==directions(i-j+1) & [lap.status]==states(i-j+1)).no]  % only the desired direction and state
                %
                %
                %                         for k=1:length(lap(l).t_cross)
                %                             plot(repmat(lap(l).t_cross(k),1,2)-lap(l).t_interest, ylim,lap(l).cross_color,'LineWidth',0.15);
                %                         end
                %                     end
                %
                %                     % Reverse the stacking order so that the patch overlays the line
                %                     chi=get(gca, 'Children');
                %                     set(gca, 'Children',flipud(chi));
                %                 end
            end
            
        end
    end
end
cluster = rmfield(cluster,'t_rel');
%toc
set(gcf, 'Position', [0 0 1800 300+200*N*no_states]);
linkaxes(a,'x')
zoom xon
xlim(timerange)
sgtitle(['Directional time histogram (spikes per bin per lap), at the time of ' t_reference]);
saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_time_histogram-cl' mat2str(cluster_no) '_' t_reference '.png']))

%% X-Z plot
figure(4000); clf; a = [];
i = 0;
for c = cluster_no
    
    for state = unique([lap.status])
        for dir = ["left" "right"]
            i = i+1;
            
            a(i) = subplot(no_states,2,i); hold on
            idx = [cluster(c).status]==state & [cluster(c).dir]==dir & abs([cluster(c).s]) >= s_thresh;
            dx = [lap(cluster(c).lap(idx)).dx_interest];
            dx = reshape(dx,[],1);
            %             p = plot(cluster(c).p(idx,1)+dx,cluster(c).p(idx,3),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
            s = scatter(cluster(c).p(idx,1)+dx,cluster(c).p(idx,3),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
            alpha(s,opacity)
            %             idx = pos.status==state & pos.dir==dir & abs(pos.s) >= s_thresh;
            %             dx = [lap(pos.lap(idx)).dx_interest];
            %             dx = reshape(dx,[],1);
            %             h = plot(pos.p(idx,1)+dx,pos.p(idx,3),'.', 'MarkerSize',0.2);
            %             set(h, 'Color', '#D0D0D0');
            
            axis equal
            %             title(['#' num2str(c) ':' cluster(c).region ', Place field (' convertStringsToChars(dir) 'ward ' convertStringsToChars(state) ')'])
            title(['Place field(s) (' convertStringsToChars(dir) 'ward ' convertStringsToChars(state) ')'])
            if i >= 2*no_states-1
                xlabel('Horizontal position (cm)')
            end
            ylabel('Elevation (cm)')
            
            if mod(i,2*no_states)==0
                linkaxes([a(i) a(i-1) a(i-2) a(i-3)])
                ylim([-35 20])
                if isZoom
                    xlim([-50 20+max([lap.gap_length])]+max([lap.dx_interest]))
                else
                    xrange = xlim;
                    xlim([min(xrange(1),-50) max(xrange(2),20+max([lap.gap_length])+max([lap.dx_interest]))])
                end
                states = repmat(unique([lap.status]),2,1); % alternating between ditch and jump for gaps
                states = states(end:-1:1); % correcting the order
                directions = ["right" "left" "right" "left"]; % alternating between left and right
                for j=i:-1:i-3
                    subplot(no_states,2,j);
                    
                    for l=[lap([lap.dir]==directions(i-j+1) & [lap.status]==states(i-j+1)).no]  % only the desired direction and state
                        % gap range
                        lim = xlim;
                        plot([lim(1) lap(l).dx_interest], [0 0],lap(l).cross_color,'LineWidth',0.25);
                        plot([lap(l).dx_interest lap(l).dx_interest], [-depth 0],lap(l).cross_color,'LineWidth',0.25);
                        plot(repmat(lap(l).gap_length+lap(l).dx_interest,1,2), [-depth 0],lap(l).cross_color,'LineWidth',0.25);
                        plot([lap(l).dx_interest lap(l).gap_length+lap(l).dx_interest], [-depth -depth],lap(l).cross_color,'LineWidth',0.25);
                        plot([lap(l).gap_length+lap(l).dx_interest lim(2)+lap(l).dx_interest], [0 0],lap(l).cross_color,'LineWidth',0.25);
                    end
                    %Reverse the stacking order so that the patch overlays the line
                    chi=get(gca, 'Children');
                    set(gca, 'Children',flipud(chi));
                end
                i=i-2*no_states; % overlay
            end
            
        end
    end
end
linkaxes(a)
%toc
set(gcf, 'Position', [0 0 1800 300+200*no_states]);
sgtitle(['Elevation versus Horizontal position, ' x_reference ' frame of ref (speed filtered: s >= ' num2str(s_thresh) ' cm/s)']);
saveas(gcf,fullfile(exp_directory, 'Analysis',['Z-X-cl' mat2str(cluster_no) '_' x_reference '.png']))

%% 3D plot
figure(4001); clf; a = [];
i = 0;
for c = cluster_no
    
    for state = unique([lap.status])
        for dir = ["left" "right"]
            i = i+1;
            
            a(i) = subplot(no_states,2,i); hold on
            idx = [cluster(c).status]==state & [cluster(c).dir]==dir & abs([cluster(c).s]) >= s_thresh;
            dx = [lap(cluster(c).lap(idx)).dx_interest];
            dx = reshape(dx,[],1);
            %plot3(cluster(c).p(idx,1)+dx,cluster(c).p(idx,2),cluster(c).p(idx,3),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
            s = scatter3(cluster(c).p(idx,1)+dx,cluster(c).p(idx,2),cluster(c).p(idx,3),'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
            alpha(s,opacity)
            %             idx = posi.status==state & posi.dir==dir & abs(posi.s) >= s_thresh;
            %             dx = [lap(posi.lap(idx)).dx_interest];
            %             dx = reshape(dx,[],1);
            %             plot3(posi.p(idx,1)+dx,posi.p(idx,2),posi.p(idx,3),'.', 'MarkerSize',0.2,'Color', '#DDDDDD');
            
            %axis equal
            %legend(['#' num2str(c) ':' cluster(c).region]);
            
            %             if i <= 2*no_states
            %                 title(['Firing field (' convertStringsToChars(dir) 'ward)'])
            %             end
            title(['Place field(s) (' convertStringsToChars(dir) 'ward ' convertStringsToChars(state) ')'])
            
            xlabel('X position (cm)')
            ylabel('Y Position (cm)')
            zlabel('Z position (cm)')
            
            if mod(i,2*no_states)==0
                Link = linkprop(a,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
                setappdata(gcf, 'StoreTheLink', Link);
                view([0.2 -1 0.5])
                zlim([-35 20])
                if isZoom
                    xlim([-50 50+max([lap.gap_length])])
                else
                    xrange = xlim;
                    xlim([min(xrange(1),-50) max(xrange(2),50+max([lap.gap_length]))])
                end
                states = repmat(unique([lap.status]),2,1); % alternating between ditch and jump for gaps
                states = states(end:-1:1); % correcting the order
                directions = ["right" "left" "right" "left"]; % alternating between left and right
                for j=i:-1:i-3
                    subplot(no_states,2,j);
                    lim = xlim;
                    [x,y] = meshgrid(-150:150,-7:7);
                    z = - depth * (x > 0) + depth * (x > max([lap.gap_length]+[lap.dx_interest])); %
                    z(x==1)=nan;
                    z(x==floor(max([lap.gap_length]+[lap.dx_interest])))=nan;
                    hold on
                    surf(x,y,z,'FaceAlpha',0.25,'EdgeAlpha',0.1,'FaceColor','#D0D0D0'); hold on
                    
                    
                end
                i=i-2*no_states; % overlay
            end
            
        end
    end
end

%Reverse the stacking order so that the patch overlays the line
N = length(cluster_no);
chi=get(gca, 'Children');
subplot(no_states,2,4);
try
    leg = [legendCell(legendCell~="") strings(N,1)]';
    legend(leg(:),'Location','southeast')
catch
end
Link = linkprop(a,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);
view([-13 7])
%toc
set(gcf, 'Position', [0 0 1800 300+200*no_states]);
sgtitle(['3D place fields, ' x_reference ' frame of ref (speed filtered: s >= ' num2str(s_thresh) ' cm/s)']);
saveas(gcf,fullfile(exp_directory, 'Analysis',['3D-cl' mat2str(cluster_no) '_' x_reference '.png']))

%% phase histogram

% binning
hist.phase_bin_size = 10; % deg
hist.phase_edges = -180:hist.phase_bin_size:180;

N = length(cluster_no);
figure(5000); clf;
i = 0;
for c = cluster_no
    ylimmax = 0;
    for state = unique([lap.status])
        for dir = ["left" "right"]
            i = i+1;
            
            idx = [cluster(c).dir] == dir & [cluster(c).status] == state & abs([cluster(c).s]) >= s_thresh;
            
            n = length([lap([lap.dir]==dir & [lap.status]==state).no]); % number of laps
            
            hist.cluster = histcounts(cluster(c).phase(idx), hist.phase_edges); % spikes in each bin
            hist.phase_hist = hist.cluster / (hist.phase_bin_size*n) ; % spikes per bin per lap
            
            a(i) = subplot(N*no_states,2,i);
            h = histogram('BinCounts', hist.phase_hist, 'BinEdges', hist.phase_edges, 'FaceColor',colors(c),'FaceAlpha',0.8); hold on; %rate map histogram (20% transparent)
            
            ylabel('Spike count');
            title(['#' num2str(c) ':' cluster(c).region ', ' convertStringsToChars(dir) 'ward ' convertStringsToChars(state)])
            if i >= 2*N*no_states-1
                xlabel('Phase (deg)')
            end
            if nnz(idx)>0
                ylimmax = max([ylimmax, ylim]);
            end
            if mod(i,2*no_states)==0
                linkaxes([a(i) a(i-1) a(i-2) a(i-3)],'y')
                ylim([0, ylimmax]);
%                 %ylim(max([0 5],ylim)); % start from 0 when no cell rate map is null
%                 states = repmat(unique([lap.status]),2,1); % alternating between ditch and jump for gaps
%                 states = states(end:-1:1); % correcting the order
%                 directions = ["right" "left" "right" "left"]; % alternating between left and right
                %                 for j=i:-1:i-3
                %                     subplot(N*no_states,2,j);
                %
                %                     for l=[lap([lap.dir]==directions(i-j+1) & [lap.status]==states(i-j+1)).no]  % only the desired direction and state
                %
                %
                %                         for k=1:length(lap(l).t_cross)
                %                             plot(repmat(lap(l).t_cross(k),1,2)-lap(l).t_interest, ylim,lap(l).cross_color,'LineWidth',0.15);
                %                         end
                %                     end
                %
                %                     % Reverse the stacking order so that the patch overlays the line
                %                     chi=get(gca, 'Children');
                %                     set(gca, 'Children',flipud(chi));
                %                 end
            end
            
        end
    end
end

%toc
set(gcf, 'Position', [0 0 1800 300+200*N*no_states]);
linkaxes(a,'x')
%zoom xon
xlim([-180 180])
sgtitle(['Directional phase histogram (spikes per bin per lap), at the time of ' t_reference '(speed filtered: s >= ' num2str(s_thresh) ' cm/s)']);
saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_phase_histogram-cl' mat2str(cluster_no) '_' t_reference '.png']))