function extractfields(exp_directory)
%% extract fields information from the cluster ratemaps
%
%   See also PLOTRATEMAP, PLOTTIME, PLOT_FIELDS.
%
%   Date 2023-01-13
%   Author Shahin G Lashkari
%
%% Selecting the appropriate files
if nargin < 1
    clc; clear; close all
    answer = inputdlg({'Rat', 'Date'},'Enter the rat number and the date',[1 30],{'1068', '2022-12-20'});
    rat_no = answer{1};
    date_str = answer{2};
    [datafile,exp_directory] = uigetfile(['E:\Rat' rat_no '\Analysis\' date_str '\processed_data.mat'],'Select Data File');
    if isequal(datafile, 0)
        error('Data file was not selected!')
    end
    mat_filename = fullfile(exp_directory, datafile);
    load(mat_filename, 'lap', 'cluster', 'colors','hist','exp');
else
    mat_filename = fullfile(exp_directory, 'processed_data.mat');
    load(mat_filename, 'lap', 'cluster', 'colors','hist','exp');
end

%%
start = tic;
field = [];
n = 10;
centers = hist.centers(1)-hist.bin_size*n:hist.bin_size: hist.centers(end)+hist.bin_size*n;
x_reference = "stationary";
for c = [cluster.no]
    figure(c); clf;
    i = 0;
    this_field = [];
    clear ax
    for dir = ["left" "right"]
        i = i + 1;
        
        ratemap = mean([cluster(c).hist.(x_reference).([convertStringsToChars(dir) '_jump']); ...
            cluster(c).hist.(x_reference).([convertStringsToChars(dir) '_ditch'])]);
        
        ax(i) = subplot(2,2,i);
        histogram('BinCounts', ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); hold on; %rate map histogram (20% transparent)
        ylabel('Firing rate (Hz)')
        ylim(max([0 5],ylim));
        title(dir)
        xlabel('Horizontal position (cm)')
        ax(i+2) = subplot(2,2,i+2);
        ratemap = [zeros(1,n) ratemap zeros(1,n)];
        findpeaks(movmean(ratemap,n), centers,'MinPeakDistance',50, 'MinPeakProminence', 0.5,'Annotate','extents','WidthReference','halfheight');
        ylim(max([0 5],ylim));
        title(dir)
        xlabel('Horizontal position (cm)')
        Ax = gca;
        Kids = Ax.Children;
        if numel(Kids) >= 5
            Borders = Kids(1:2);
            Line = Borders.XData;
            border = Line(1:3:end);
        else
            Line = [];
            border = exp.xmax;
        end
        [pks,locs,w] = findpeaks(movmean(ratemap,n), centers,'MinPeakDistance',50, 'MinPeakProminence', 0.5,'Annotate','extents','WidthReference','halfheight');
        
        xlim([0 exp.xmax])
        ylabel('Firing rate (Hz)')
        no_fields.(dir) = numel(pks);
        if no_fields.(dir) >= 1
            
            
            range= zeros(no_fields.(dir), 2);
            for f=1:no_fields.(dir)
                range(f,:) = locs(f) + [-w(f) w(f)];
            end
            ax(i) = subplot(2,2,i);
            
            for f = 2:no_fields.(dir)
                if range(f-1,2) > range(f,1)
                    for j=1:numel(border)
                        if border(j) > locs(f-1) && border(j) < locs(f)
                            range(f-1,2) = min(border(j), range(f-1,2));
                            range(f,1) = max(border(j), range(f,1));
                        end
                    end
                end
            end
            
            data = [];
            
            for j = 1:numel(ratemap)
                for k = 1:round(4*ratemap(j)) % 4x precesion
                    data = [data centers(j)];
                end
            end
            
            dash = ["-." "--" ":"];
            for f=1:no_fields.(dir)
                plot(locs(f) * [1 1],ylim,'k' + dash(f))
                plot(range(f,1) * [1 1],ylim, 'g'+dash(f))
                plot(range(f,2) * [1 1],ylim, 'g'+dash(f))
                p50(f) = prctile(data(data>range(f,1) & data<range(f,2)),50);
                p2(f) = prctile(data(data>range(f,1) & data<range(f,2)),2);
                p98(f) = prctile(data(data>range(f,1) & data<range(f,2)),98);
                plot(p50(f) * [1 1],ylim,'b' + dash(f))
                plot(p2(f) * [1 1],ylim,'m' + dash(f))
                plot(p98(f) * [1 1],ylim,'m' + dash(f))
            end
            
            if dir == "left"
                for f = 1:no_fields.left
                    this_field(f).cluster = c;
                    this_field(f).center = p50(f);
                    this_field(f).range = [p2(f) p98(f)];
                    this_field(f).peak = pks(f);
                    this_field(f).dir = dir;
                end
            else
                for f = 1:no_fields.right
                    this_field(no_fields.left+f).cluster = c;
                    this_field(no_fields.left+f).center = p50(f);
                    this_field(no_fields.left+f).range = [p2(f) p98(f)];
                    this_field(no_fields.left+f).peak = pks(f);
                    this_field(no_fields.left+f).dir = dir;
                end
            end
            
        end
    end
    linkaxes(ax, 'x')
    xlim([0 exp.xmax])
    
    % check if the rightward and leftward fields are the same field
    common.right = zeros(no_fields.right, 1);
    common.left = zeros(no_fields.left, 1);
    for fl = 1:no_fields.left
        for fr = 1:no_fields.right
            range_r = this_field(no_fields.left+fr).range;
            range_l = this_field(fl).range;
            common_range = [max(range_r(1), range_l(1)) min(range_r(2), range_l(2))];
            if diff(common_range)/diff(range_r) >= 0.5 && diff(common_range)/diff(range_l) >= 0.5
                disp(['common field for left and right ' num2str(diff(common_range)/diff(range_r)) ', ' num2str(diff(common_range)/diff(range_l))]);
                common.right(fr) = 1;
                common.left(fl) = 1;
            elseif diff(common_range)/diff(range_r) > 1/3 || diff(common_range)/diff(range_l) > 1/3
                disp(['No common field for left and right ' num2str(diff(common_range)/diff(range_r)) ', ' num2str(diff(common_range)/diff(range_l))]);
            end
        end
    end
    
    for f = 1:no_fields.left
        this_field(f).iscommon = boolean(common.left(f));
        if diff(this_field(f).range) < exp.xmax / 2
            field = [field this_field(f)];
        end
    end
    for f = 1:no_fields.right
        this_field(no_fields.left+f).iscommon = boolean(common.right(f));
        if diff(this_field(no_fields.left+f).range) < exp.xmax / 2
            field = [field this_field(no_fields.left+f)];
        end
    end
    
    sgtitle(['Cluster #' num2str(c) ' ' convertStringsToChars(cluster(c).region) ' : combined stationary Ratemap (jumping and ditching)']);
    set(gcf, 'Position', [100 100 1600 1100]);
    saveas(gcf,[exp_directory filesep 'Analysis' filesep 'field_extraction_cluster_' num2str(c) '.jpg'])
end
% field = rmfield(field, 'combined');

%% choosing x_ref (stationary or moving) ==> detectstatus
close all
detectstatus;

% sort fields based on the location of the peak
[~, idx] = sort([field.center]);
field = field(idx);

for f = 1:numel(field)
    field(f).no = f;
end

%% plot 2 plots (jump, ditch) for each field
no_states = length(unique([lap.status]));
for f = [field.no]
    clear a
    c = field(f).cluster;
    x_reference = "stationary";
    dir = field(f).dir;
    i = 0;
    for state = sort(unique([lap.status]), 'descend')
        i = i+1;
        a(i) = subplot(no_states,1,i);
        ratemap = cluster(c).hist.stationary.([convertStringsToChars(dir) '_' convertStringsToChars(state)]);
        histogram('BinCounts', ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c),'FaceAlpha',0.1); hold on; %rate map histogram (90% transparent)
        ratemap = ratemap .* (hist.centers >= field(f).range(1) & hist.centers <= field(f).range(2));
        histogram('BinCounts', ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c),'FaceAlpha',0.8); hold on; %rate map histogram (20% transparent)
        field(f).hist.(state) = ratemap;
        ylabel({['F#' num2str(f) ' (C#' num2str(c) ')'];cluster(c).region})
        xlim([0 exp.xmax]);
        
        
        
        % calculating area under the curve
        field(f).area.(state) = sum(ratemap);
    end
    
    field(f).hist.combined = mean([field(f).hist.jump; field(f).hist.ditch]);
    field(f).area.combined = sum(field(f).hist.combined);
    
    linkaxes(a,'y')
    ylim(max([0 5],ylim)); % start from 0 when no cell rate map is null
    i = 0;
    for state = sort(unique([lap.status]), 'descend')
        i = i+1;
        subplot(no_states,1,i);
        for l=[lap([lap.dir]==dir & [lap.status]==state).no]  % only the desired direction and state
            % gap range
            dx = lap(l).corr  * (x_reference == "moving");
            plot(repmat(lap(l).gap(1),1,2) + dx, ylim,lap(l).cross_color,'LineWidth',0.25);
            plot(repmat(lap(l).gap(2),1,2) + dx, ylim,lap(l).cross_color,'LineWidth',0.25);
        end
        plot(field(f).range(1) * [1 1],ylim,'m--')
        plot(field(f).range(2) * [1 1],ylim,'m--')
        %Reverse the stacking order so that the patch overlays the line
        chi=get(gca, 'Children');
        set(gca, 'Children',flipud(chi));
    end
    sgtitle(['field ' num2str(f) ' ratemap, ' convertStringsToChars(x_reference) ' frame of ref']);
    set(gcf, 'Position', [100 100 1600 1100]);
    saveas(gcf,[exp_directory filesep 'Analysis' filesep convertStringsToChars(x_reference) '_field' num2str(f) '.jpg'])
    close all
end

% field direction (leftward/rightward) and status (jump/ditch)
for f = [field.no]
    field(f).ratio_jump = field(f).area.jump / (field(f).area.jump + field(f).area.ditch + eps);
    
    if field(f).ratio_jump > 2/3
        field(f).state = "jump";
    elseif field(f).ratio_jump < 1/3
        field(f).state = "ditch";
    else
        field(f).state = "both";
    end
end

%% Saving data (fields)
mat_filename = fullfile(exp_directory,'processed_data.mat');
save(mat_filename, 'field', '-append');
fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);
disp(['File ' mat_filename ' has been updated!'])