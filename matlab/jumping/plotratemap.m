function cluster = plotratemap(exp_directory, posi, lap, cluster, colors,hist,exp)
%%% Rate map calculations and plots
%
%   See also ANALYZEDATA, EXTRACTFIELDS.
%
%   Date 2023-01-02 (originally 2022-01-31)
%   
%   Author Shahin G Lashkari

%% Selecting the appropriate files
if nargin < 1
    clc; close all
    answer = inputdlg({'Rat', 'Date'},'Enter the rat number and the date',[1 30],{'1068', '2022-12-20'});
    rat_no = answer{1};
    date_str = answer{2};
    [datafile,exp_directory] = uigetfile(['E:\Rat' rat_no '\Analysis\' date_str '\processed_data.mat'],'Select Data File');
    if isequal(datafile, 0)
        error('Data file was not selected!')
    end
    mat_filename = fullfile(exp_directory, datafile);
    load(mat_filename, 'posi', 'lap', 'cluster', 'colors','hist','exp');
else
    mat_filename = fullfile(exp_directory, 'processed_data.mat');
end
%%
start = tic;
cluster_no = [cluster.no];
showlapbylap = false;

%% Plotting the spikes on the path
for l = 1:length(lap)
    if ~showlapbylap
        break;
    end
    % Figures of the rat in the mid-jump
    figure(10+l)
    set(gcf, 'Position', [100 100 1770 1000]);
    ax1 = subplot(5,1,1:2);
    img = imread([exp_directory filesep 'frames' filesep 'frame-' num2str(lap(l).frame) '.pgm']);
    img = imlocalbrighten(img);
    img = medfilt2(img);
    imshow(img);
    figure(10+l)
    hold on
    idx = posi.lap == l;
    quiver(posi.x(idx) * exp.ppcm,posi.y(idx) * exp.ppcm,posi.vx(idx) * exp.ppcm,posi.vy(idx) * exp.ppcm,2);
    
    % occupancy histogram
    hist.posi = histcounts(posi.x(posi.lap==l), hist.edges) * posi.dt; % seconds in each bin
    ax2 = subplot(5,1,3);
    histogram('BinCounts', hist.posi, 'BinEdges', hist.edges);
    ylabel('Occupancy (sec)')
    xlim([0 exp.xmax])
    
    for c=[cluster.no] % each cluster in cluster_no
        ax1 = subplot(5,1,1:2);
        h = plot(cluster(c).x([cluster(c).lap]==l) * exp.ppcm, cluster(c).y([cluster(c).lap]==l) * exp.ppcm,...
            'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
        title([cluster(c).region ': cluster no. ' num2str(c) ', lap ' num2str(l) ', ' ...
            convertStringsToChars(lap(l).dir) 'ward ' convertStringsToChars(lap(l).status)], 'Interpreter', 'none');
        
        % spike histogram
        ax3 = subplot(5,1,4);
        hist.cluster = histcounts(cluster(c).x([cluster(c).lap]==l), hist.edges); % spikes in each bin
        histogram('BinCounts', hist.cluster, 'BinEdges', hist.edges);
        ylim(max([0 5],ylim))
        ylabel('Number of spikes')
        
        % Rate map
        ax4 = subplot(5,1,5);
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
        ylim(max([0 5],ylim))
        ylabel('Rate Map (Hz)')
        xlabel('Horizontal position (cm)')
        
        linkaxes([ax2 ax3 ax4],'x')
        xlim([0 exp.xmax])
        
        saveas(gcf,[exp_directory filesep 'Analysis' filesep 'cluster' num2str(c) '_lap' num2str(l) '.jpg'])
        %set(h, 'Visible','off')
        delete(h);
    end
end

%% Rat-map with speed filter (for stationary or moving frame)
N = numel([cluster.no]);
M = min([10 N]); % groups of 10 max
no_states = length(unique([lap.status]));

for x_reference = ["stationary" "moving"]
    clear a
    i = 0;
    j = 1;
    if x_reference == "stationary"
        figure(2001); clf;
    else
        figure(3001); clf;
    end
    for c = [cluster.no]
        
        for state = sort(unique([lap.status]), 'descend')
            for dir = ["left" "right"]
                i = i+1;
                idx = posi.status==state & posi.dir==dir & abs(posi.s) >= hist.s_thresh;
                dx = posi.dx * (x_reference == "moving");
                hist_posi = histcounts(posi.x(idx) + dx(idx), hist.edges) * posi.dt; % seconds in each bin
                idx = [cluster(c).status]==state & [cluster(c).dir]==dir & abs([cluster(c).s]) >= hist.s_thresh;
                dx = [cluster(c).dx] * (x_reference == "moving");
                hist_cluster = histcounts(cluster(c).x(idx) + dx(idx), hist.edges); % spikes in each bin
                rate_map = hist_cluster ./ (hist_posi + eps); % adding eps to avoid division by zero
                a(i) = subplot(M*no_states,2,i);
                histogram('BinCounts', rate_map, 'BinEdges', hist.edges, 'FaceColor',colors(c),'FaceAlpha',0.8); hold on; %rate map histogram (20% transparent)
                
                cluster(c).hist.(x_reference).([convertStringsToChars(dir) '_' convertStringsToChars(state)]) = rate_map;
                % for large sets
                set(gca,'XTick', []);
                ylabel({['#' num2str(c)];cluster(c).region},'FontSize', 8)
                xlim([0 exp.xmax]);
                
                if mod(i,2*no_states)==0
                    linkaxes([a(i) a(i-1) a(i-2) a(i-3)],'xy')
                    ylim(max([0 5],ylim)); % start from 0 when no cell rate map is null
                    states = repmat(sort(unique([lap.status]), 'descend'),2,1); % alternating between ditch and jump for gaps
                    states = states(end:-1:1); % correcting the order
                    directions = ["right" "left" "right" "left"]; % alternating between left and right
                    for k=i:-1:i-3
                        subplot(M*no_states,2,k);
                        for l=[lap([lap.dir]==directions(i-k+1) & [lap.status]==states(i-k+1)).no]  % only the desired direction and state
                            % gap range
                            dx = lap(l).corr  * (x_reference == "moving");
                            plot(repmat(lap(l).gap(1) + dx,1,2), ylim,lap(l).cross_color,'LineWidth',0.25);
                            plot(repmat(lap(l).gap(2) + dx,1,2), ylim,lap(l).cross_color,'LineWidth',0.25);
                        end
                        %Reverse the stacking order so that the patch overlays the line
                        chi=get(gca, 'Children');
                        set(gca, 'Children',flipud(chi));
                    end
                end
                
                if i == 2*M*no_states || (c == cluster_no(end) && i+(j-1)*2*M*no_states == 2*N*no_states)
                    subplot(M*no_states,2,i);
                    xlabel('Horizontal position (cm)')
                    set(gca,'XTick', 0:50:exp.xmax);
                    subplot(M*no_states,2,i-1);
                    xlabel('Horizontal position (cm)')
                    set(gca,'XTick', 0:50:exp.xmax);
                    set(gcf, 'Position', [0 0 1800 300+200*N*no_states]);
                    sgtitle(['Directional ratemap ' num2str(j) '/' num2str(ceil(N/10)) ', ' convertStringsToChars(x_reference) ' frame of ref (speed filtered: s >= ' num2str(hist.s_thresh) ' cm/s)']);
                    saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_ratemap_' num2str(j) '_' convertStringsToChars(x_reference) '.jpg']));
                    i = 0;
                    j = j + 1;
                    if c < cluster_no(end)
                        if x_reference == "stationary"
                            figure(2000+j); clf;
                        else
                            figure(3000+j); clf;
                        end
                    end
                    a = [];
                end
                
            end
        end
    end
end

%% save file if the functiona in called, otherwise display the info
if nargout == 0
    save(mat_filename, 'cluster', '-append');
    fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);
    disp(['File ' mat_filename ' has been updated!'])
    clear cluster
end
