%%% Rate map calculations and plots
% SGL 2022-01-31
clc; clear; close all
[datafile,exp_directory] = uigetfile(fullfile('D:\Analysis', 'analyzed_data.mat'), 'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename, 'pos', 'posi', 'lap', 'cluster','ppcm', 'colors','xmax','hist');
start = tic;
cluster_no =[3 13 17 18 22 23 28 29 31]; %[cluster.no];
showlapbylap = true;
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
    quiver(posi.x(idx) * ppcm,posi.y(idx) * ppcm,posi.vx(idx) * ppcm,posi.vy(idx) * ppcm,2);
    
    % occupancy histogram
    hist.posi = histcounts(posi.x(posi.lap==l), hist.edges) * posi.dt; % seconds in each bin
    ax2 = subplot(5,1,3);
    histogram('BinCounts', hist.posi, 'BinEdges', hist.edges);
    ylabel('Occupancy (sec)')
    xlim([0 xmax])
    
    for c=cluster_no % each cluster in cluster_no
        ax1 = subplot(5,1,1:2);
        h = plot(cluster(c).x([cluster(c).lap]==l) * ppcm, cluster(c).y([cluster(c).lap]==l) * ppcm,...
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
        xlim([0 xmax])
        
        saveas(gcf,[exp_directory filesep 'Analysis' filesep 'cluster' num2str(c) '_lap' num2str(l) '.jpg'])
        %set(h, 'Visible','off')
        delete(h);
    end
end

%% Directional rate map for all the laps
% groups of 14
for i = 1:ceil(length(cluster)/14)
    figure(200+i); clf;
    for j=1:14
        c = (i-1) * 14 + j;
        if c > length(cluster)
            continue;
        end
        % leftward rate mapj
        hist.posi = histcounts(posi.x(posi.dir=="left"), hist.edges) * posi.dt; % seconds in each bin
        hist.cluster = histcounts(cluster(c).x([cluster(c).dir]=="left"), hist.edges); % spikes in each bin
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        
        a(2*j-1) = subplot(14,2,14*2-(2*j-1));
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
        ylabel(['#' num2str(c) ':' cluster(c).region]);
        if j == 14
            title('Rate Map (leftward)')
        elseif j ==1
            xlabel('Horizontal position (cm)')
        end
        
        % rightward rate map
        hist.posi = histcounts(posi.x(posi.dir=="right"), hist.edges) * posi.dt; % seconds in each bin
        hist.cluster = histcounts(cluster(c).x([cluster(c).dir]=="right"), hist.edges); % spikes in each bin
        hist.ratemap = hist.cluster ./ (hist.posi + eps); % adding eps to avoid division by zero
        
        a(2*j) = subplot(14,2,14*2-(2*j-1)+1);
        histogram('BinCounts', hist.ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); %rate map histogram
        ylabel(['#' num2str(c) ':' cluster(c).region]);
        linkaxes([a(2*j-1) a(2*j)],'y')
        ylim(max([0 5],ylim))
        if j == 14
            title('Rate Map (rightward)')
        elseif j ==1
            xlabel('Horizontal position (cm)')
        end
        
        % edge of gap
        hold on
        plot(repmat(lap(l).gap(1),2,1), ylim,'b','LineWidth',2);
        plot(repmat(lap(l).gap(2),2,1), ylim,'b','LineWidth',2);
        a(2*j-1) = subplot(14,2,14*2-(2*j-1));
        hold on
        plot(repmat(lap(l).gap(1),2,1), ylim,'b','LineWidth',2);
        plot(repmat(lap(l).gap(2),2,1), ylim,'b','LineWidth',2);
    end
    
    set(gcf, 'Position', [100 100 1600 1100]);
    linkaxes(a,'x')
    zoom xon
    xlim([0 xmax])
    saveas(gcf,fullfile(exp_directory, 'Analysis',['Directional_ratemap-' num2str(i) '.jpg']))
end
%%
fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);