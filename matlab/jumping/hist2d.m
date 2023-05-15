%   See also ANALYZEDATA, PLOTTIME.
%
%   Date 2023-01-03
%
%   Author Shahin G Lashkari
%
clc; clear; close all;
answer = inputdlg({'Rat', 'Date'},'Enter the rat number and the date',[1 30],{'1068', '2022-12-20'});
rat_no = answer{1};
date_str = answer{2};
%% Selecting the appropriate files
[datafile,exp_directory] = uigetfile(['E:\Rat' rat_no '\Analysis\' date_str '\processed_data.mat'],'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename);
c = 3;
%%
figure(1); clf
i = 0;
ax = zeros(4,1);
cmax = 0;
X = posi.p(:, 1) + posi.dx';
Z = posi.p(:, 3);
[~, Xedges, Zedges] = histcounts2(X, Z, 'BinWidth', [3 3]);

for state = ["jump" "ditch"]
    for dir = ["left" "right"]
        i = i + 1;
        ax(i) = subplot(2,2,i);
        
        idx = posi.status==state & posi.dir==dir;
        X = posi.p(idx, 1) + posi.dx(idx)';
        Z = posi.p(idx, 3);
        occupancy_counts = histcounts2(X, Z,'XBinEdges',Xedges,'YBinEdges',Zedges);
        occupancy_counts = occupancy_counts * posi.dt;  % seconds in each bin
        
        idx = [cluster(c).status]==state & [cluster(c).dir]==dir;
        x = cluster(c).p(idx, 1) + posi.dx(idx)';
        z = cluster(c).p(idx, 3);
        counts = histcounts2(x, z,'XBinEdges',Xedges,'YBinEdges',Zedges); % spikes in each bin
        
        rate = counts ./ (occupancy_counts + eps);
        % rate(rate > 100) = 0;
        rate(occupancy_counts < 0.1) = 0;
        
        histogram2('DisplayStyle','tile','ShowEmptyBins','on', 'XBinEdges',Xedges,'YBinEdges',Zedges,'BinCounts',rate);
        colorbar
        cmax = max(cmax, max(caxis));
        title('Rate map')
        xlabel('Horizontal Position (cm)');
        ylabel('Elevation (cm)');
        axis equal
    end
end
i = 0;
for state = ["jump" "ditch"]
    for dir = ["left" "right"]
        i = i + 1;
        subplot(2,2,i);
        caxis([0 cmax]);
    end
end
linkaxes(ax, 'xy')
sgtitle(['#' num2str(c) ':' cluster(c).region]);
%%
figure(2); clf
i = 0;
ax = zeros(4,1);
cmax = 0;

X = posi.x + posi.dx;
PHASE = posi.phase;
[~, Xedges, Zedges] = histcounts2(X, PHASE, 'BinWidth', [3 18]);

for state = ["jump" "ditch"]
    for dir = ["left" "right"]
        i = i + 1;
        ax(i) = subplot(2,2,i);
        idx = posi.status==state & posi.dir==dir;
        
        X = posi.x(idx) + posi.dx(idx);
        PHASE = posi.phase(idx);
        occupancy_counts = histcounts2(X, PHASE,'XBinEdges',Xedges,'YBinEdges',Zedges);
        occupancy_counts = occupancy_counts * posi.dt;  % seconds in each bin
        
        idx = [cluster(c).status]==state & [cluster(c).dir]==dir;
        x = cluster(c).x(idx) + cluster(c).dx(idx);
        phase = cluster(c).phase(idx);
        counts = histcounts2(x, phase,'XBinEdges',Xedges,'YBinEdges',Zedges); % spikes in each bin
        
        rate = counts ./ (occupancy_counts + eps);
        % rate(rate > 100) = 0;
        rate(occupancy_counts < 0.1) = 0;
        
        histogram2('DisplayStyle','tile','ShowEmptyBins','on', 'XBinEdges',Xedges,'YBinEdges',Zedges,'BinCounts',rate); hold on;
        histogram2('DisplayStyle','tile','ShowEmptyBins','on', 'XBinEdges',Xedges,'YBinEdges',Zedges+360,'BinCounts',rate);
        colorbar
        title('Rate map');
        xlabel('Horizontal Position (cm)');
        ylabel('Phase (angle)');
        cmax = max(cmax, max(caxis));
    end
end
i = 0;
for state = ["jump" "ditch"]
    for dir = ["left" "right"]
        i = i + 1;
        subplot(2,2,i);
        caxis([0 cmax]);
    end
end
linkaxes(ax, 'xy')
sgtitle(['#' num2str(c) ':' cluster(c).region]);