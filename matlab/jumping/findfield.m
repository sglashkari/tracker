% function range = findfield(data, hist)
clc; clear; close all;
exp_directory = 'E:\Rat1068\Analysis\2022-12-20\';
% exp_directory = 'E:\Rat1055\Analysis\2022-11-09\';
% exp_directory = 'E:\Rat980\Analysis\2021-12-21\';
load(fullfile(exp_directory,'processed_data.mat'), 'pos', 'posi', 'lap', 'cluster', 'colors','hist','exp');
%%
hist.centers = movmean(hist.edges, 2, 'Endpoints', 'discard');
n = 10;

hist.centers = hist.centers(1)-hist.bin_size*n:hist.bin_size: hist.centers(end)+hist.bin_size*n;

for c = [cluster.no]
    figure(c); clf;
    
    ratemap = mean([cluster(c).hist.stationary.left_jump; cluster(c).hist.stationary.right_jump;  ...
        cluster(c).hist.stationary.left_ditch;  cluster(c).hist.stationary.right_ditch]);
    
    a(1) = subplot(2,1,1);
    histogram('BinCounts', ratemap, 'BinEdges', hist.edges, 'FaceColor',colors(c)); hold on; %rate map histogram (20% transparent)
    ylim(max([0 5],ylim));
    a(2) = subplot(2,1,2);
    ratemap = [zeros(1,n) ratemap zeros(1,n)];
    findpeaks(movmean(ratemap,n), hist.centers,'MinPeakDistance',50, 'MinPeakProminence', 0.5,'Annotate','extents','WidthReference','halfheight');
    ylim(max([0 5],ylim));
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
    [pks,locs,w,p] = findpeaks(movmean(ratemap,n), hist.centers,'MinPeakDistance',50, 'MinPeakProminence', 0.5,'Annotate','extents','WidthReference','halfheight');
    linkaxes(a, 'x')
    xlim([0 1.1*exp.xmax])
    
    no_fields = numel(pks);
    if no_fields < 1
        continue;
    end
    
    range(no_fields,2) = 0;
    for f=1:no_fields
        range(f,:) = locs(f) + [-w(f) w(f)];
    end
    a(1) = subplot(2,1,1);
    
    for f = 2:no_fields
        if range(f-1,2) > range(f,1)
            for i=1:numel(border)
                if border(i) > locs(f-1) && border(i) < locs(f)
                    range(f-1,2) = min(border(i), range(f-1,2));
                    range(f,1) = max(border(i), range(f,1));
                end
            end
        end
    end
    
    data = [];
    
    for i = 1:numel(ratemap)
        for j = 1:round(4*ratemap(i)) % 4x precesion
            data = [data hist.centers(i)];
        end
    end
    
    dash = ["-." "--" ":"];
    for f=1:no_fields
        plot(locs(f) * [1 1],ylim,'k' + dash(f))
        plot(range(f,1) * [1 1],ylim, 'g'+dash(f))
        plot(range(f,2) * [1 1],ylim, 'g'+dash(f))
        p50(f) = prctile(data(data>range(f,1) & data<range(f,2)),50);
        p2(f) = prctile(data(data>range(f,1) & data<range(f,2)),2);
        p98(f) = prctile(data(data>range(f,1) & data<range(f,2)),98);
        plot(p50(f) * [1 1],ylim,'b' + dash(f))
        plot(p2(f) * [1 1],ylim,'m' + dash(f))
        plot(p98(f) * [1 1],ylim,'m' + dash(f))
        area(f) = sum(data(data>range(f,1) & data<range(f,2)));
    end
    
end