%%% Hysteresis calculations and plots
% SGL 2022-02-16
clc; clear; close all
[datafile,exp_directory] = uigetfile(fullfile('D:\Analysis', 'analyzed_data.mat'), 'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename, 'lap','pos');
start = tic;
min_passage = 4;
%% 
lap(1).gap_status = "increase";
for l=2:length(lap)
    if lap(l).gap(2) > lap(l-1).gap(2) + 0.5
        lap(l).gap_status = "increase";
    elseif lap(l).gap(2) < lap(l-1).gap(2) - 0.5
        lap(l).gap_status = "decrease";
    else
        lap(l).gap_status = lap(l-1).gap_status;
    end
end

figure(1); clf
%plot([lap.t_jump],[lap.gap_length]);
plot(pos.t,pos.len)
xlim([0 max(pos.t)])
ylabel('Gap Length (cm)')
xlabel('Time (sec)')
set(gcf, 'Position', [100 100 1800 600]);
saveas(gcf,fullfile(exp_directory, 'Analysis',['GapLength-Time.png']))

%% all 4 combinations
Legend=cell(4,1);
i = 0;
hist.edges = 10:2:40;
for increment = ["increase" "decrease"]
    for dir = ["right" "left"]
        i = i + 1;
        idx = ([lap.gap_status]==increment) & ([lap.status] == 'jump') & ([lap.dir] == dir);
        hist.jump = histcounts([lap(idx).gap_length], hist.edges); % jumps in each bin
        idx = ([lap.gap_status]==increment) & ([lap.status] == 'ditch') & ([lap.dir] == dir);
        hist.ditch = histcounts([lap(idx).gap_length], hist.edges); % ditches in each bin
        hist.passage = hist.jump+hist.ditch;
        figure(10+i);
        hist.probability = round(hist.jump ./ (hist.passage + eps),1);
        histogram('BinCounts', hist.probability, 'BinEdges', hist.edges); % probability histogram
        title(increment+"  "+dir);
        figure(2);
        subplot(2,1,mod(i-1,2)+1)
        title(dir+"ward")
        plot(hist.edges(hist.passage>=min_passage),hist.probability(hist.passage>=min_passage),'-o'); hold on
        hist
        %Legend{mod(i-1,2)+1}=increment+" & "+dir+"ward"; %[convertStringsToChars(increment) ' & ' convertStringsToChars(dir) 'ward'] ;
        %legend(Legend);
        sum(hist.passage)
    end
end
for i=1:2
    subplot(2,1,mod(i-1,2)+1)
ylim([-0.1 1.1]);
xlim([15 36])
legend(["increase" "decrease"]);
end

%% lumped version
figure(3); clf
i = 0;
hist.edges = 10:2:40;
for increment = ["increase" "decrease"]
    i = i + 1;
    idx = ([lap.gap_status]==increment) & ([lap.status] == 'jump');
    hist.jump = histcounts([lap(idx).gap_length], hist.edges); % jumps in each bin
    idx = ([lap.gap_status]==increment) & ([lap.status] == 'ditch');
    hist.ditch = histcounts([lap(idx).gap_length], hist.edges); % ditches in each bin
    hist.passage = hist.jump+hist.ditch;
%     figure(10+i);
    hist.probability = round(hist.jump ./ (hist.passage + eps),1);
%     histogram('BinCounts', hist.probability, 'BinEdges', hist.edges); % probability histogram
    title(increment+"  "+dir);
%     figure(3);
    title('lumped')
    plot(hist.edges(hist.passage>=min_passage),hist.probability(hist.passage>=min_passage),'-o'); hold on
    hist
    sum(hist.passage)
end
ylim([-0.1 1.1]);
xlim([14 35])
legend(["increase" "decrease"]);