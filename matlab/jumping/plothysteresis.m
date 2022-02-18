%%% Hysteresis calculations and plots
% SGL 2022-02-16
clc; clear; close all
[datafile,exp_directory] = uigetfile(fullfile('D:\Analysis', 'analyzed_data.mat'), 'Select Data File');
if isequal(datafile, 0)
    error('Data file was not selected!')
end
mat_filename = fullfile(exp_directory, datafile);
load(mat_filename, 'lap');
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

figure(1);
plot([lap.t_jump],[lap.gap_length]);

%% increase & rightward
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
        figure(2);
        plot(hist.edges(hist.passage>=min_passage),hist.probability(hist.passage>=min_passage),'-.'); hold on
        hist
        Legend{i}=[convertStringsToChars(increment) ' & ' convertStringsToChars(dir) 'ward'] ;
        sum(hist.passage)
    end
end
ylim([-0.1 1.1]);
xlim([10 40])
legend(Legend);