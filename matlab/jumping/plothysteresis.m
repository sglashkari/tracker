%%% Hysteresis calculations and plots
% SGL 2022-02-10
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
hist.edges = 10:2:40;
idx = ([lap.gap_status]=='increase') & ([lap.status] == 'jump') & ([lap.dir] == 'right');
hist.jump = histcounts([lap(idx).gap_length], hist.edges); % jumps in each bin
idx = ([lap.gap_status]=='increase') & ([lap.status] == 'ditch') & ([lap.dir] == 'right');
hist.ditch = histcounts([lap(idx).gap_length], hist.edges); % ditches in each bin
hist.passage = hist.jump+hist.ditch;
figure(2);
hist.probability = round(hist.jump ./ (hist.passage + eps),1);
histogram('BinCounts', hist.probability, 'BinEdges', hist.edges); % probability histogram
figure(4);
plot(hist.edges(hist.passage>=min_passage),hist.probability(hist.passage>=min_passage),'--'); hold on
hist
sum(hist.passage)
%% increase & leftward
idx = ([lap.gap_status]=='increase') & ([lap.status] == 'jump') & ([lap.dir] == 'left');
hist.jump = histcounts([lap(idx).gap_length], hist.edges); % jumps in each bin
idx = ([lap.gap_status]=='increase') & ([lap.status] == 'ditch') & ([lap.dir] == 'left');
hist.ditch = histcounts([lap(idx).gap_length], hist.edges); % ditches in each bin
hist.passage = hist.jump+hist.ditch;
figure(3);
hist.probability = round(hist.jump ./ (hist.passage + eps),1);
histogram('BinCounts', hist.probability, 'BinEdges', hist.edges); % probability histogram
figure(4);
plot(hist.edges(hist.passage>=min_passage),hist.probability(hist.passage>=min_passage),'--'); hold on
hist
sum(hist.passage)
%% decrease & rightward
idx = ([lap.gap_status]=='decrease') & ([lap.status] == 'jump') & ([lap.dir] == 'right');
hist.jump = histcounts([lap(idx).gap_length], hist.edges); % jumps in each bin
idx = ([lap.gap_status]=='decrease') & ([lap.status] == 'ditch') & ([lap.dir] == 'right');
hist.ditch = histcounts([lap(idx).gap_length], hist.edges); % ditches in each bin
hist.passage = hist.jump+hist.ditch;
figure(12);
hist.probability = round(hist.jump ./ (hist.passage + eps),1);
histogram('BinCounts', hist.probability, 'BinEdges', hist.edges); % probability histogram
figure(4);
plot(hist.edges(hist.passage>=min_passage),hist.probability(hist.passage>=min_passage),'--'); hold on
hist
sum(hist.passage)
%% decrease & leftward
idx = ([lap.gap_status]=='decrease') & ([lap.status] == 'jump') & ([lap.dir] == 'left');
hist.jump = histcounts([lap(idx).gap_length], hist.edges); % jumps in each bin
idx = ([lap.gap_status]=='decrease') & ([lap.status] == 'ditch') & ([lap.dir] == 'left');
hist.ditch = histcounts([lap(idx).gap_length], hist.edges); % ditches in each bin
hist.passage = hist.jump+hist.ditch;
figure(13);
hist.probability = round(hist.jump ./ (hist.passage + eps),1);
histogram('BinCounts', hist.probability, 'BinEdges', hist.edges); % probability histogram
figure(4);
plot(hist.edges(hist.passage>=min_passage),hist.probability(hist.passage>=min_passage),'--');
hist
sum(hist.passage)
ylim([-0.1 1.1]);
%xlim([10 40])