%% force
%
%   See also PLOTRATEMAP, PLOTTIME.
%
%   Date 2023-01-27
%   Author Shahin G Lashkari
%
clc; clear; close all;
exp_directory = 'E:\Rat1068\Analysis\2022-12-20\';
% exp_directory = 'E:\Rat1055\Analysis\2022-11-09\';
% exp_directory = 'E:\Rat980\Analysis\2021-12-21\';
load(fullfile(exp_directory,'processed_data.mat'), 'daq', 'lap', 'cluster', 'colors','hist','csc');


% for phase calculations
csc_filename = fullfile(exp_directory, 'LFP','LFP.ncs'); % use apbin2lfp and plotcsc to optimize it

%% Take-off
timerange = [-0.3 0.2]; % [-0.5 0.5] 0.5 sec before to 0.5 sec after
for l = 1:numel(lap)
    lap(l).t_max_grf_takeoff = nan;
    peakValue = nan;
    maxPeakIndex = nan;
%     figure(100+l); clf; hold on
%     clear a
    if lap(l).status == "ditch"
        idx = daq.t >= lap(l).t_land + timerange(1) & daq.t <= lap(l).t_land + timerange(2);
%         plot(daq.t(idx)-lap(l).t_land, daq.filt.loadcell(:,idx));
%         ylabel('Force (N)')
        [peakValues, peakIndices] = findpeaks(daq.filt.loadcell(2,idx), daq.t(idx)-lap(l).t_land,'MinPeakDistance',0.05, 'MinPeakProminence', 5,'Annotate','extents','WidthReference','halfheight');
%         xlim(timerange)
        if ~isempty(peakValues)
            [peakValue, maxPeakIndex] = max(peakValues);
            maxPeakIndex = peakIndices(maxPeakIndex);
            lap(l).t_max_grf_takeoff = lap(l).t_land + maxPeakIndex;
        end
    else
        idx = daq.t >= lap(l).t_jump_exact + timerange(1) & daq.t <= lap(l).t_jump_exact + timerange(2);
%         plot(daq.t(idx)-lap(l).t_jump_exact, daq.filt.loadcell(:,idx));
%         ylabel('Force (N)')
        [peakValues, peakIndices] = findpeaks(daq.filt.loadcell(1,idx), daq.t(idx)-lap(l).t_jump_exact,'MinPeakDistance',0.05, 'MinPeakProminence', 1,'Annotate','extents','WidthReference','halfheight');
        if isempty(peakValues)
            [peakValues, peakIndices] = findpeaks(daq.filt.loadcell(3,idx), daq.t(idx)-lap(l).t_jump_exact,'MinPeakDistance',0.05, 'MinPeakProminence', 1,'Annotate','extents','WidthReference','halfheight');
        end
        
%         xlim(timerange)
        if ~isempty(peakValues)
            [peakValue, maxPeakIndex] = max(peakValues);
            maxPeakIndex = peakIndices(maxPeakIndex);
            lap(l).t_max_grf_takeoff = lap(l).t_jump_exact + maxPeakIndex;
        end
    end
%     plot(maxPeakIndex, peakValue, 'or')
end
%% Landing
close all
timerange = [-0.05 0.4]; % [-0.5 0.5] 0.5 sec before to 0.5 sec after
for l = 1:numel(lap)
    lap(l).t_max_grf_landing = nan;
    peakValue = nan;
    maxPeakIndex = nan;
%     figure(100+l); clf; hold on
    if lap(l).status == "ditch"
        idx = daq.t >= lap(l).t_jump_exact + timerange(1) & daq.t <= lap(l).t_jump_exact + timerange(2);
%         plot(daq.t(idx)-lap(l).t_jump_exact, daq.filt.loadcell(:,idx));
%         ylabel('Force (N)')
        [peakValues, peakIndices] = findpeaks(daq.filt.loadcell(2,idx), daq.t(idx)-lap(l).t_jump_exact,'MinPeakDistance',0.05, 'MinPeakProminence', 5,'Annotate','extents','WidthReference','halfheight');
%         xlim(timerange)
        if ~isempty(peakValues)
            [peakValue, maxPeakIndex] = max(peakValues);
            maxPeakIndex = peakIndices(maxPeakIndex);
            lap(l).t_max_grf_landing = lap(l).t_jump_exact + maxPeakIndex;
        end
    else
        idx = daq.t >= lap(l).t_land + timerange(1) & daq.t <= lap(l).t_land + timerange(2);
%         plot(daq.t(idx)-lap(l).t_land, daq.filt.loadcell(:,idx));
%         ylabel('Force (N)')
        [peakValues, peakIndices] = findpeaks(daq.filt.loadcell(1,idx), daq.t(idx)-lap(l).t_land,'MinPeakDistance',0.05, 'MinPeakProminence', 1,'Annotate','extents','WidthReference','halfheight');
        if isempty(peakValues)
            [peakValues, peakIndices] = findpeaks(daq.filt.loadcell(3,idx), daq.t(idx)-lap(l).t_land,'MinPeakDistance',0.05, 'MinPeakProminence', 1,'Annotate','extents','WidthReference','halfheight');
        end
        
%         xlim(timerange)
        if ~isempty(peakValues)
            [peakValue, maxPeakIndex] = max(peakValues);
            maxPeakIndex = peakIndices(maxPeakIndex);
            lap(l).t_max_grf_landing = lap(l).t_land + maxPeakIndex;
        end
    end
    %plot(maxPeakIndex, peakValue, 'or')
end
%%
phases1 = interp1(csc.t, csc.phase, [lap([lap.status]=="jump").t_jump_exact]);
phases2 = interp1(csc.t, csc.phase, [lap([lap.status]=="jump").t_max_grf_takeoff]);
phases3 = interp1(csc.t, csc.phase, [lap([lap.status]=="jump").t_max_grf_landing]);
phases4 = interp1(csc.t, csc.phase, [lap([lap.status]=="ditch").t_jump_exact]);
phases5 = interp1(csc.t, csc.phase, [lap([lap.status]=="ditch").t_max_grf_takeoff]);
phases6 = interp1(csc.t, csc.phase, [lap([lap.status]=="ditch").t_max_grf_landing]);
%% phase histogram
clear a
% binning
hist.phase_bin_size = 5; % deg
hist.phase_edges_deg = -180:hist.phase_bin_size:180;

figure(1001); clf
a(1) = subplot(2,3,1);
histogram(phases1, hist.phase_edges_deg, 'FaceColor', 'g','FaceAlpha',1 ,'Normalization', 'probability');
xlim([-180 180])
xlabel('Phase (deg)')
ylabel('Normalized count')
title('phase at the time of jump take-off')
a(2) = subplot(2,3,2);
histogram(phases2, hist.phase_edges_deg, 'FaceColor', 'g','FaceAlpha',1 ,'Normalization', 'probability');
xlim([-180 180])
xlabel('Phase (deg)')
ylabel('Normalized count')
title('phase at the time of peak GRF (jump take-off)')
a(2) = subplot(2,3,3);
histogram(phases3, hist.phase_edges_deg, 'FaceColor', 'g','FaceAlpha',1 ,'Normalization', 'probability');
xlim([-180 180])
xlabel('Phase (deg)')
ylabel('Normalized count')
title('phase at the time of peak GRF (jump landing)')
linkaxes(a, 'xy')

a(1) = subplot(2,3,4);
histogram(phases4, hist.phase_edges_deg, 'FaceColor', 'g','FaceAlpha',1 ,'Normalization', 'probability');
xlim([-180 180])
xlabel('Phase (deg)')
ylabel('Normalized count')
title('phase at the time of ditch take-off')
a(2) = subplot(2,3,5);
histogram(phases5, hist.phase_edges_deg, 'FaceColor', 'g','FaceAlpha',1 ,'Normalization', 'probability');
xlim([-180 180])
xlabel('Phase (deg)')
ylabel('Normalized count')
title('phase at the time of peak GRF (ditch take-off)')
a(3) = subplot(2,3,6);
histogram(phases5, hist.phase_edges_deg, 'FaceColor', 'g','FaceAlpha',1 ,'Normalization', 'probability');
xlim([-180 180])
xlabel('Phase (deg)')
ylabel('Normalized count')
title('phase at the time of peak GRF (ditch landing)')
linkaxes(a, 'xy')

set(gcf, 'Position', [100 100 1600 1100]);
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'phase_grf.jpg'])

figure(1002); clf
a(1) = subplot(2,3,1);
polarplot(phases1, ones(size(phases1)),'og'); hold on
m = mean(exp(1i*phases1),'omitnan');
polarplot(angle(m), abs(m),'or', 'MarkerFaceColor', 'r')
title('phase at the time of jump take-off')
a(2) = subplot(2,3,2);
polarplot(phases2, ones(size(phases2)),'og'); hold on
m = mean(exp(1i*phases2),'omitnan');
polarplot(angle(m), abs(m),'or', 'MarkerFaceColor', 'r')
title('phase at the time of peak GRF (jump take-off)')
a(3) = subplot(2,3,3);
polarplot(phases3, ones(size(phases3)),'og'); hold on
m = mean(exp(1i*phases3),'omitnan');
polarplot(angle(m), abs(m),'or', 'MarkerFaceColor', 'r')
title('phase at the time of peak GRF (jump landing)')

a(1) = subplot(2,3,4);
polarplot(phases4, ones(size(phases4)),'og'); hold on
m = mean(exp(1i*phases4),'omitnan');
polarplot(angle(m), abs(m),'or', 'MarkerFaceColor', 'r')
title('phase at the time of ditch take-off')
a(2) = subplot(2,3,5);
polarplot(phases5, ones(size(phases5)),'og'); hold on
m = mean(exp(1i*phases5),'omitnan');
polarplot(angle(m), abs(m),'or', 'MarkerFaceColor', 'r')
title('phase at the time of peak GRF (ditch take-off)')
a(3) = subplot(2,3,6);
polarplot(phases6, ones(size(phases6)),'og'); hold on
m = mean(exp(1i*phases6),'omitnan');
polarplot(angle(m), abs(m),'or', 'MarkerFaceColor', 'r')
title('phase at the time of peak GRF (ditch landing)')

set(gcf, 'Position', [100 100 1600 1100]);
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'phase_grf_polar.jpg'])
figure(1003); clf;
subplot(2,1,1)
plot([lap([lap.status]=="jump").no],[lap([lap.status]=="jump").t_max_grf_takeoff]-[lap([lap.status]=="jump").t_jump_exact],'ob'); hold on;
plot([lap([lap.status]=="ditch").no],[lap([lap.status]=="ditch").t_max_grf_takeoff]-[lap([lap.status]=="ditch").t_jump_exact],'or')
ylabel('Time difference between take-off and max GRF (sec)')
ylim([0 2])
subplot(2,1,2)
plot([lap([lap.status]=="jump").no],[lap([lap.status]=="jump").t_max_grf_landing]-[lap([lap.status]=="jump").t_land],'ob'); hold on;
plot([lap([lap.status]=="ditch").no],[lap([lap.status]=="ditch").t_max_grf_landing]-[lap([lap.status]=="ditch").t_jump_exact],'or')
ylabel('Time difference between landing and max GRF (sec)')
xlabel('Lap number')
ylim([0 2])
set(gcf, 'Position', [100 100 1600 1100]);
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'grf_time.jpg'])