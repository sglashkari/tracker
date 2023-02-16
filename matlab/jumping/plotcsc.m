%%PLOTCSC plots all CSC channels for comparison
%
%   See also ANALYZEDATA, APBIN2LFP.
%
%   Date 2023-01-01 (originally 2021-01-31, 2022-01-25)
%   
%   Author Shahin G Lashkari
%
clc; clear; close all;
answer = inputdlg({'Rat', 'Date'},'Enter the rat number and the date',[1 30],{'1068', '2022-12-13'});
rat_no = answer{1};
date_str = answer{2};
%% Selecting the appropriate files
[lfp_filename,lfp_directory] = uigetfile(['E:\Rat' rat_no '\Analysis\' date_str '\LFP\*.ncs'],'Select One or More CSC Files','MultiSelect','on');
if isa(lfp_filename,'double')
    return;
end

tic
m = zeros(size(lfp_filename));
for ch=1:length(lfp_filename)
    [timecsc,lfp] = read_bin_csc(fullfile(lfp_directory, lfp_filename{ch}));
    idx = timecsc >= 100 & timecsc < 200;
    timecsc = timecsc(idx);
    lfp = lfp(idx);
    
    [theta, phase, mag] = filterlfp(timecsc, lfp, 'theta');
    
    fprintf('%.f%% Average magnitude of CSC %d is %.2f\n', 100*ch/length(lfp_filename), ch, mean(mag)*1e6);
    m(ch) = mean(mag) *1e6;
end
[sorted_m, idx] = sort(m, 'descend');
figure(1); clf
plot(sorted_m); hold on
plot(idx-1,'.')
toc
disp('Highest amplitude channels are (in order):')
fprintf('%d\n', idx(1:20)-1)
disp('')
plot(100*sorted_m/max(m));
set(gcf, 'Position', [100 100 1600 1100]);
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'theta_comparison.fig'])
% 
%%
figure(2); clf
tic
cscs = idx(1:50:end);
for ch=cscs
    
    [timecsc,lfp] = read_bin_csc(fullfile(lfp_directory, lfp_filename{ch}));
    idx = timecsc >= 100 & timecsc < 200;
    timecsc = timecsc(idx);
    lfp = lfp(idx);
    [theta, phase, mag] = filterlfp(timecsc, lfp, 'theta');
    
    ax1 = subplot(2,1,1);
    %plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc,theta *1e6);%,'r');
%     ylim([-300 300])
    ylabel('Theta (\muV)')
    %title(['lap ' num2str(l)])

    ax2 = subplot(2,1,2);
    hold on
    plot(timecsc,mag*1e6);%,'b');
    %ylim([-200 200])
    ylabel('Magnitude')
    
end
legend(num2str(cscs'-1))
linkaxes([ax1 ax2],'x')
toc