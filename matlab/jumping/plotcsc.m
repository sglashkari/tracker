%%PLOTCSC plots all CSC channels for comparison
%
%   See also ANALYZEDATA.
%
%   SGL 2022-01-25 (originally 2021-01-31)
%
clc; clear; close all

% [csc_filename,csc_directory] = uigetfile('D:\Analysis\*.ncs','Select One or More CSC Files','MultiSelect','on');
% if isa(csc_filename,'double')
%     return;
% elseif isa(csc_filename,'char')
%     csc_filename = {csc_filename};
% end
exp_directory = 'D:\Analysis\2021-12-10';
csc_filename = cell(4,1);
for sh=1:4
    csc_filename{sh}=fullfile(exp_directory, ['CA1-Shank' num2str(sh)],'B','LFP32.ncs'); % use plotcsc to optimize it
end

for ch=1:length(csc_filename)
    [timecsc,lfp] = read_bin_csc(csc_filename{ch});
    idx = timecsc >= 100 & timecsc < 200;
    timecsc = timecsc(idx);
    lfp = lfp(idx);
    
    [theta, phase, mag] = filtertheta(timecsc,lfp);
    
    ax1 = subplot(3,1,1);
    %plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
    hold on;
    plot(timecsc,theta *1e6);%,'r');
    ylim([-300 300])
    ylabel('Theta (\muV)')
    %title(['lap ' num2str(l)])
    
    ax2 = subplot(3,1,2);
    hold on
    plot(timecsc,phase);%,'b');
    ylim([-200 200])
    ylabel('Phase (Angle)')
    
    ax3 = subplot(3,1,3);
    hold on
    plot(timecsc,mag*1e6);%,'b');
    %ylim([-200 200])
    ylabel('Magnitude')
    
    
    disp(['Average magnitude of CSC' num2str(ch) ' is ' num2str(mean(mag) *1e6,'%.2f')]);
end
legend(num2str(1:length(csc_filename)))
linkaxes([ax1 ax2 ax3],'x')

%%
% figure(2)
% cscs = 11:12;
% for ch=cscs
%     
%     csc_filename= fullfile(exp_directory,'Neuralynx',['CSC' num2str(ch) '.ncs']);
%     [timecsc,lfp] = readcsc(csc_filename, lap(l).t * 1e6); % microseconds
%     [theta, phase, mag] = filtertheta(timecsc,lfp);
%     
%     ax1 = subplot(3,1,1);
%     %plot(timecsc,lfp * 1e6,'Color','#D0D0D0')
%     hold on;
%     plot(timecsc,theta *1e6);%,'r');
%     ylim([-300 300])
%     ylabel('Theta (\muV)')
%     title(['lap ' num2str(l)])
%     
%     ax2 = subplot(3,1,2);
%     hold on
%     plot(timecsc,phase);%,'b');
%     ylim([-200 200])
%     ylabel('Phase (Angle)')
%     
%     ax3 = subplot(3,1,3);
%     hold on
%     plot(timecsc,mag*1e6);%,'b');
%     %ylim([-200 200])
%     ylabel('Magnitude')
%     
% end
% legend(num2str(cscs'))
% linkaxes([ax1 ax2 ax3],'x')