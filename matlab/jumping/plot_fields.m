function plot_fields(exp_directory)
%% plot fields based on location
%
%   See also EXTRACTFIELDS, PLOTRATEMAP, PLOTTIME.
%
%   Date 2023-01-13
%   Author Shahin G Lashkari
%
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
    load(mat_filename, 'field', 'exp', 'lap');
else
    mat_filename = fullfile(exp_directory, 'processed_data.mat');
    load(mat_filename, 'field', 'exp', 'lap');
end
%%
figure(500); clf
i = 0;
clear a
for dir = ["left" "right"]
    i = i + 1;
%     idx = [field.dir] == dir;
    a(i) = subplot(1,2,i); hold on
    
    xlim([0 exp.xmax]);
    ylim([-0.01 1.01])
    
    for l=[lap.no]  % only the desired direction and state
        % gap range
        plot(repmat(lap(l).gap(1),1,2), ylim,'k','LineWidth',0.1);
        plot(repmat(lap(l).gap(2),1,2), ylim,'k','LineWidth',0.1);
    end
    
    for idx = find([field.dir] == dir)
        if field(idx).x_ref == "moving"
            color = 'magenta';
        elseif field(idx).x_ref == "stationary"
            color = 'green';
        else
            color = 'white';
        end
        scatter([field(idx).center], [field(idx).ratio_jump], 100*(log([field(idx).peak])+1), 'MarkerEdgeColor','k', 'MarkerFaceColor', color, 'MarkerFaceAlpha', min(1, 2*field(idx).stability),'MarkerEdgeAlpha',min(1, 2*field(idx).stability)); hold on;
        text([field(idx).center], [field(idx).ratio_jump], string([field(idx).no]), 'VerticalAlignment','middle','HorizontalAlignment','center', 'FontSize', 2*(log([field(idx).peak])+1)+4 ,'Color' , 'black')
    end
    title(['field direction: ' convertStringsToChars(dir)], 'Interpreter', 'none');
    ylabel('Jumping field probability')
    xlabel('Horizontal position (cm)')
    box on
end
set(gcf, 'Position', [100 100 1600 1100]);
sgtitle({['Rat No. ' num2str(exp.rat_no) ];['Day: ' char(exp.date)]});
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'fields.jpg'])
saveas(gcf,[exp_directory filesep 'Analysis' filesep 'fields.fig'])

end