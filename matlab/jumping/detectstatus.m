% clc; clear; close all;
% % exp_directory = 'E:\Rat1068\Analysis\2022-12-20\';
% % exp_directory = 'E:\Rat1055\Analysis\2022-11-09\';
% exp_directory = 'E:\Rat980\Analysis\2021-12-21\';
% load(fullfile(exp_directory,'processed_data.mat'), 'daq', 'lap', 'cluster', 'colors','hist','exp', 'field');
%%
% for f= 1:numel(field)
%     clear a
%     figure(f)
%     i = 0;
%     xmax = 1.1 * exp.xmax;
%     for x_reference = ["stationary" "moving"]
%         i = i + 1;
%         a(i) = subplot(2,2,i);
%         c = field(f).cluster;
%         firingRates = mean([cluster(c).hist.(x_reference).left_jump; ...
%             cluster(c).hist.(x_reference).right_jump; ...
%             cluster(c).hist.(x_reference).left_ditch; ...
%             cluster(c).hist.(x_reference).right_ditch]);
%         histogram('BinCounts', firingRates, 'BinEdges', hist.edges, 'FaceColor',lap(l).cross_color + "o",'FaceAlpha',0.1); hold on; %rate map histogram (90% transparent)
%         firingRates = mean([field(f).(x_reference).hist.left_jump; ...
%             field(f).(x_reference).hist.right_jump; ...
%             field(f).(x_reference).hist.left_ditch; ...
%             field(f).(x_reference).hist.right_ditch]);
%         histogram('BinCounts', firingRates, 'BinEdges', hist.edges, 'FaceColor',colors(c),'FaceAlpha',0.8); hold on; %rate map histogram (20% transparent)
%         title(x_reference)
%         xlim([0 xmax])
%         idx = cluster(c).x >= field(f).(x_reference).range(1) & cluster(c).x <= field(f).(x_reference).range(2);
%         b(i) = subplot(2,2,i+2);
%         
%         data = [];
%         
%         for j = 1:numel(firingRates)
%             for k = 1:round(4*firingRates(j)) % 4x precesion
%                 data = [data; hist.centers(j)];
%             end
%         end
%         p.(x_reference) = fitdist(data,'Normal');
%         %p.(x_reference) = pd;
%         plot(0:0.1:xmax, normpdf(0:0.1:xmax,p.(x_reference).mu,p.(x_reference).sigma));
%         xlim([0 xmax])
%         % Calculate the probability density function (PDF) of the firing rates
%         pdf = firingRates / sum(firingRates * hist.bin_size);
%         % Calculate the entropy of the PDF
%         entropy = -sum(pdf .* log2(pdf+eps));
%         % Calculate the spatial information score
%         skaggsScore.(x_reference) = entropy / log2(length(hist.centers));
% 
%     end
% 
%     if field(f).moving.center < lap(1).gap(2)+lap(1).corr - 10
%         red_black_stationary{1} = 'red';
%         red_black_moving{1} = 'black';
%     else
%         red_black_stationary{1} = 'black';
%         red_black_moving{1} = 'red';
%     end
%     
%     
%     if skaggsScore.stationary > skaggsScore.moving
%         red_black_stationary{2} = 'red';
%         red_black_moving{2} = 'black';
%     else
%         red_black_stationary{2} = 'black';
%         red_black_moving{2} = 'red';
%     end
%     
%     
%     if p.stationary.sigma < p.moving.sigma
%         red_black_stationary{3} = 'red';
%         red_black_moving{3} = 'black';
%     else
%         red_black_stationary{3} = 'black';
%         red_black_moving{3} = 'red';
%     end
%     
%     subplot(2,2,3);
% %     title(['Stationary: {\color{' red_black_stationary{1} '}location \color{' red_black_stationary{2} '}Skaggs ' num2str(skaggsScore.stationary) '\color{' red_black_stationary{3} '}std dev ' num2str(p.stationary.sigma) '}'])
%     title(['Stationary: {\color{' red_black_stationary{1} '}location \color{' red_black_stationary{2} '}Skaggs \color{' red_black_stationary{3} '}std dev}'])
%     
% subplot(2,2,4);
% %     title(['Moving: {\color{' red_black_moving{1} '}location \color{' red_black_moving{2} '}Skaggs ' num2str(skaggsScore.moving) '\color{' red_black_moving{3} '}std dev ' num2str(p.moving.sigma) '}'])
%     title(['Moving: {\color{' red_black_moving{1} '}location \color{' red_black_moving{2} '}Skaggs \color{' red_black_moving{3} '}std dev}'])
%     
%     linkaxes(a,'xy')
%     linkaxes(b,'xy')
%     set(gcf, 'Position', [100 100 1600 1100]);
%     saveas(gcf,[exp_directory filesep 'Analysis' filesep 'status_field' num2str(f) '.jpg'])
% end


%%
for f = 1:numel(field)
    figure(400+f); clf; hold on;
    c = field(f).cluster;
    dir = field(f).dir;
    field(f).lap_mean = nan(numel(lap),1);
    for l = [lap([lap.dir] == dir).no]
        plot(lap(l).gap(1), l, lap(l).cross_color + "|");
        plot(lap(l).gap(2), l, lap(l).cross_color + "|");
        idx = cluster(c).lap == l & cluster(c).x >= field(f).range(1) & cluster(c).x <= field(f).range(2);
        field(f).lap_mean(l) = mean(cluster(c).x(idx));
%         cc = cluster(c).x(idx);
%         if ~isempty(cc)
%             field(f).lap_mean(l) = mean([cc(1), cc(end)]);
%         end
        plot(cluster(c).x(idx), l*ones(size(cluster(c).x(idx))) ,'.','Color', colors(c));
    end
    plot([field(f).lap_mean], [lap.no] ,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c));
%     field(f).err.stationary = mean(abs([field(f).lap_mean] - mean([field(f).lap_mean], 'omitnan')), 'omitnan');
    field(f).err.stationary = std([field(f).lap_mean], 'omitnan');
%         field(f).err.stationary = std(x_all, 'omitnan');
%     field(f).err.moving = mean(abs([field(f).lap_mean] - [lap.gap_length]' - mean([field(f).lap_mean] - [lap.gap_length]', 'omitnan')), 'omitnan');
    field(f).err.moving = std([field(f).lap_mean] - [lap.gap_length]', 'omitnan');
    disp(field(f).err.stationary);
    disp(field(f).err.moving);
    stability = sum(~isnan([field(f).lap_mean]))/numel([lap([lap.dir] == dir).no]);
    ratio = field(f).err.stationary / (field(f).err.stationary + field(f).err.moving + eps);
    if stability < 0.3
        if field(f).err.stationary < field(f).err.moving
            disp('not enough laps')
%             field(f).x_ref = "stationary";
%             title(['{\color{green} ' convertStringsToChars(dir) '}, not enough laps: {\color{red} stationary} - stability: ' num2str(stability) ' - ratio: ' num2str(ratio)])
        else
            disp('not enough laps')
%             field(f).x_ref = "moving";
%             title(['{\color{green} ' convertStringsToChars(dir) '}, not enough laps: {\color{red} moving} - stability: ' num2str(stability) ' - ratio: ' num2str(ratio)])
        end
        field(f).x_ref = "ambiguous"; 
                title(['{\color{green} ' convertStringsToChars(dir) '}, not enough laps: {\color{red} ambiguous} - stability: ' num2str(stability) ' - ratio: ' num2str(ratio)])

    elseif ratio < 0.45
        disp('stationary');
        field(f).x_ref = "stationary";
        title(['{\color{green} ' convertStringsToChars(dir) '}, {\color{red} stationary} - stability: ' num2str(stability) ' - ratio: ' num2str(ratio)])
    elseif ratio > 0.55
        disp('moving');
        field(f).x_ref = "moving";
        title(['{\color{green} ' convertStringsToChars(dir) '}, {\color{red} moving} - stability: ' num2str(stability) ' - ratio: ' num2str(ratio)])
    else
        field(f).x_ref = "ambiguous";        
        title(['{\color{green} ' convertStringsToChars(dir) '}, {\color{red} ambiguous} - stability: ' num2str(stability) ' - ratio: ' num2str(ratio)])

    end
    field(f).stability = stability;
    xlim([0 exp.xmax])
    set(gcf, 'Position', [100 100 1600 1100]);
    saveas(gcf,[exp_directory filesep 'Analysis' filesep 'status_field' num2str(f) '.jpg'])
end
