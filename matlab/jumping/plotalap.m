clc; clear; close all

exp_directory = 'D:\Analysis\2021-12-21';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'pos','posi', 'lap', 'cluster','ppcm', 'colors','xmax','hist','daq');


cluster_no = 29; %[17 18 23 28 31 39]; %[7 28 29 40]; %24 ; %[23 28]; %[17 18 23 28]; % [22 29]; [3 13 17 18 22 23 28 29 31]; [cluster.no];
colors(18)="#f58231";
colors(28)="#3cb44b";

v_thresh = 0;
direction = 'left';
showCSC = true;
timerange = [-2 2]; % 2 sec before to 2 sec after
x_reference = 'moving';
t_reference = 'take-off';
%l = 74; % 76; 90 10 22 (( 88, 90, 109, 119, ((( 
l= 60; %183;
lap_no = 62;

%%
legendCell = strings(size(cluster))';
for c=[cluster.no]
    if ismember(c,cluster_no)
        legendCell(c) = string(num2str(c, 'cluster #%-d'));
    else
        legendCell(c) = "";
    end
end

% time of interest: landing or take-off
for l = 1:length(lap)
    lap(l).t_land = lap(l).t_cross(end);
    if strcmp(x_reference,'stationary')
        lap(l).dx_interest = 0;
    elseif strcmp(x_reference,'moving')
        lap(l).dx_interest = lap(l).corr;
    end
    if strcmp(t_reference,'take-off')
        lap(l).t_interest = lap(l).t_jump_exact;
    elseif strcmp(t_reference,'landing')
        lap(l).t_interest = lap(l).t_land;
    end
end

depth = lap(l).gap_depth;
no_states = length(unique([lap.status]));

%% X-Z plot
figure(1); clf;
l = lap_no;
for c = cluster_no

    hold on
    idx = [cluster(c).lap]==l;
    dx = [lap(cluster(c).lap(idx)).dx_interest];
    dx = reshape(dx,[],1);
    plot(cluster(c).p(idx,1)+dx+rand(sum(idx),1)/10,cluster(c).p(idx,3)+rand(sum(idx),1)/10,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
    
%     idx = pos.lap == l;
%     dx = [lap(pos.lap(idx)).dx_interest];
%     dx = reshape(dx,[],1);
%     h = plot(pos.p(idx,1)+dx,pos.p(idx,3),'.', 'MarkerSize',0.2);
%     set(h, 'Color', '#D0D0D0');
    axis equal
    title(['Lap No.' num2str(l) ', place field(s), ' convertStringsToChars(lap(l).dir) 'ward ' convertStringsToChars(lap(l).status)])
    xlabel('Horizontal position (cm)')
    ylabel('Elevation (cm)')
    
    ylim([-35 20])
    xlim([-50 20+max([lap.gap_length])]+max([lap.dx_interest]))
    
    % gap range
    lim = xlim;
    plot([lim(1) lap(l).dx_interest], [0 0],lap(l).cross_color,'LineWidth',0.25);
    plot([lap(l).dx_interest lap(l).dx_interest], [-depth 0],lap(l).cross_color,'LineWidth',0.25);
    plot(repmat(lap(l).gap_length+lap(l).dx_interest,1,2), [-depth 0],lap(l).cross_color,'LineWidth',0.25);
    plot([lap(l).dx_interest lap(l).gap_length+lap(l).dx_interest], [-depth -depth],lap(l).cross_color,'LineWidth',0.25);
    plot([lap(l).gap_length+lap(l).dx_interest lim(2)+lap(l).dx_interest], [0 0],lap(l).cross_color,'LineWidth',0.25);
    %Reverse the stacking order so that the patch overlays the line
    chi=get(gca, 'Children');
    set(gca, 'Children',flipud(chi));
    
end
%toc

set(gcf, 'Position', [0 0 800 200*no_states]);
saveas(gcf,fullfile(exp_directory, 'Analysis',['Z-X-lap' mat2str(l) '-cl' mat2str(cluster_no) '_' x_reference '.png']))
%sgtitle(['Elevation versus Horizontal position, ' x_reference ' frame of ref (velocity filtered: v >= ' num2str(v_thresh) ' cm/s)']);

%% 3D plot
f=figure(2); clf;
% f.WindowState = 'maximized';
% table surface
[x,y] = meshgrid(-150:150,-7:7);
z = - depth * (x > 0) + depth * (x > max([lap.gap_length]+[lap.dx_interest])); %
z(x==1)=nan;
z(x==floor(max([lap.gap_length]+[lap.dx_interest])))=nan;
hold on
surf(x,y,z,'FaceAlpha',0.25,'EdgeAlpha',0.1,'FaceColor','#D0D0D0'); hold on

axis equal
title(['Lap No.' num2str(l) ', place field(s), ' convertStringsToChars(lap(l).dir) 'ward ' convertStringsToChars(lap(l).status)])

xlabel('X position (cm)')
ylabel('Y Position (cm)')
zlabel('Z position (cm)')

for c = cluster_no
    idx = [cluster(c).lap]==l;
    dx = reshape([lap(cluster(c).lap(idx)).dx_interest],[],1);
    
    for nr = 1:nnz(idx) % plotting step by step
        iidx = find(idx,nr);
        iidx = iidx(end);
        s = scatter3(cluster(c).p(iidx,1)+dx(nr)+rand/10,cluster(c).p(iidx,2)+rand/10,cluster(c).p(iidx,3)+rand/10,'o','MarkerEdgeColor','black', 'MarkerFaceColor', colors(c)); hold on
        alpha(s,0.5)
        
        %     idx = pos.lap == l;
        %     dx = [lap(pos.lap(idx)).dx_interest];
        %     dx = reshape(dx,[],1);
        %     h = plot(pos.p(idx,1)+dx,pos.p(idx,3),'.', 'MarkerSize',0.2);
        %     set(h, 'Color', '#D0D0D0');
        
        
        ylim([-35 20])
        view([0.2 -1 0.5])
        zlim([-35 20])
        xlim([-50 50+max([lap.gap_length])])
        
        if nr < nnz(idx)
            iiidx = find(idx,nr+1);
            iiidx = iiidx(end);
            %             pause(cluster(c).p(iiidx)-cluster(c).p(iidx))
            pause(0.1)
        end
    end
    
end
%toc

%set(gcf, 'Position', [0 0 800 200*no_states]);
        ylim([-35 20])
        %view([0.2 -1 0.5])
        zlim([-35 20])
        xlim([-50 50+max([lap.gap_length])])
saveas(gcf,fullfile(exp_directory, 'Analysis',['Z-X-lap' mat2str(l) '-cl' mat2str(cluster_no) '_' x_reference '.png']))
%sgtitle(['Elevation versus Horizontal position, ' x_reference ' frame of ref (velocity filtered: v >= ' num2str(v_thresh) ' cm/s)']);