% want to plot seperate clusters but combine all laps each maze with
% heatmat , on Sept 24,2020 updated to no velocity, on Sept 30,2020 updated
% to seconds per bin for occupancy and heat map
clc; close all;
if ~exist('pos_clust','var')
    loaddata;
end

%% selecting a portion of the cells and plotting
close all


%%
% ...
%
%
Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
listing = dir(fullfile(Nlx_directory,'**','TT*.ntt'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing);
absolue_paths = folders+'\'+names;

%% maze
for maze = 1:5
    
    TimeEV = readevent;
    StartTime = TimeEV(1:2:end);
    EndTime = TimeEV(2:2:end);
    
    %% Pos.p
    [TimePos,x,y] = readposp(fullfile(Nlx_directory,'pos.p'), StartTime(maze), EndTime(maze));
    y = 480 - y;
    AngularPosition = rad2deg(atan2(y-240, x-320));
    AngularPosition = wrapTo360(AngularPosition);
    
    %% lap detector
    lap = 3;
    
    TimePosF = TimePos(~isnan(AngularPosition));
    AngularPositionF = AngularPosition(~isnan(AngularPosition));
    AngularPositionF = AngularPositionF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));
    TimePosF =  TimePosF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));
    
    LapTime = TimePosF(diff(AngularPositionF) < -340);
    TimeRange = [LapTime(lap) LapTime(lap+1)]; % seconds
    
    
    %% selecting a portion of the cells and plotting
    
    start_t = seconds(sort(unique(start_time))*1e-6)';
    end_t = seconds(sort(unique(end_time))*1e-6)';
    end_t(3) = [];
    
    S = timerange(start_t(maze),end_t(maze));
    pos_col = [false(2,1); true];
    clust_col = (maze_no == maze) & (tt_no < 18);
    maze1 = pos_clust(S, [pos_col; clust_col]);
    maze_table = timetable2table(maze1);
    maze_table.time = seconds(maze1.time);
    maze_array = table2array(maze_table);
    no_cells = size(maze1,2)-1;
    
    all_cells = repmat(1:no_cells,size(maze1,1),1);
    modified_framerate = maze_array(:,3:end);
    g = modified_framerate.*all_cells;
    h=~isnan(g);
    c1 = find(h);
    [r,c] = find(h);
    
    array_combined = c; %g(c2)
    
    array_combined(array_combined==0) = NaN;
    array_combined(1)=0;
    array_combined(2)=no_cells+1;
    
    t = seconds(maze_array(r,1));
    
    %% selecting a portion
    
    
    %% Histogram ..
    
    
    frameRate = 30; %frame rate in Hz
    edges = linspace(0, 360, 360);
    Edges = movmean(edges,2);
    t_spike = seconds(t);    
    
    %% interpolation
    modified_framerate = 3000; % interpolation rate 3 kHz
    TimePosInterp = TimePos(1):1/modified_framerate: TimePos(end);
    
    xInterp = interp1(TimePos,x,TimePosInterp);
    yInterp = interp1(TimePos,y,TimePosInterp);
    AngularPositionInterp = rad2deg(atan2(yInterp-240, xInterp-320));
    AngularPositionInterp = wrapTo360(AngularPositionInterp);
    
    %% histogram
    edges = linspace(0, 360, 360);
    HistOcc = histcounts(AngularPosition, edges)/frameRate;
    HistOccInterp = histcounts(AngularPositionInterp, edges)/modified_framerate;
    
    %% velocity
    % 5 ft = 1524 mm = 480 pixels
    % each pixel = 3.175 mm
    
    vx = gradient(x)./gradient(TimePos); % pixels/sec
    vy = gradient(y)./gradient(TimePos); % pixels/sec
    
    velocity = sqrt(vx.^2+vy.^2); % pixels/sec
    velocity = velocity * 0.3175; % cm/sec
    
    velocityF = gradient(AngularPositionF)./gradient(TimePosF); % deg/sec (1 deg ~ 1.3 cm)
    
    %% rate map of cells (histogram of occupancy corrected, and velocity)
    E=zeros(17,1);
    RateMap = zeros(17,length(edges)-1);
    for i=1:17
        t_spike_cluster = t_spike(array_combined==i);
        SpikePos = interp1(TimePosF,AngularPositionF,t_spike_cluster);
        SpikeVelocity = interp1(TimePos,velocity,t_spike_cluster);
        
        HistClust = histcounts(SpikePos(SpikeVelocity>-1), edges); % speed filter (no speed filter)
        
        RateMap(i,:) = (HistClust./HistOccInterp);
        
        [~, argmaxRateMap] = max(RateMap(i,:));
        E(i)=Edges(argmaxRateMap);
    end
    
    % no sorting
    
    for i=[2 3 9 12]
        fig = figure(maze*100+i);
        subplot(4,1,1);
        imagesc(HistOccInterp)
        colormap('parula')
        cb=colorbar;
        cb.Position = cb.Position + [cb.Position(1)*0.12 0 -cb.Position(3)*0.3 0];
        title(['Occupancy (Interpolated - cluster no. ' num2str(i) '- maze ' num2str(maze)])
        
        
        subplot(4,1,2);
        bar(Edges(2:end),RateMap(i,:));
        ylim([0 18])
        xlim([0 360])
        title(['Rate Map (Interpolated - Maze ' num2str(maze)])
        
        subplot(4,1,3);
        imagesc(RateMap(i,:))
        colormap('parula')
        cb=colorbar;
        cb.Position = cb.Position + [cb.Position(1)*0.12 0 -cb.Position(3)*0.3 0];
        xlim([0 360])
        title(['Rate Map (Interpolated - Maze ' num2str(maze)])
        
        subplot(4,1,4)
        depth = [-1.5 -1.5; -3 -1.5; -1.5 -3 ; -3 -3 ; -3 -3];
        if maze ~= 5
            x = [0  89  89  112 112 242 242 269 269 360];
        else
            x = [0  90  90  118 118 237 237 269 269 360];
        end
        y = [0  0   depth(maze,1) depth(maze,1) 0   0   depth(maze,2) depth(maze,2)  0   0];
        line(x,y,'LineWidth',2)
        ylim([-4 1])
        xlim([0 360])
        
        % save as pdf
        set(fig,'PaperOrientation','landscape','PaperSize',[8.5 11], 'PaperUnits', 'inches','PaperPosition', [0 0 8.5 11]);
        saveas(fig, ['Cluster' num2str(i) '-maze' num2str(maze) '.pdf'])
        
    end
end