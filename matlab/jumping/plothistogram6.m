% same as 4 
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
counter = 0;
for maze = 1:5
    close all
    TimeEV = readevent;
    StartTime = TimeEV(1:2:end);
    EndTime = TimeEV(2:2:end);
    
    %% Pos.p
    [TimePos,x,y] = readposp(fullfile(Nlx_directory,'pos.p'), StartTime(maze), EndTime(maze));
    y = 480 - y;
    AngularPosition = rad2deg(atan2(y-240, x-320));
    AngularPosition = wrapTo360(AngularPosition);
    
    %% lap detector
    TimePosF = TimePos(~isnan(AngularPosition));
    AngularPositionF = AngularPosition(~isnan(AngularPosition));
    AngularPositionF = AngularPositionF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));
    TimePosF =  TimePosF(TimePosF>StartTime(maze) & TimePosF<EndTime(maze));
    
    LapTime = TimePosF(diff(AngularPositionF) < -340);
    
    if maze == 2
        maxlap = 4;
    else
        maxlap = 10;
    end
    for lap = 1:maxlap
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
        f = maze_array(:,3:end);
        g = f.*all_cells;
        h=~isnan(g);
        c1 = find(h);
        [r,c] = find(h);
        
        array_combined = c; %g(c2)
        
        array_combined(array_combined==0) = NaN;
        array_combined(1)=0;
        array_combined(2)=no_cells+1;
        
        t = seconds(maze_array(r,1));
        
        %% Histogram ..
        
        frameRate = 30; %frame rate in Hz
        edges = linspace(0, 360, 360);
        Edges = movmean(edges,2);
        t_spike = seconds(t);
        HistOcc = histcounts(AngularPositionF, edges);
        
        
        %% interpolation
        f = 3000; % interpolation 3 kHz
        TimePosInterp = TimePos(1):1/f: TimePos(end);
        
        xInterp = interp1(TimePos,x,TimePosInterp);
        yInterp = interp1(TimePos,y,TimePosInterp);
        AngularPositionInterp = rad2deg(atan2(yInterp-240, xInterp-320));
        AngularPositionInterp = wrapTo360(AngularPositionInterp);
        
        %% histogram
        edges = linspace(0, 360, 360);
        HistOcc = histcounts(AngularPosition, edges);
        HistOccInterp = histcounts(AngularPositionInterp, edges);
        
        %% velocity
        % 5 ft = 1524 mm = 480 pixels
        % each pixel = 3.175 mm
        
        vx = gradient(x)./gradient(TimePos); % pixels/sec
        vy = gradient(y)./gradient(TimePos); % pixels/sec
        
        velocity = sqrt(vx.^2+vy.^2); % pixels/sec
        velocity = velocity * 0.3175; % cm/sec
        
        %% rate map of cells (histogram of occupancy corrected, and velocity)
        E=zeros(17,1);
        RateMap = zeros(17,length(edges)-1);
        for i=1:17
            t_spike_cluster = t_spike(array_combined==i);
            t_spike_cluster = t_spike_cluster(t_spike_cluster >= LapTime(lap) & t_spike_cluster <=LapTime(lap+1));
            SpikePos = interp1(TimePosF,AngularPositionF,t_spike_cluster);
            SpikeVelocity = interp1(TimePos,velocity,t_spike_cluster);
            
            HistClust = histcounts(SpikePos(SpikeVelocity>5), edges); % speed filter
            
            RateMap(i,:) = (HistClust./HistOccInterp)*f; % interpolation rate
            
        end
        
        % no sorting
        
        for i=1:17
            fig = figure(i);
            subplot(10,1,lap);
            rectangle('Position',[160,0,40,1],'FaceColor','g','EdgeColor','g')
            hold on
            rectangle('Position',[100,0,10,1],'FaceColor','c','EdgeColor','c')
            hold on
            rectangle('Position',[260,0,10,1],'FaceColor','c','EdgeColor','c')
            hold on
            bar(Edges(2:end),RateMap(i,:));
            ylabel(['lap ' num2str(lap)])
            ylim([0 1])
            if lap == 1
                counter = (maze-1)*17+i;
                title(['Maze ' num2str(maze) ': Rate Map (Interpolated) for cluster no. ' num2str(i) ' - TT' num2str(tt_no(counter)) ' - cl' num2str(cluster_no(counter))]);
            elseif lap == maxlap
                % save as PDF
                set(fig,'Visible','on');
                set(fig,'PaperUnits', 'inches','PaperPosition', [0 0 8.5 11]);
                saveas(fig,[ 'cluster' num2str(i) '- maze' num2str(maze) '.pdf'])
            end
        end
        
    end
end