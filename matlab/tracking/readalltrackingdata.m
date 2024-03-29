%%READSIDETRACKINGDATA reads tracking.dat file which is the output of the
%%trackrat softeware and creates a tracking.csv as the output
% update of the original version to read from any directory
% SGL 2022-01-02
% For calculating the length of the ditch

clc; clear; close all

%% Side View
path = uigetdir('/media/dome3tracking/','Select an experiment directory');
idx = find(path=='/', 1, 'last');
date_str = path(idx(end)+1:end);
rat_no_str = char(extractBetween(path,'Rat_',filesep));
%%
Filename = fullfile(path,'side','tracking.dat');
if ~isfile(Filename)
    fprintf([Filename ' does not exist!\n']);
    return;
end

FileID = fopen(Filename,'r');

format(:,1) = {'uint32';'uint32';'double'; 'uint32'; 'uint32'; 'uint32'; 'uint32'; 'single'; 'single'};
format(:,2) = mat2cell(ones(9,2),ones(1,9));
format(:,3) = split('frame,flag,time,p1,p2,p3,p4,x,y',",");

m = memmapfile(Filename,'Format',format);
fclose(FileID);

Samples = m.Data;
xr = [Samples.x]';
xl = [Samples.y]';
%l_side = (0.0308 * xr - 17.114 - (0.0308 * xl - 4.4788)) * 2.54; % inch to cm (Rat 980)
%l_side = (0.0257 * (xr - xl) - 13.01) * 2.54; % inch to cm (Rat 1024)
%l_side = (0.0291 * (xr - xl) - 12.298) * 2.54;  % inch to cm (Rat 1055 old)
l_side = (0.0297 * (xr - xl) - 12.421) * 2.54;  % inch to cm (Rat 1055, 1068)
frame_no_side = [Samples.frame]';

figure(10); clf;
plot(frame_no_side, movmedian(xr - xl, 100));

figure(1); clf; hold on;
l_side = movmedian(l_side,500); % length of ditch
yyaxis left
plot(frame_no_side,l_side,'b');
ylabel('length (cm)')
yyaxis right
plot(frame_no_side,l_side / 2.54,'r');
ylabel('length (in)')
%plot(frame_no_side,movmedian((0.0308 * xl - 4.4788) * 2.54,2000),'g'); % (Rat 980)
%plot(frame_no_side,movmedian(0.0257 * xl * 2.54,2000),'g'); % (Rat 1024)
plot(frame_no_side,movmedian(0.0297 * xl * 2.54,500),'g'); % (Rat 1055, 1068)
%% Top view
Filename = fullfile(path,'top','tracking.dat');
if ~isfile(Filename)
    fprintf([Filename ' does not exist!\n']);
    return;
end
FileID = fopen(Filename,'r');
m = memmapfile(Filename,'Format',format);
fclose(FileID);

Samples = m.Data;
x = [Samples.x]';
y = [Samples.y]';
p1 = [Samples.p1]';
p2 = [Samples.p2]';
p3 = [Samples.p3]';
p4 = [Samples.p4]';
p = [p1 p2 p3 p4];
flag = [Samples.flag]';
frame_no = [Samples.frame]';
time = [Samples.time]'; % seconds (wrapped)
t = unwrap((time-64)/64*pi)/pi*64; % range 0 .. 128

%
figure(2)
plot(t,x,'.')
hold on
plot(t(flag==0 & x>-1),x(flag==0 & x>-1),'or');
figure(3)
plot(x,y,'.')
axis equal
axis([0 2048 0 400])
figure(4)
plot(frame_no,x,'.')

%% Time diff between the cameras
Filename = fullfile(path,'side','times.csv');
copyfile(Filename, fullfile(path, 'side_times.csv'));
T_side = readmatrix(Filename);
Filename = fullfile(path,'top','times.csv');
copyfile(Filename, fullfile(path, 'top_times.csv'));
copyfile(fullfile(path,'side','video.avi'), fullfile(path, ['Rat' rat_no_str '_' date_str '.avi']));
T_top = readmatrix(Filename);
t_side = t(1) + round(T_side(frame_no_side+1,5) - T_top(1,5),3); % seconds (with milliseconds precision)

% interpolation
ditch_length = round(interp1(t_side,l_side,t,'linear','extrap'),2);

%% Pose
Filename = fullfile(path,'pose.csv');
if isfile(Filename)
    P = readtable(Filename);
    if height(P)~=length(t)
        disp([Filename ' does not have the right number of elements!'])
        pos = zeros(length(t),3);
        rot = [zeros(length(t),3) ones(length(t),1)];
        success = zeros(length(t),1);
    else
        head(P,10)
        try % new
            pos = [P.trans_mark_wrt_ref_x P.y_2 P.z_2];
            rot = [P.rot_mark_wrt_ref_x P.y_3 P.z_3 P.w_1];
            clust_pos = [P.clust_pos_on_wrld_plane_x P.y_4];
        catch % old (before 2022-03-05 update)
            pos = [P.trans_mark_wrt_base_x P.y_2 P.z_2];
            rot = [P.rot_mark_wrt_base_x P.y_3 P.z_3 P.w_1];
        end
        success = P.success_code;
        fprintf('Accuracy of 3D tracking is %.2f%%\n',100*sum(success>0)/length(t))
        % figure
        figure(5); clf
        ax1 = subplot(2,1,1);
        plot(frame_no,x,'.b')
        ax2 = subplot(2,1,2);
        plot(frame_no(success>0),clust_pos(success>0,1),'.r'); hold on
        plot(frame_no(success>0),pos(success>0,1),'.b');
        legend('x: center of cluster','x: crown ref frame')
        xlabel('Frame No')
        ylabel('Horizontal position (m)')
        linkaxes([ax1 ax2],'x');
        figure(6); clf
        plot3(pos(:,1), pos(:,2), pos(:,3), '.')
        axis equal
    end
else
    pos = zeros(length(t),3);
    rot = [zeros(length(t),3) ones(length(t),1)];
    success = zeros(length(t),1);
    disp([Filename ' does not exist!'])
end

%%
T = table(frame_no,t,flag,x,y,ditch_length, success, pos, rot);
fprintf('Accuracy of 2D tracking is %.2f%%\n',100*sum(x>-1)/length(t))
if max(flag == 1)
    fprintf('Accuracy of high confidence 2D tracking is %.2f%%\n',100*sum(flag>0)/length(t))
end
Filename = fullfile(path,'full-tracking.csv');
writetable(T,Filename,'Delimiter',',')

%%
clear T
T = readtable(Filename);
head(T,10)