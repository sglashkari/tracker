%%READSIDETRACKINGDATA reads tracking.dat file which is the output of the
%%trackrat softeware and creates a tracking.csv as the output
% update of the original version to read from any directory
% SGL 2022-01-02
% For calculating the length of the ditch

clc; clear; close all

%% Side View
[trackingFile,path] = uigetfile('tracking.dat','Select a tracking data file');
Filename = fullfile(path,'..','side','tracking.dat');
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
l_side = (0.0308 * xr - 17.114 - (0.0308 * xl - 4.4788)) * 2.54; % inch to cm
frame_no_side = [Samples.frame]';

figure(1); clf;
%plot(frame_no_side,l_side); hold on
l_side = movmedian(l_side,1000); % length of ditch
hold on
plot(frame_no_side,l_side,'r');
hold on
plot(frame_no_side,movmedian((0.0308 * xl - 4.4788) * 2.54,1000),'g');

%% Top view
Filename = fullfile(path,'..','top','tracking.dat');
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
Filename = fullfile(path,'..','side','times.csv');
T_side = readmatrix(Filename);
Filename = fullfile(path,'..','top','times.csv');
T_top = readmatrix(Filename);
t_side = t(1) + round(T_side(frame_no_side+1,5) - T_top(1,5),3); % seconds (with milliseconds precision)

% interpolation
ditch_length = round(interp1(t_side,l_side,t,'linear','extrap'),2);

%% Pose
Filename = fullfile(path,'..','pose.csv');
if isfile(Filename)
    P = readtable(Filename);
    if height(P)~=length(t)
        disp([Filename ' does not have the right number of elements!'])
        return;
    end
    head(P,10)
    pos = [P.trans_mark_wrt_base_x P.y_2 P.z_2];
    rot = [P.rot_mark_wrt_base_x P.y_3 P.z_3 P.w_1];
    success = P.success_code;
    % figure
    figure(5)
    ax1 = subplot(2,1,1);
    plot(frame_no,x,'.b')
    ax2 = subplot(2,1,2);
    plot(frame_no,P.trans_mark_wrt_base_x,'.r')
    linkaxes([ax1 ax2],'x');
else
    pos = zeros(length(t),3);
    rot = [zeros(length(t),3) ones(length(t),1)];
    disp([Filename ' does not exist!'])
end

%%
T = table(frame_no,t,flag,x,y,ditch_length, success, pos, rot);
fprintf('Accuracy of 2D tracking is %.2f%%\n',100*sum(x>-1)/length(t))
if max(flag == 1)
    fprintf('Accuracy of high confidence 2D tracking is %.2f%%\n',100*sum(flag>0)/length(t))
end
Filename = fullfile(path,'..','full-tracking.csv');
writetable(T,Filename,'Delimiter',',')

%% 
T = readtable(Filename);
head(T,10)