%pos_p_gen_jumping generates Pos.p for jumping rig experiments
% October 27, 2020
% heavily adapted from pos_p_gen2 and pos_p_gen
%
% Author Shahin G Lashkari
clear; 
clc; close all

exp_directory = 'D:\OneDrive - Johns Hopkins\JHU\913_Jumping_Recording\2020-10-25_Rat913-01\'; %day 1 

% %a new folder will be created % tmp
Nlx_directory = fullfile(exp_directory,'Neuralynx'); %tmp
% Nlx_directory = 'C:\Users\dome3neuralynx\Desktop\Day4 - Copy' %tmp
% exp_directory = Nlx_directory; %tmp
pixels_width = 640;
pixels_height = 480;
frame_rate = 30; % fps
%% Read data from file

tracking_file_name = ls(fullfile(exp_directory,'Videos','*.mat')); %tmp
load(fullfile(exp_directory,'Videos',tracking_file_name),'position');

flag = (position(:,2) > 0);

x = position(flag,2);
y = position(flag,3);
t = position(flag,1)/frame_rate; % in seconds
angle = ones(size(x)) * -99;

x = x - min(x);
y = y - min(y);

camera_image_width = max(x);
camera_image_height = max(y);

% calculating offset by software sychronization between ROS and Neuralynx
event_file_name = fullfile(Nlx_directory,'Events.nev');
% [T5, T6] = soft_sync_cam_nlx(event_file_name);

T5 = 29.056608915328979;
T6 = 29.660408973693848;

% Day1 530 542
% Day2 9038 9050
% Day3 3976 3988
% Day4 1165 1177 = 231 243
% 913: 
% Day1 626 643

light_on_cam = 626/frame_rate; % in sec
light_off_cam = 643/frame_rate; % in sec

offset_1 = T5 - light_on_cam * 1e6;
offset_2 = T6 - light_off_cam * 1e6;

T6-T5
light_off_cam-light_on_cam
offset = (offset_1 + offset_2)/2;


%% Calculations
ts = t * 1e6 + offset; % in microseconds

ratio_width = pixels_width / camera_image_width;
ratio_height = pixels_height / camera_image_height;
[ratio, argmin] = min([ratio_width,ratio_height]);

pos_x = x * ratio + (argmin == 2) * (pixels_width - camera_image_width * ratio) / 2 ;
pos_y = y * ratio + (argmin == 1) * (pixels_height - camera_image_height * ratio) / 2 ;

%crop
pos_x = max(0,pos_x);
pos_y = max(0,pos_y);
pos_x = min(pixels_width,pos_x);
pos_y = min(pixels_height,pos_y);

figure; plot(pos_x, pixels_height - pos_y,'.')
axis equal
axis([0 pixels_width 0 pixels_height])

figure; plot(ts,pos_x,'ob')
figure; plot(ts,pos_x,'ob')

% interpolation
tsq = (ts(1):1e4:ts(end))';
pos_xq = interp1(ts,pos_x,tsq);
pos_yq = interp1(ts,pos_y,tsq);
ts = tsq;
pos_x = pos_xq;
pos_y = pos_yq;
angle = ones(size(pos_x)) * -99;

hold on; plot(ts,pos_x,'.r')

%% Writing header of the files
event_file_info = dir(event_file_name);
start_time = datetime(event_file_info.date)
header = {"%%BEGINHEADER";...
       "%Stadium-style rig for jumping rats: Created from 2D head tracking data";...
       strcat("%Date of the Experiment: ",datestr(start_time,"mmmm dd, yyyy"));...
       strcat("%Time of the Experiment: ",datestr(start_time,"HH:MM AM"));...
       "%File Type: Binary";...
       "%Format:";...
       "%--------------------------";...
       "%|  Timestamp   (Double)  |";...
       "%|------------------------|";...
       "%|  X coord     (Single)  |";...
       "%|------------------------|";...
       "%|  Y coord     (Single)  |";...
       "%|------------------------|";...
       "%|  Direction   (Single)  |";...
       "%--------------------------";...
       "%";...
       "%Use following record type to read in:";...
       "%Public Type PosRecord";...
       "%    TS As Double";...
       "%    X As Single";...
       "%    Y As Single";...
       "%    Dir As Single";...
       "%End Type";...
       "%%ENDHEADER"};

pos_p_file_name = fullfile(Nlx_directory,'Pos.p');
pos_p_ascii_file_name = fullfile(Nlx_directory,'Pos.p.ascii');

if isfile(pos_p_file_name)
    answer=questdlg('The Pos.p file already exists! Do you want to overwrite it?', 'Warning!', 'Yes', 'No', 'Yes');
    if ~isequal(answer,'Yes')
        disp('Pos.p was not overwritten!')
        return
    end
end

file = fopen(pos_p_file_name,'w');
file_ascii = fopen(pos_p_ascii_file_name,'w');
format_spec = '%s\r\n';

for row = 1:size(header,1)
    fprintf(file,format_spec,header{row,:});
    fprintf(file_ascii,format_spec,header{row,:});
end

%% Writing data to the file in binary and ascii format
positions = [ts, pos_x, pos_y, angle];
format_spec_ascii = '%.2f, %.4f, %.4f, %.4f\r\n';

for row = 1:size(positions,1)
    fwrite(file,ts(row),'float64');
    fwrite(file,positions(row,2:4),'float32');
    fprintf(file_ascii,format_spec_ascii,positions(row,:));
end

fclose(file);
disp('Pos.p created!')
fclose(file_ascii);
disp('Pos.p.ascii created!')