function [t, x, y, p, frame, flag] = readtrackingdata(exp_directory)
%%READTRACKINGDATA reads tracking.dat file which is the output of the
% trackrat (git) and extracts the data
%
% if only 3 outputs are requested (t, x, y) then the result only shows the 
% correctly tracked data points (flag == 1)
% otherwise all the extracted data will spit out
% if no output is requested, then the 
% SGL 2020-12-01

if nargin == 0
    exp_directory = '/home/shahin/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02';
    exp_directory = uigetdir(exp_directory,'Select Experiment Folder');
    if exp_directory == 0
        return;
    end
end

Filename = fullfile(exp_directory,'Videos','tracking.dat');
FileID = fopen(Filename,'r');

format(:,1) = {'uint32';'uint32';'double'; 'uint32'; 'uint32'; 'uint32'; 'uint32'; 'single'; 'single'};
format(:,2) = mat2cell(ones(9,2),ones(1,9));
format(:,3) = split('frame,flag,time,p1,p2,p3,p4,x,y',",");

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
frame = [Samples.frame]';
time = [Samples.time]'; % seconds (wrapped)
t = unwrap((time-64)/64*pi)/pi*64; % range 0 .. 128

% If flag is not requested discard the data that are not been correctly tracked
if nargout < 6 
    x = x (flag == 1);
    y = y (flag == 1);
    t = t (flag == 1);
    p = p (flag == 1, :);
    frame = frame (flag == 1);
end

% If no output only plot the results:
if nargout == 0 
    close all
    figure(1)
    plot(x,y,'.');
    axis equal
    figure(2)
    ax1 = subplot(2,1,1);
    plot(t, x,'.');
    ax2 = subplot(2,1,2);
    plot(t, y,'.');
    linkaxes([ax1 ax2],'x');
    clear t
else
        

end