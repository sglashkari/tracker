function [t, x, y, p, frame, flag] = readtrackingdata(exp_directory)
%%READBINARY reads tracking.dat file with the header
% StartTime is the vector of initial times (one value for each maze)
% EndTime is the vector of final times (one value for each maze)
% Times are in seconds
% SGL 2020-12-01

if nargin == 0
    exp_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
    Filename = fullfile(exp_directory,'Videos','tracking.dat');
end
if nargin < 3
    TimeEV = readevent(fullfile(Nlx_directory,'Events.nev'));
    StartTime = TimeEV(1:2:end);
    EndTime = TimeEV(2:2:end);
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
time = [Samples.time]'; %seconds
t = unwrap((time-64)/64*pi)/pi*64; % range 0 .. 128


plot((time-64)/64*pi)

hold on
t = unwrap((time-64)/64*pi)/pi*64;
plot(t)
figure(2)
subplot(2,1,1)
plot(frame,x,'.')
subplot(2,1,2)
plot(frame,y,'.')
% % Time Range
% idx = [];
% for i=1:length(StartTime)
%     idx = [idx; find(Time >= StartTime(i) & Time <= EndTime(i))];
% end
% 

if nargout == 0 % If no output only plot the results:
    close all
    figure(1)
    plot(X, 480 - Y,'.');
    axis equal
    axis([0 640 0 480])
    figure(2)
    plot(Time, Y,'.',Time(idx), Y(idx),'*');
    clear t
elseif narout == 3 % discard the ones that are not been correctly tracked
    x = x (flag == 1);
    y = y (flag == 1);
    t = t (flag == 1);
end
        

end