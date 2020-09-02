function [Time,X,Y, Dir, Header] = readposp(Filename, StartTime, EndTime)
%%READBINARY reads pos.p file with the header
% StartTime is the vector of initial times (one value for each maze)
% EndTime is the vector of final times (one value for each maze)
% Times are in seconds
% sgl 2020-08-30

if nargin == 0
    Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
    Filename = fullfile(Nlx_directory,'pos.p');
end
if nargin < 3
    TimeEV = readevent(fullfile(Nlx_directory,'Events.nev'));
    StartTime = TimeEV(1:2:end);
    EndTime = TimeEV(2:2:end);
end
FileID = fopen(Filename,'r');
Header = arrayfun(@(~) fgetl(FileID), 1:24, 'UniformOutput', false)';
header_length = ftell(FileID);
format(:,1) = mat2cell(['double'; repmat('single',[3 1])],ones(1,4));
format(:,2) = mat2cell(ones(4,2),ones(1,4));
format(:,3) = split('Timestamp,XPos,YPos,Dir',",");

m = memmapfile(Filename,'Format',format, 'Offset', header_length);
fclose(FileID);

Samples = m.Data;
X = [Samples.XPos]';
Y = [Samples.YPos]';
Dir = [Samples.Dir]';
Time = [Samples.Timestamp]'* 1e-6; %seconds

% Time Range
idx = [];
for i=1:length(StartTime)
    idx = [idx; find(Time >= StartTime(i) & Time <= EndTime(i))];
end

% If no output only plot the results:
if nargout == 0
    close all
    figure(1)
    plot(X, 480 - Y,'.');
    axis equal
    axis([0 640 0 480])
    figure(2)
    plot(Time, Y,'.',Time(idx), Y(idx),'*');
    clear Time
else
    X = X(idx);
    Y = Y(idx);
    Dir = Dir(idx);
    Time = Time(idx);
end

end