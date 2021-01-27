function [Time,X,Y, Dir, Header] = readposp(Filename, StartTime, EndTime)
%%READBINARY reads pos.p file with the header
% StartTime is the vector of initial times (one value for each maze)
% EndTime is the vector of final times (one value for each maze)
% Times are in seconds
% sgl 2020-08-30

if nargin == 0
    Nlx_directory = '/home/shahin/Downloads/Day3/';
    [posp,Nlx_directory] = uigetfile(fullfile(Nlx_directory,'Pos.p'), 'Select Pos.p File');
    if isequal(posp, 0)
        error('Pos.p file was not selected!')
    end
    Filename = fullfile(Nlx_directory, posp);
end
if nargin < 3
    [Time,Data,Header,EventIDs,TTLs] = readevent(fullfile(Nlx_directory,'Events.nev'));
    
    % for Day 3 rat 913
    idx = EventIDs==4;
    StartTime = Time(idx);
    EndTime = Time(find(idx)+1);
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