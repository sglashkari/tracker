function [cluster, timetable, table, header] = readnttparam(file_name)
%%READBINARY reads binary file with the header
%
%% old method
% fileID = fopen(file_name,'r');
% header = arrayfun(@(~) fgetl(fileID), 1:header_length, 'UniformOutput', false);
% row = 0;
% while true
%     row = row + 1;
%     try
%         data(row,1) = fread(fileID,1,'double');
%         data(row,2:data_size) = fread(fileID,data_size-1,'single');
%     catch
%         break;
%     end
% end
% fclose(fileID);
%% new method
fileID = fopen(file_name,'r');
header = arrayfun(@(~) fgetl(fileID), 1:48, 'UniformOutput', false);
header_length = ftell(fileID);
format(:,1) = mat2cell(['double'; repmat('single',[26 1])],ones(1,27));
format(:,2) = mat2cell(ones(27,2),ones(1,27));
format(:,3) = split('Timestamp,MPeakX,MPeakY,MPeakA,MPeakB,PreVallEyX,PreValleyY,PreValleyA,PreValleyB,SpikeID,EnergyX,EnergyY,EnergyA,EnergyB,MaxHeight,MaxWidth,XPos,YPos,Time,PeakX,PeakY,PeakA,PeakB,ValleyX,ValleyY,ValleyA,ValleyB',",");

m = memmapfile(file_name,'Format',format, 'Offset', header_length);
cluster = m.Data;
table = struct2table(cluster);
table.Time = seconds(table.Time);
timetable = table2timetable(table);
head(timetable)
fclose(fileID);
end