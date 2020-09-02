clc; clear; close all;
addpath('..'); %adding readbinary to the path
file_name = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx\TT14\TT14.ntt.parms';
[cluster, timetable, table, header] = readnttparam(file_name);
size(data)

%% these are the format of the file for ntt.param
% format(:,1) = mat2cell(['double'; repmat('single',[26 1])],ones(1,27));
% format(:,2) = mat2cell(ones(27,2),ones(1,27));
% format(:,3) = split('Timestamp,MPeakX,MPeakY,MPeakA,MPeakB,PreVallEyX,PreValleyY,PreValleyA,PreValleyB,SpikeID,EnergyX,EnergyY,EnergyA,EnergyB,MaxHeight,MaxWidth,XPos,YPos,Time,PeakX,PeakY,PeakA,PeakB,ValleyX,ValleyY,ValleyA,ValleyB',",");
