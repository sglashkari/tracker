function [Time,Data,Header] = readcsc(Filename, TimeRange)
if nargin == 0
     exp_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04';
     exp_directory = '/home/shahin/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02';
     exp_directory = '/home/shahin/Desktop/20-12-09';
     Filename = fullfile(exp_directory, 'Neuralynx', 'CSC4.ncs')
end
FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples
HeaderExtractionFlag = 1;

if nargin < 2
    ExtractionMode = 1;
    ExtractionModeVector = [];
else
    ExtractionMode = 4;
    ExtractionModeVector = TimeRange;
end

if ispc
    addpath('../../pkgs/MatlabImportExport_v6.0.0'); % Neuralynx packages for Windows
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
        Samples, Header] = Nlx2MatCSC(Filename, FieldSelectionFlags,...
        HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
else
    addpath('../../pkgs/releaseDec2015/binaries'); % Neuralynx packages for Linux/Mac packages
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
        Samples, Header] = Nlx2MatCSC_v3(Filename, FieldSelectionFlags,...
        HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
end

InputFormat = 'yyyy/MM/dd HH:mm:ss';
StartTimeString = extractAfter(Header{8},'-TimeCreated ');
StartTime = datetime(StartTimeString,'InputFormat',InputFormat);
FinishTimeString = extractAfter(Header{9},'-TimeClosed ');
FinishTime = datetime(FinishTimeString,'InputFormat',InputFormat);

SamplingFrequencyString = extractAfter(Header{15},'-SamplingFrequency ');
SamplingFrequency = str2double(SamplingFrequencyString);

ADBitVoltsString = extractAfter(Header{17},'-ADBitVolts ');
ADBitVolts = str2double(ADBitVoltsString);

Data = Samples(:)*ADBitVolts;
N = size(Samples);
s = 1:N(1):N(1)*N(2);
sq = 1:1:N(1)*N(2);

Time = interp1(s,Timestamps,sq); % check accuracy of this method

Data(isnan(Time))=[];
Time(isnan(Time))=[];

Time = Timestamps(1) + (1:length(Data)) * 1e6 / SamplingFrequency; % micrseconds
Time = (Time * 1e-6)'; %seconds

if nargout == 0
    plot(Time,Data);
    clear Time;
end
    
end