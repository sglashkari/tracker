function [Time,Data,Header] = readtt(Filename, ExtractionModeVector)
if nargin == 0
     Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
     Filename = fullfile(Nlx_directory,'TT2',"TT2.ntt");
 end
FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, Spike Channel Numbers, Cell Numbers, Spike Features, Samples
HeaderExtractionFlag = 1;

if nargin < 2
    ExtractionMode = 1;
    ExtractionModeVector = [];
else
    ExtractionMode = 4;
end

if ispc
    addpath('../../pkgs/MatlabImportExport_v6.0.0'); % Neuralynx packages for Windows
    [Timestamps, ScNumbers, CellNumbers, Features, Samples, ...
        Header] =  Nlx2MatSpike(Filename, FieldSelectionFlags, ...
        HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
else
    addpath('../../pkgs/Nlx2Mat_relDec15/binaries'); % Neuralynx packages for Linux/Mac packages
    [Timestamps, ScNumbers, CellNumbers, Features, Samples, ...
        Header] =  Nlx2MatSpike_v3(Filename, FieldSelectionFlags, ...
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
ADBitVolts = str2num(ADBitVoltsString);

N = length(Samples);
ADBitVoltsMatrix = repmat(ADBitVolts,32,1,N);

Data = Samples.*ADBitVoltsMatrix;

Time = Timestamps; % micrseconds
Time = (Time * 1e-6)'; %seconds

if nargout == 0
    plot((0:31)/SamplingFrequency,Data(:,:,randi(N)));
    xlim([0 31/SamplingFrequency])
end
    
end