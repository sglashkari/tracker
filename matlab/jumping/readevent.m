function [Time,Data,Header,EventIDs,TTLs] = readevent(Filename)
if nargin == 0
    Nlx_directory = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Neuralynx';
    Filename = fullfile(Nlx_directory,'Events.nev');
end
FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, Event IDs, TTLs, Extras, Event Strings
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];


if ispc
    addpath('../../pkgs/MatlabImportExport_v6.0.0'); % Neuralynx packages for Windows
    [Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = ...
        Nlx2MatEV( Filename, FieldSelectionFlags, HeaderExtractionFlag, ...
        ExtractionMode, ExtractionModeVector );
else
    addpath('../../pkgs/releaseDec2015/binaries'); % Neuralynx packages for Linux/Mac packages
    [Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = ...
        Nlx2MatEV_v3( Filename, FieldSelectionFlags, HeaderExtractionFlag, ...
        ExtractionMode, ExtractionModeVector );
end

% InputFormat = 'yyyy/MM/dd HH:mm:ss';
% StartTimeString = extractAfter(Header{7},'-TimeCreated ');
% StartTime = datetime(StartTimeString,'InputFormat',InputFormat);
% FinishTimeString = extractAfter(Header{8},'-TimeClosed ');
% FinishTime = datetime(FinishTimeString,'InputFormat',InputFormat);

Data = EventStrings;%(EventIDs==4);
Time = Timestamps;%(EventIDs==4); % micrseconds

% % only for one day, modify events file for future (SGL 9 Aug, 2020)
% if ~isempty(strfind(Filename, '200329'))
%     Data(5:end+1)=Data(4:end);
%     Data(1:4)={'start maze 1', 'stop maze 1','start maze 2','stop maze 2'};
%     Time(5:end+1)=Time(4:end);
%     Time(4)= 5100000000;
% end

Time = (Time * 1e-6)'; %seconds

if nargout == 0
    Data
end

end