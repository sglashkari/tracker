function combinentt(Filenames)
%%% combine multiple ntt files into a single ntt file.
% Example:
%       Filenames = {'/home/shahin/test/matlab/TT6_0000.ntt', '/home/shahin/test/matlab/TT6_0001.ntt'}
% 
% Shahin 2021-01-20
%
if nargin == 0
    [file,path] = uigetfile('*.ntt', 'Select Two or More NTT Files', 'MultiSelect', 'on');
    Filenames = fullfile(path,file);
end
%   INPUT ARGUMENTS:
FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, Spike Channel Numbers, Cell Numbers, Spike Features, Samples
HeaderExtractionFlag = 1; % header import is desired
ExtractMode = 1; % Extract All
ExtractionModeVector = []; % Extract All (The vector value is ignored)

%   OUTPUT VARIABLES:
%   Timestamps: A 1xN integer vector of timestamps.
%   ScNumbers: A 1xN integer vector of spike channel numbers. This is the order that
%              the spike AEs were created and have nothing to do with the
%              AD channel number.
%   CellNumbers: A 1xN integer vector of classified cell numbers. If no cell was
%                classified for this spike, this value will be zero.
%   Features: A 8xN integer vector of the features (e.g. Peak, Valley, etc.) calculated
%             by Cheetah.
%   Samples: A 32x4xN integer matrix of the data points. Where M is the number of
%            subchannels in the spike file (NTT M = 4).
%            These values are in AD counts.
%   Header: A Mx1 string vector of all the text from the Neuralynx file header, where
%           M is the number of lines of text in the header.

Timestamps = [];

for i = 1:length(Filenames)
    Filename = Filenames{i};
    [Timestamps1, ScNumbers1, CellNumbers1, Features1, Samples1, Header1] = ...
        Nlx2MatSpike( Filename, FieldSelectionFlags, ...
        HeaderExtractionFlag, ExtractMode, ExtractionModeVector);
    
    if i == 1
        Timestamps = Timestamps1;
        ScNumbers = ScNumbers1;
        CellNumbers = CellNumbers1;
        Features = Features1;
        Samples = Samples1;
        Header = Header1;
    else
        Timestamps = [Timestamps; Timestamps1];
        ScNumbers = [ScNumbers;ScNumbers1];
        CellNumbers = [CellNumbers;CellNumbers1];
        Features = [Features;Features1];
        Samples = cat(3,Samples,Samples1);
    end
end

Filename = strcat(extractBefore(Filename,'.ntt'), '_combined.ntt');
AppendToFileFlag = 0; % Delete file if exists
ExportMode = 1; % Export All
ExportModeVector = []; % (Extract All): The vector value is ignored.

Mat2NlxSpike( Filename, AppendToFileFlag, ExportMode, ExportModeVector, ...
    FieldSelectionFlags, Timestamps, ScNumbers, CellNumbers, ...
    Features, Samples, Header);

end