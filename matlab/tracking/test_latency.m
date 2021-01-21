%%% Testing open arena latency
% between software PTP and 
% Shahin 2021-01-20

clc; clear;
event_file = 'home/shahin/Desktop/test_latency/Events.nev';

FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, Event IDs, TTLs, Extras, Event Strings
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];
addpath('../../pkgs/releaseDec2015/binaries'); % Neuralynx packages for Linux/Mac packages

[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV_v3( event_file,...
        FieldSelectionFlags, HeaderExtractionFlag, ExtractionMode, ExtractionModeVector );
    
% addpath('../jumping');
% [Timestamps,~,~,EventIDs,TTLs] = readevent(event_file);

T1 = Timestamps(TTLs==1 & EventIDs==11);
T2 = Timestamps(TTLs==2 & EventIDs==11);
T3 = Timestamps(TTLs==3 & EventIDs==11);
T4 = Timestamps(TTLs==4 & EventIDs==11);
T5 = Timestamps(TTLs==5 & EventIDs==11);
T6 = Timestamps(TTLs==6 & EventIDs==11);