%% Bharath's code
% Load cl-files
% addpackagepath('Enhanced_rdir');
% clFiles = rdir(fullfile(expFolder,'**','*.*'),['regexp(name,''',filesep,'cl-maze\d+\w?\.\d+'')'],true);
% clFiles = {clFiles.name};
% nClust = length(clFiles);

%[expFolders,epochs] = processArgs(varargin{:});
offsetRange = [15, 25];
%offsetRange = [12.85, 13.85];
currentOffset = 20.7345;
offsets = linspace(offsetRange(1), offsetRange(2), 100);
skgInfVals = []; %array to store different values of Skaggs information score for different offsets
offsetDifferences = []; %to store different values of offset differences

C = importdata('cl-maze1.1');
clustXPos = C.data(18:18:end);
clustYPos = C.data(19:18:end);
clustTs = C.data(20:18:end);
% cX = (clustXPos*0.14)/640;
% cY = (clustYPos*0.14)/480;

P = importdata('../Pos.p.ascii');
occTs = P.data(:, 1);
occXPos = P.data(:, 2);
occYPos = P.data(:, 3);
% oX = (occXPos*0.14)/640;
% oY = (occYPos*0.14)/480;

for i = 1:length(offsets)
disp(i)
offsetDiff = offsets(i) - currentOffset;
shiftedTs = clustTs + offsetDiff*1e6; %these are the shifted timestamps for the cluster (in microseconds)
clustXNew = [];
clustYNew = []; %these are new arrays created for shifted X and Y positions
inds = [];
for j=1:length(shiftedTs)
    [~, ind] = min(abs(occTs - shiftedTs(j)));
    if ind == 1 | ind == length(occTs)
        disp('true')
        continue
    else
        inds = [inds ind];
        clustXNew = vertcat(clustXNew, occXPos(ind));
        clustYNew = vertcat(clustYNew, occYPos(ind));
    end
end

frameRate = 30; %frame rate in Hz, this is an approximation to the frame rate of Balazs' algorithm 
box_dim_in_meter = 0.14; % use 1.6 for dome, 0.14 for box
pixels_width = 640;
pixels_height = 480;
%Xedges = linspace(0, 0.14, 35);
%Yedges = linspace(0, 0.14, 35);
Xedges = linspace(0, 640, 64);
Yedges = linspace(0, 480, 48);
histClust = histcounts2(clustXNew, clustYNew, Xedges, Yedges);
histOcc = histcounts2(occXPos, occYPos, Xedges, Yedges);
rateMap = (histClust./histOcc)*frameRate;
rateMap(isnan(rateMap)) = 0;
rateMap = imgaussfilt(rateMap, 1);
hold on 
imagesc(rateMap)
colorbar
caxis([0, 5])
figName = strcat('fig_', num2str(i), '.png');
filename = fullfile('C:\Users\natub\OneDrive - Johns Hopkins University\KnierimLab\Presentations\Jim weekly updates\ratemaps with offset', figName);
saveas(gcf, filename);
close all
% skgInf = skaggsInformation(histOcc, histClust, frameRate);
% skgInfVals = [skgInfVals, skgInf];
offsetDifferences = [offsetDifferences, offsetDiff];
end
% plot(offsetDifferences, skgInfVals)
% title('Skaggs information score - Cell 1, TT10, Day 20');
% xlabel('Offset difference (secs)');
% ylabel('Skaggs information score');

