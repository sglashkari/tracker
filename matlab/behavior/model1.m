%% This model assumes that rat
%   1. uses psychometric curve 1 (reverse unit step function with
%   transition at a) after jumping
%   2. uses psychometric curve 2 (reverse unit step function with
%   transition at b) after ditching
%
%   Shahin 2021-07-18
%
clc;
clear;
close all;
load('data');

trial = 3;

data = exp{trial};
%data = cat(1, exp{:}); % ALL TRIALS
similarity_count = zeros(25);

a = 6:0.5:18;
b = 6:0.5:18;

for i = 1:length(a)
    for j = 1:length(b)
        counter = 0;
        for n=2:length(data)
            if data(n-1,3)  % jump
                if data(n,1) <= b(j)   % ditch length
                    prediction = 1;
                else
                    prediction = 0;
                end
            else % ditch
                if data(n,1) <= a(i)   % ditch length
                    prediction = 1;
                else
                    prediction = 0;
                end                
            end
            counter = counter + (prediction == data(n,3));
        end
        similarity_count(i,j) = counter;
    end
end

% imagesc(similarity_count);
% caxis([0 length(data)]);
% colorbar;


%[px,py] = gradient(similarity_count);


% similarity_count==max(max(similarity_count));

%figure
%contour(a,b,similarity_count./length(a)./length(b))
% hold on
% quiver(a,b,px,py)
% hold off


figure
h = heatmap(a,b,100*round(similarity_count./length(data),2));

h.Title = ['Percentage of successful predictions for a pair of thresholds at trial ' num2str(trial)];
h.XLabel = 'Ditch Length Thresh 1';
h.YLabel = 'Ditch Length Thresh 2';
h.Colormap = jet;

set(gcf, 'Position', [100 100 1000 1000]);