clc; clear; close all
exp_directory = 'D:\Analysis\2021-12-21';
mat_filename = fullfile(exp_directory,'analyzed_data.mat');
load(mat_filename,'cluster');

for i=1:length(cluster)
    for j=i+1:length(cluster)
        t1 = cluster(i).t;
        t2 = cluster(j).t;
        ist1t2 = ismember(t1,t2);
        ist2t1 = ismember(t2,t1);
        match_ratio = 100*max(mean(ist1t2),mean(ist2t1));
        if match_ratio > 10
            fprintf('cluster #%d (cluster %d from shank %d) matches cluster #%d (cluster %d from shank %d) by %.2f%%.\n',i,cluster(i).cl,cluster(i).sh,j,cluster(j).cl,cluster(j).sh,100*max(mean(ist1t2),mean(ist2t1)));
        end
    end
end