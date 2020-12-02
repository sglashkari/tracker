%% Extract data from the output of winclust (cluster mazes) and pos.p (ascii) 
% and transform it to MAT format

%% Neural data
exp_directory = '~/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02';
exp_directory = uigetdir(exp_directory,'Select Experiment Folder');
if exp_directory == 0
    return;
end
Nlx_directory = fullfile(exp_directory,'Neuralynx');

listing = dir(fullfile(Nlx_directory,'**','cl-maze*.*'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing); % number of maze-clusters

absolue_paths = mat2cell(fullfile(folders,names),ones(N,1));

tic
A = cellfun(@(x) importdata(x,',',13), absolue_paths, 'UniformOutput', false);

%% Postion data
pos_p_filename = fullfile(Nlx_directory,'Pos.p.ascii');
B = importdata(pos_p_filename,',',24);

toc
mat_filename = fullfile(exp_directory,'data.mat');
save(mat_filename,'A','B')
disp(['File ' mat_filename ' has been created!'])