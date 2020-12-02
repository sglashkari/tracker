%
Nlx_directory = '~/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02/Neuralynx';
Nlx_directory = uigetdir(Nlx_directory,'Select Neuralynx Folder');

listing = dir(fullfile(Nlx_directory,'TT*','*.ntt'));
names = string({listing.name}');
folders = string({listing.folder}');
N = length(listing); % number of maze-clusters

absolue_paths = mat2cell(fullfile(folders,names),ones(N,1));

% Compelete undo
if isfolder(fullfile(Nlx_directory,'old_files'))
    disp('old_files directory already exists.')
    return;
else
    mkdir(fullfile(Nlx_directory,'old_files'));
end
for i=1:N
    movefile(fullfile(folders{i},'*.ntt'),Nlx_directory)
    if isfile(fullfile(folders{i},'*.ntt.parms'))
        delete(fullfile(folders{i},'*.ntt.parms')) 
    end
    if isfile(fullfile(folders{i},'Defaults'))
        delete(fullfile(folders{i},'Defaults'))
    end
    disp(['folder ' folders{i} ' is undone!'])
end
movefile(fullfile(Nlx_directory,'TT*'),fullfile(Nlx_directory,'old_files'))
delete(fullfile(Nlx_directory,'Defaults'))
 
%copyfile(fullfile(Nlx_directory,'old_files','TT*'),Nlx_directory)
