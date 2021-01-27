Nlx_directory = '/home/shahin/Downloads/Day3/';
[defaults,Nlx_directory] = uigetfile('*.*', 'Select Defaults File',fullfile(Nlx_directory,'Defaults'));
if isequal(defaults, 0)
    error('Defaults file was not selected!')
end
Filename = fullfile(Nlx_directory, defaults);

fileID=fopen(Filename);
datacell = textscan(fileID,'%s', 'HeaderLines', 6,'Delimiter',',');
fclose(fileID);
A = string(datacell{1});
idx = find(contains(A,'"Epoch"'));
A(idx+1)