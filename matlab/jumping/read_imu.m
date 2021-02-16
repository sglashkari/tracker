function imu = read_imu(Nlx_directory)
% Reading IMU data from Neuralynx directory
% Sept 28, 2019
% Author Shahin G. Lashkari

addpath(genpath('..'));

if nargin < 1
    Nlx_directory = '/Users/shahin/Desktop/test/2019-09-18/';
end
IMU_files_info = dir(fullfile(Nlx_directory,'BASE_*.ncs'));
IMU_files = fullfile(Nlx_directory,{IMU_files_info.name}');

if ~isequal(length(IMU_files),6)
    error('IMU files not found!\n%s','Please enter a valid Neuralynx directory.');
end
sensor = extractBetween({IMU_files_info.name},'BASE_','.ncs');FieldSelectionFlags = [1 0 0 0 1];

Mlx2MatCSC_both = {
    @(x) Nlx2MatCSC_v3(x, FieldSelectionFlags,0,1,[]); % 1 for Linux or Mac
    @(x) Nlx2MatCSC(x, FieldSelectionFlags,0,1,[]);  % 2 for Windows
    };
 
if isunix
    addpath('../../pkgs/releaseDec2015/binaries');  % Linux or Mac
    [Timestamps,Samples] = cellfun(Mlx2MatCSC_both{1}, IMU_files, 'UniformOutput', false);
else
    addpath('../../pkgs/MatlabImportExport_v6.0.0'); % Windows
    [Timestamps,Samples] = cellfun(Mlx2MatCSC_both{2}, IMU_files, 'UniformOutput', false);
end

% vectorization
samples_modified = cellfun(@(x) x(:),Samples, 'UniformOutput', false);
ns = 1:length(samples_modified{1}); % number of samples

% NLX IMU sampling rate is 3kHz so the interval is 333 microseconds 
% but its timestamping is at the rate of 5.859Hz and the
% interval of 170667 microseconds
ts = Timestamps{1};
n = length(ts) - 1; % number of timestamps minus the last one
nt = 1:512:(512*n+1) ; % timestamps happen every 512 samples

p = polyfit(nt,ts, 1);
timestamp_modified = polyval(p,ns)'; % modified NLX timestamp

sensor = cellfun(@lower,sensor,'UniformOutput',false); %lowercase

imu = cell2struct(samples_modified, sensor, 1);

imu.a = [imu.lax imu.lay imu.laz];
imu.w = [imu.avx imu.avy imu.avz];
% nan for saturation
% imu.avx(imu.avx==min(imu.avx) | imu.avx==max(imu.avx))=nan;
% imu.avy(imu.avy==min(imu.avy) | imu.avy==max(imu.avy))=nan;
% imu.avz(imu.avz==min(imu.avz) | imu.avz==max(imu.avz))=nan;
% imu.lax(imu.lax==min(imu.lax) | imu.lax==max(imu.lax))=nan;
% imu.lay(imu.lay==min(imu.lay) | imu.lay==max(imu.lay))=nan;
% imu.laz(imu.laz==min(imu.laz) | imu.laz==max(imu.laz))=nan;

imu.time = timestamp_modified * 1e-6; %time in seconds
imu.t = imu.time - imu.time(1); % time from start in seconds
end