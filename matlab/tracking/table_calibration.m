%%TABLE_CALIBRATION reads pose.csv file which is the output of Balazs's
%%tracker and calculates the 
%%tracking
%
% SGL 2022-11-27

clc; clear; close all

[Filename,path] = uigetfile('/media/dome3tracking/pose.csv','Select a tracking data file');
Filename = fullfile(path,Filename);
if ~isfile(Filename)
    fprintf([Filename ' does not exist!\n']);
    return;
end

%% Pose
clc
P = readtable(Filename);
head(P,10)
x = P.trans_mark_wrt_ref_x;
y = P.y_2;
z = P.z_2;
rot_x = P.rot_mark_wrt_ref_x;
rot_y = P.y_3;
rot_z = P.z_3;
rot_w = P.w_1;
success = P.success_code;
frame_no = P.frame_num;

% X-Z tilt
figure(1); clf; hold on
plot(x,z,'.b')

xz_coefficients = polyfit(x, z, 1);
xFit = linspace(min(x), max(x), 1000);
zFit = polyval(xz_coefficients , xFit);
plot(xFit, zFit,'.r')

% X-Y tilt
figure(2); clf; hold on
plot(x,y,'.b')

xy_coefficients = polyfit(x, y, 1);
yFit = polyval(xy_coefficients , xFit);
plot(xFit, yFit,'.r')

% quaternion (tilt correction)
q_cam_wrt_ref = [0.999817556, 0.017256834, 0.007451162, 0.003457459]; %[P.rot_cam_wrt_ref_x P.y_1 P.z_1 P.w];
q_cam_wrt_ref = [0.999963667, 0.002368844, 0.007501642, 0.003346137];
q_cam_wrt_ref = [0.999965995, 0.002408396, 0.007170988, 0.003347238];
q_cam_wrt_ref = [q_cam_wrt_ref(end) q_cam_wrt_ref(1:3)]; % MATLAB form
yaw = atan(xy_coefficients(1)); 
pitch = atan(xz_coefficients(1));
roll = 0;
q_correction = angle2quat(yaw, pitch, roll, 'ZYX');
%q_correction = angle2quat(roll, pitch, yaw, 'XYZ');
q_cam_wrt_ref_corrected = quatmultiply(q_cam_wrt_ref, q_correction);
q_cam_wrt_ref_corrected = [q_cam_wrt_ref_corrected(2:4) q_cam_wrt_ref_corrected(1)]; % Balazs form
format longG
fprintf("q_cam_wrt_ref_corrected:\n[%.9f, %.9f, %.9f, %.9f]\n", q_cam_wrt_ref_corrected(end,:))
fprintf("Error: pitch %.3f deg, yaw %.3f deg.\n", rad2deg(pitch), rad2deg(yaw))

% position
p_cam_wrt_ref = [ -0.215771515260324 , 0.289453220825853, 1.45478007719928 ]; %[P.trans_cam_wrt_ref_x P.y P.z];
format short
p_correction = mean([x y z])
figure(1)
plot(p_correction(1), p_correction(3), '+m', 'MarkerSize',48)
ylabel('Z (m)')
xlabel('X (m)')
figure(2)
plot(p_correction(1), p_correction(2), '+m', 'MarkerSize',48)
ylabel('Y (m)')
xlabel('X (m)')
p_cam_wrt_ref_corrected = p_cam_wrt_ref - p_correction;
format longG
fprintf("[%.9f, %.9f, %.9f]\n", p_cam_wrt_ref_corrected(end,:))
fprintf("Error is %.1f mm.\n", norm(p_correction)*1000)

