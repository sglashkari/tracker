%% Reads a video frame by frame
% Updated Sept 16, 2020
% Showing the frame given the neuralynx time
% Author Shahin G Lashkari
clc; clear; close all
vid_filename = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Videos\2020-03-29_20-22-00.mp4';
videoReader  = VideoReader(vid_filename)

framerate = 30;

time = 5966.75;

offset = 3887091990.1667* 1e-6 ; %seconds

frame_no = round((time - offset) * framerate);

datestr(seconds(frame_no/30),'HH:MM:SS')

while true
    time = round(frame_no/framerate+offset,3);
    dim = [0.12 0.83 0.1 0.1];
    frame_image = read(videoReader,frame_no);
    imshow(frame_image);
    annotation('textbox',dim,'String',['frame = ', num2str(frame_no)],'FitBoxToText','on','Color','white');
    annotation('textbox',dim+[0 0 0 0.03],'String',['time = ' num2str(time)],'FitBoxToText','on','Color','white');
    
    try
        w = waitforbuttonpress;
    catch
        break;
    end
    
    click_type = get(gcf, 'SelectionType');
 
    if strcmp(click_type,'normal') %right click
        frame_no = frame_no + 1;
    elseif strcmp(click_type,'alt') %left click
        frame_no = frame_no - 1;
    end

end

tic
S = zeros(100,1);
Z = zeros(1080, 1920, 3, 100, 'uint8');
for frame_counter=1:length(S)
    R = read(videoReader,frame_counter+frame_no);
    S(frame_counter)=sum(sum(R(:,:,1)));
    Z(:,:,:,frame_counter) = R;
end

toc;

imgmode = mode(Z);
toc;

imshow(imgmode)

datacursormode on

%% Top Cam | Side Cam
% Day1 530 542
% Day2 9038 9050
% Day3 3976 3988
% Day4 1165 1177    offset  3887091990.1667 |  231 243 offset 3918225323.5