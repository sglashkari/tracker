%% Reads a video frame by frame
% Updated March 28, 2020, modified Oct 27, 2020
% detecting the LED light on/off
% Author Shahin G Lashkari
close all
default_filename = 'D:\OneDrive - Johns Hopkins\JHU\913_Jumping_Recording\2020-10-25_Rat913-01\Videos\2020-10-25_15-59-09.mp4'; %day 1
quit
try
[file, path] = uigetfile({'*.avi;*.mpg;*.mpeg;*.mp4','Video Files (*.avi,*.mpg,*.mpeg,*.mp4)'; '*.*',  'All Files (*.*)'},...
    'Select a video file',default_filename);
    videoReader  = VideoReader(fullfile(path,file))
catch
    videoReader  = VideoReader(default_filename)
end
% 
% tic
% N = 1000;
% R = zeros(N,1);
% G = zeros(N,1);
% B = zeros(N,1);
% for frame_no=1:N
%     frames = read(videoReader,frame_no);
%     R(frame_no)=max(max(frames(:,:,1))); % R = 1, G = 2, B = 3
% end
% plot(1:N,R,'r',1:N,G,'g',1:N,B,'b');
% % [~,argmax] = min(R(R>1.1*min(R)));
% toc

frame_no = input('What is the frame number? ');

figure;
while true
    dim = [0.15 0.75 0.1 0.1];
    frames = read(videoReader,frame_no);
    ax = imshow(frames);
    fig = ancestor(ax, 'figure');
    annotation('textbox',dim,'String',num2str(frame_no),'FitBoxToText','on','Color','white');
    
    try
        waitforbuttonpress;
    catch
        break;
    end
    
    button = double(get(gcf,'CurrentCharacter'));
    % 28 leftarrow
    % 29 rightarrow
    % 30 uparrow
    % 31 downarrow
    try
        if (button == 28) || (button == 31)
            frame_no = frame_no - 5;
        elseif (button == 29) || (button == 30)
            frame_no = frame_no + 1;
        end
    catch
        frame_no = frame_no + 1;
    end
end

% Rat 883
% Day1 530 542
% Day2 9038 9050
% Day3 3976 3988
% Day4 1165 1177 = 231 243
%----------------------------
% Rat 913
% Day1 626 643