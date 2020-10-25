%% Reads a video frame by frame
% Updated March 28, 2020, modified Oct 22, 2020
% detecting the Red LED light
% Author Shahin G Lashkari
close all
%vid_filename = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Videos\2020-03-29_20-22-00.mp4';
vid_filename = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\913_Jumping_Recording\2020-10-23-Day0\2020-10-23-133050.mp4';
videoReader  = VideoReader(vid_filename)

% tic
% S = zeros(1e4,1);
% for frame_no=1:length(S)
%     R = read(videoReader,frame_no);
%     S(frame_no)=sum(sum(R(:,:,1)));
%     Z(:,:,:,frame_no) = R;
% end
% plot(1:length(S),S);
% [~,argmax] = max(S);
% toc
% frame_no = argmax - 10;

frame_no = 15;

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
    
    button = double(get(gcf,'CurrentCharacter'))
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


% Day1 530 542
% Day2 9038 9050
% Day3 3976 3988
% Day4 1165 1177 = 231 243