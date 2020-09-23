%% Reads a video frame by frame
% Updated March 28, 2020
% detecting the Red LED light
% Author Shahin G Lashkari
close all
vid_filename = 'C:\Users\Shahin\OneDrive - Johns Hopkins University\JHU\883_Jumping_Recording\200329_Rat883-04\Videos\2020-03-29_20-22-00.mp4';
videoReader  = VideoReader(vid_filename);

tic
S = zeros(1e4,1);
for frame_no=1:length(S)
    R = read(videoReader,frame_no);
    S(frame_no)=sum(sum(R(:,:,1)));
    Z(:,:,:,frame_no) = R;
end
plot(1:length(S),S);
[~,argmax] = max(S);
toc
frame_no = argmax - 10;
    
% frame = 232;

figure;
while true
    dim = [0.15 0.75 0.1 0.1];
    frames = read(videoReader,frame_no);
    ax = imshow(frames);
    fig = ancestor(ax, 'figure');
    annotation('textbox',dim,'String',num2str(frame_no),'FitBoxToText','on','Color','white');
    
    try
        w = waitforbuttonpress;
    catch
        break;
    end
    
    sel = get(fig, 'SelectionType');
    if strcmpi(sel, 'alt')
        frame_no = frame_no - 1;
    else
        frame_no = frame_no + 1;
    end
    
end


% Day1 530 542
% Day2 9038 9050
% Day3 3976 3988
% Day4 1165 1177 = 231 243