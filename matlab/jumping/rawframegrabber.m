%% Reads a pgm frame by frame
% Updated Nov 15, 2020
% detecting the LED light on/off
% Author Shahin G Lashkari
close all
default_path = '/home/dome3tracking/Videos/';%'D:\OneDrive - Johns Hopkins\JHU\913_Jumping_Recording\2020-10-25_Rat913-01\Videos\2020-10-25_15-59-09.mp4'; %day 1
% try
path = uigetdir(default_path);


frame_no = input('What is the frame number? ');

figure;
while true
    dim = [0.15 0.75 0.1 0.1];
    
    row=740;  col=1536;
    fin=fopen('b1.raw','r');
    I=fread(fin,row*col,'uint8=>uint8'); 
    Z=reshape(I,row,col);
    Z=Z';
    k=imshow(Z)

    filename = ['frame-' num2str(frame_no) '.raw'];
    
    frame = imread(fullfile(path,filename));
    ax = imshow(frame);
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
imshow(frame);