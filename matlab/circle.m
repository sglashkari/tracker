clc; close all;
vid_filename = 'D:\domeLatencyMeasurement\200704_latencyTest\200704_latencyTest_45fps_noRecording_video.avi';
videoReader  = VideoReader(vid_filename)
videoWriter = VideoWriter(insertBefore(vid_filename,".","_tracked"));
num_frames = videoReader.NumFrame; % 1683
frame_rate = videoReader.FrameRate;
videoWriter.FrameRate = frame_rate;
videoWriter.Quality = 100;
c=zeros(num_frames,4);
open(videoWriter);
h = figure;
tic
for frame_no = 1:num_frames%videoReader.NumFrames
    % Read image (i.e. a frame) from a video
    I = read(videoReader,frame_no);
    
    hold off
    ax = imshow(I);
    hold on
    
    %% Projected image
    % ref: https://www.mathworks.com/matlabcentral/answers/495720-how-to-fit-an-ellipse-to-an-image-in-matlab#answer_405863
    xmin = 0;
    ymin = 450;
    I2 = imcrop(I,[xmin ymin 900 400]);
    
    % Binarize
    Igray = rgb2gray(I2);
    BW = imbinarize(Igray);
    
    % Extract the maximum area
    BW = imclearborder(BW,8);
    BW = bwareafilt(BW,1);
    
    % Calculate centroid, orientation and major/minor axis length of the ellipse
    s = regionprops('table', BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
    x_proj = s.Centroid(1)+xmin;
    y_proj = s.Centroid(2)+ymin;
    
    %% crown
    [centers, radii, metric] = imfindcircles(I,[11 20]);
    
    % picking one cirle
    if size(centers,1)>1
        fprintf('check frame %d\n', frame_no)
        if abs(centers(2,1)-centers(1,1))<100
            if centers(2,1)>centers(1,1)
                centers = centers(2,:);
            else
                centers = centers(1,:);
            end
        else
            if centers(2,1)>centers(1,1)
                centers = centers(1,:);
            else
                centers = centers(2,:);
            end
        end
    end
    
    if isempty(centers)
        centers = [nan nan];
    end
    x_crown = centers(1);
    y_crown = centers(2);
    
    plot(x_proj, y_proj, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    plot(x_crown, y_crown, 'b+', 'MarkerSize', 10, 'LineWidth', 2);
    
    if mod(frame_no,100)==0
        disp(frame_no);
        toc
    end
    
    frame = getframe(gcf);
    writeVideo(videoWriter,frame);
    pause(1/frame_rate);
    c(frame_no,:) = [x_proj, y_proj, x_crown, y_crown];
end

close(videoWriter)

k=1:length(c);
t=(k-1)/frame_rate;

save(strrep(vid_filename,'avi','mat'),'c','t','frame_rate')