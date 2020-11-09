#!/bin/bash

used_webcams=$(fuser /dev/video0 | wc -w);
if [ $used_webcams -gt 0 ]
then
	echo $used_webcams 'webcam(s) is currently not available!'
	exit
fi

dt=$(date +%Y-%m-%d_%H-%M-%S);
echo $dt

mkdir ~/Videos/$dt;

gnome-terminal -- bash -c "cd ~/Videos/${dt}; ~/flycapture/bin/FlyCapture2Test"

# Top webcam
top_cam_available=$(ls /dev/video* | grep video0 | wc -l);
if [ $top_cam_available -eq 1 ]; then
	# settings
	v4l2-ctl \
	--device=0 \
	--set-ctrl=brightness=180 \
	--set-ctrl=contrast=128 \
	--set-ctrl=saturation=225 \
	--set-ctrl=white_balance_temperature_auto=1 \
	--set-ctrl=gain=255 \
	--set-ctrl=power_line_frequency=2 \
	--set-ctrl=white_balance_temperature=5000 \
	--set-ctrl=sharpness=130 \
	--set-ctrl=backlight_compensation=1 \
	--set-ctrl=exposure_auto=0 \
	--set-ctrl=exposure_absolute=350 \
	--set-ctrl=exposure_auto_priority=0 \
	--set-ctrl=pan_absolute=-16912 \
	--set-ctrl=tilt_absolute=-29721 \
	--set-ctrl=focus_absolute=0 \
	--set-ctrl=focus_auto=0 \
	--set-ctrl=zoom_absolute=130 \
	--set-ctrl=led1_mode=1
	gnome-terminal -- bash -c "ffmpeg -f v4l2 -framerate 30 -video_size 1280x720 -input_format mjpeg -i /dev/video0 ~/Videos/${dt}/${dt}_top.mkv"
fi

# Side webcam
side_cam_available=$(ls /dev/video* | grep video2 | wc -l);
if [ $side_cam_available -eq 1 ]; then

	# settings
	v4l2-ctl \
	--device=2 \
	--set-ctrl=brightness=200 \
	--set-ctrl=contrast=128 \
	--set-ctrl=saturation=128 \
	--set-ctrl=white_balance_temperature_auto=1 \
	--set-ctrl=gain=255 \
	--set-ctrl=power_line_frequency=2 \
	--set-ctrl=white_balance_temperature=5000 \
	--set-ctrl=sharpness=200 \
	--set-ctrl=backlight_compensation=0 \
	--set-ctrl=exposure_auto=0 \
	--set-ctrl=exposure_absolute=1000 \
	--set-ctrl=exposure_auto_priority=0 \
	--set-ctrl=focus_absolute=0 \
	--set-ctrl=focus_auto=0 \
	--set-ctrl=zoom_absolute=150 \
	--set-ctrl=led1_mode=1 

	gnome-terminal -- bash -c "ffmpeg -f v4l2 -framerate 30 -video_size 1280x720 -input_format mjpeg -i /dev/video2 ~/Videos/${dt}/${dt}_side.mkv"
fi