#!/bin/bash
:'
counter=000
for f in *.raw; do 
    mv $f frame-$((counter)).raw && ((counter++))	# mv -- "$f" "${f%.raw}.png"
done
'
:'
for f in frame-{590154..591351}.raw; do
	#mv -- "$f" "${f%.raw}.png"
	ffmpeg -f rawvideo -r 1 -s 1600x580 -pix_fmt gray -i "$f" "${f%.raw}.pgm" # 1536x740, 1600x580 (Day 3), 1440x708, 
	cp $f frame-$((counter)).raw && ((counter++))	# mv -- "$f" "${f%.raw}.png"
done

'
counter=000
for f in frame-{0..1000}.raw; do
	ffmpeg -f rawvideo -r 1 -s 2048x400 -pix_fmt gray -i "$f" frame-$((counter)).pgm # 1536x740, 1600x580 (Day 3), 1440x708, 2048x400 (2021-11-19 top), 1920x600
	#ffmpeg -i input -vf "scale='min(1280,iw)':min'(720,ih)':force_original_aspect_ratio=decrease,pad=1280:720:-1:-1:color=black" output

    ((counter++))	# mv -- "$f" "${f%.raw}.png"
done
#ffmpeg -f rawvideo -s 1536x740 -pix_fmt gray -i frame-1.raw frame-1.pgm

#ffmpeg -framerate 240 -f rawvideo -s 1536x740 -i frame-%d.raw -vcodec libx264 -crf 0 -pix_fmt gray -r 240 test.mp4
#ffmpeg -framerate 30 -i frame-%d.pgm -c:v libx264 -crf 0 -r 30 output.avi



#ffmpeg -framerate 240 -s 1536x740 -pix_fmt gray -i frame-%d.raw -f rawvideo -c:v libx264 -crf 0 -r 240 output.avi

#ffmpeg -framerate 250 -s 1440x708 -pix_fmt gray -i frame-%d.raw -f rawvideo -c:v libx264 -crf 0 -r 250 uncompressed.avi \
#																			-c:v libx264 -crf 20 -r 250 compressed.avi

#ffmpeg -framerate 300 -s 1600x580 -pix_fmt gray -i frame-%d.raw -f rawvideo -c:v libx264 -crf 0 -r 300 uncompressed.avi \
#																			-c:v libx264 -crf 20 -r 300 compressed.avi


#ffmpeg -framerate 300 -s 1600x580 -i frame-%d.ppm -i audio.wav -f rawvideo -c:v libx264 -crf 0 -r 300 uncompressed.avi \
#																			-c:v libx264 -crf 20 -r 300 compressed.avi

#ffmpeg -framerate 30 -i frame-%d.pgm -c:v libx264 -crf 0 -r 30 uncompressed.avi \
#									 -c:v libx264 -crf 20 -r 30 compressed.avi

ffmpeg -framerate 30 -i frame-%d.pgm -crf 18 -r 30 output.mp4
ffmpeg -i output.mp4 -i background.jpg -filter_complex "[0:v]scale=2048:2048:force_original_aspect_ratio=decrease[fg];[1][fg]overlay=0:976" output2.mp4


#rm frame-{0..1200}.raw
'