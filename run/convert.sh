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
for f in frame-{605734..606931}.raw; do
	#ffmpeg -f rawvideo -r 1 -s 1600x580 -pix_fmt gray -i "$f" frame-$((counter)).ppm # 1536x740, 1600x580 (Day 3), 1440x708, 
    ((counter++))	# mv -- "$f" "${f%.raw}.png"
done
#ffmpeg -f rawvideo -s 1536x740 -pix_fmt gray -i frame-1.raw frame-1.pgm

#ffmpeg -framerate 240 -f rawvideo -s 1536x740 -i frame-%d.raw -vcodec libx264 -crf 0 -pix_fmt gray -r 240 test.mp4
#ffmpeg -framerate 240 -i frame-%d.pgm -c:v libx264 -crf 0 -r 240 output.avi



#ffmpeg -framerate 240 -s 1536x740 -pix_fmt gray -i frame-%d.raw -f rawvideo -c:v libx264 -crf 0 -r 240 output.avi

#ffmpeg -framerate 250 -s 1440x708 -pix_fmt gray -i frame-%d.raw -f rawvideo -c:v libx264 -crf 0 -r 250 uncompressed.avi \
#																			-c:v libx264 -crf 20 -r 250 compressed.avi

#ffmpeg -framerate 300 -s 1600x580 -pix_fmt gray -i frame-%d.raw -f rawvideo -c:v libx264 -crf 0 -r 300 uncompressed.avi \
#																			-c:v libx264 -crf 20 -r 300 compressed.avi


ffmpeg -framerate 300 -s 1600x580 -i frame-%d.ppm -i audio.wav -f rawvideo -c:v libx264 -crf 0 -r 300 uncompressed.avi \
																			-c:v libx264 -crf 20 -r 300 compressed.avi

#rm frame-{0..1200}.raw
'