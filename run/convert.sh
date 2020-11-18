#!/bin/bash

counter=000
for f in *.raw; do 
    mv $f frame-$((counter)).raw && ((counter++))	# mv -- "$f" "${f%.raw}.png"
done

:'
for f in *.raw; do
	#mv -- "$f" "${f%.raw}.png"
	ffmpeg -f rawvideo -r 1 -s 1536x740 -pix_fmt gray -i "$f" "${f%.raw}.pgm"
done

#ffmpeg -f rawvideo -s 1536x740 -pix_fmt gray -i frame-1.raw frame-1.pgm

#ffmpeg -framerate 240 -f rawvideo -s 1536x740 -i frame-%d.raw -vcodec libx264 -crf 0 -pix_fmt gray -r 240 test.mp4
ffmpeg -framerate 240 -i frame-%d.pgm -c:v libx264 -crf 0 -r 240 output.avi



# ffmpeg -framerate 240 -s 1536x740 -pix_fmt gray -i frame-%d.raw -f rawvideo -c:v libx264 -crf 0 -r 240 output.avi
'
ffmpeg -framerate 240 -s 1536x740 -pix_fmt gray -i frame-%d.raw -f rawvideo -c:v libx264 -crf 0 -r 240 uncompressed.avi \
																			-c:v libx264 -crf 20 -r 240 compressed.avi