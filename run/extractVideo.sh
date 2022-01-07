#!/bin/bash
~/test/run/raw2order.sh
sleep 2

max_no=$(ls | sed -e '/raw/!d' | wc -l);
max_no=$((max_no-1))
counter=000
for f in *.raw; do
	ffmpeg -f rawvideo -r 1 -s 2048x400 -pix_fmt gray -i "$f" "${f%.raw}.pgm" && ((counter++))
	printf '%d%%\n' $((counter*100/max_no))
done

ffmpeg -framerate 30 -i frame-%d.pgm -crf 0 -pix_fmt gray -r 30 "lossless-video.mp4"

pose="../pose.csv"
if [ ! -f "$pose" ]; then
	touch "$pose"
	echo "$pose created!"
fi