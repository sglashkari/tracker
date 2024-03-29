#!/bin/bash
~/test/run/raw2order.sh

if  [[ $# -eq 0 ]]
then
	size="2048x400"
	
else
	size=$1
fi

echo "size is $size"
sleep 2

if  [[$# -eq 1]]; then
	size="2048x400"
else
	size=$1
fi

max_no=$(ls | sed -e '/raw/!d' | wc -l);
max_no=$((max_no-1))
counter=000

lastfilename="frame-$((max_no-1)).pgm"
if [ -f "$lastfilename" ]; then
    echo "$lastfilename already exists!"
    ls -t | sed -e '/pgm/!d' | (head;tail)
    echo "The pgm files are already available!"
else
	for f in *.raw; do
		ffmpeg -f rawvideo -r 1 -s $size -pix_fmt gray -i "$f" "${f%.raw}.pgm" && ((counter++))
		printf '%d%%\n' $((counter*100/max_no))
	done
fi

ffmpeg -framerate 30 -i frame-%d.pgm -crf 0 -pix_fmt gray -preset veryslow -r 30 lossless-video.mp4

pose="../pose.csv"
if [ ! -f "$pose" ]; then
	touch "$pose"
	echo "$pose created!"
fi

du -h lossless-video.mp4