#!/bin/bash
~/test/run/raw2order.sh

if  [[ $# -eq 0 ]]
then
	size="1200x700"
	
else
	size=$1
fi

echo "size is $size"
sleep 2

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

ffmpeg -framerate 20 -i frame-%d.pgm -crf 18 -pix_fmt gray -preset veryslow -r 20 side_video_crf18.mp4

du -h side_video_crf18.mp4