#!/bin/bash
if  [[ $# -eq 0 ]]
then
	file="/home/dome3tracking/Desktop/frames.csv"
else
	file=$1
fi
d=$(dirname $file)
echo $d
if  [[ $# -lt 2 ]]
then
	echo $file
#	sed 's/,/\t/g' $file > $file
	no_lines=$(cat $file | wc -l);
	echo $no_lines
	for ((i=1; i<=$no_lines; i++))
	do
		frame_initial=$(awk -v I=$i 'FNR == I {print $1}' $file | sed 's/,//')
		frame_final=$(awk -v I=$i 'FNR == I {print $2}' $file | sed 's/,//')
		frame_final=${frame_final:0:${#frame_final}-1}
		echo "line $i: from frame no $frame_initial to frame no $frame_final, in total $((frame_final-frame_initial)) frames"
		#sleep 5
		ffmpeg -start_number $frame_initial -framerate 30 -i frame-%d.pgm -crf 24 -pix_fmt gray -preset veryslow -r 30 -vframes $((frame_final-frame_initial)) "$d/frames-$frame_initial-$frame_final.mp4"
	done
else
	ffmpeg -start_number $1 -framerate 30 -i frame-%d.pgm -crf 24 -pix_fmt gray -preset veryslow -r 30 -vframes $(($2-$1)) "$d/frames-$1-$2.mp4"
fi