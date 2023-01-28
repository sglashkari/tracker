#!/bin/bash
if  [[ $# -eq 0 ]]
then
	file="/home/dome3tracking/Desktop/frames.csv"
else
	file="$PWD/times.csv"
fi
d=$(dirname $file)
text="$d/list.txt"
echo -e "#list of videos" > $text
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
		ffmpeg -i video.avi -vf "select=between(n\,$frame_initial\,$frame_final)" -vsync 0 "$d/frames-$frame_initial-$frame_final.mp4"
		echo -e "file '$d/frames-$frame_initial-$frame_final.mp4'" >> $text
	done
	ffmpeg -f concat -safe 0 -i $text -c copy "$d/concat.mp4"
	echo "$d/concat.mp4 is created!"
else
	ffmpeg -i video.avi -vf "select=between(n\,$1\,$2)" -vsync 0 "$d/frames-$1-$2.mp4"
	echo "$d/frames-$1-$2.mp4 created!"
fi