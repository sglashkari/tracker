#!/bin/bash

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <video.avi> [-n <start_frame>,<end_frame>] [<frames.csv>]"
    exit 1
else
    video="$1"
    if [[ $2 == "-n" ]]; then
        if [[ $# -lt 3 ]]; then
            echo "Usage: $0 <video.avi> [-n <start_frame>,<end_frame>] [<frames.csv>]"
            exit 1
        else
            start_frame=$(echo "$3" | awk -F ',' '{print $1}')
            end_frame=$(echo "$3" | awk -F ',' '{print $2}')
            no_lines=1
        fi
        file=""
    else
        file="$2"
        if [[ ! -f $file ]]; then
            echo "$file not found"
            exit 1
        fi
        no_lines=$(cat "$file" | wc -l);
        start_frame=""
        end_frame=""
    fi
fi

if [[ ! -f $video ]]; then
    echo "$video not found"
    exit 1
fi

if [[ $no_lines -lt 1 ]]; then
    echo "No valid frame ranges found in $file"
    exit 1
fi

d=$(dirname "$video")
text="$d/list.txt"
echo -e "#list of videos" > "$text"

for ((i=1; i<=$no_lines; i++))
do
    if [[ $no_lines -eq 1 ]]; then
        frame_initial=$start_frame
        frame_final=$end_frame
    else
        frame_initial=$(awk -v I=$i 'FNR == I {print $1}' "$file" | sed 's/,//')
        frame_final=$(awk -v I=$i 'FNR == I {print $2}' "$file" | sed 's/,//')
    fi
    echo "Processing row $i: From frame no $frame_initial to frame no $frame_final, in total $((frame_final-frame_initial)) frames"
    ffmpeg -i "$video" -vf "select=between(n\,$frame_initial\,$frame_final)" -vsync 0 "$d/frames-$frame_initial-$frame_final.mp4"
    echo -e "file '$d/frames-$frame_initial-$frame_final.mp4'" >> "$text"
done

if [[ $no_lines -ne 1 ]]; then
    ffmpeg -f concat -safe 0 -i "$text" -c copy "$d/concat.mp4"
    echo "$d/concat.mp4 is created!"
fi
