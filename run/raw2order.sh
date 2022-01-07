#!/bin/bash
~/test/run/extractTimes.sh
max_no=$(ls | sed -e '/raw/!d' | wc -l);
lastfilename="frame-$((max_no-1)).raw"
timefile=times.csv;

if [ -f "$lastfilename" ]; then
    echo "$lastfilename already exists!"
    ls -t | sed -e '/raw/!d' | (head;tail)
    echo "The raw files are already in order!"
else
    if [[ $(wc -l < "$timefile") -eq $max_no ]]; then
        counter=000
        printf 'Creating frames in order:\n'
        for f in *.raw; do 
            mv $f frame-$((counter)).raw && ((counter++))
            if ! ((counter % 1000)); then
                printf '\r%d%%' $((counter*100/max_no))
            fi
        done
        printf '\rDone!\n'
    else
        echo "The number of elements in $timefile ($(wc -l < "$timefile")) is different from the number of raw files ($max_no)!"
    fi
fi