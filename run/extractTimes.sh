#!/bin/bash

set -o noclobber
max_no=$(ls | sed -e '/raw/!d' | wc -l);
timefile=times.csv;
if [ ! -f "$timefile" ]; then
	ls | sed -e '/raw/!d' -e 's/_/, /g' -e 's/.raw//g' | awk '{printf "%2d, %2d, %2d, %3d, %8.3f\n", $1, $2, $3, $4, $1 * 3600 + $2 * 60+ $3 + $4 * 0.001}' > times.csv;
	echo "$timefile is created!"
else
    if [[ $(wc -l < "$timefile") -eq $max_no ]]; then
    	echo "$timefile already exists!"
    else
        echo "$timefile has less elements ($(wc -l < "$timefile")) than the number of raw files ($max_no)!"
        exit 1
    fi
fi
set +o noclobber
cat "$timefile" | (head;tail)