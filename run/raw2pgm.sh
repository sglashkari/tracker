#!/bin/bash
files=(*.{jpg,pgm})
non_raw_file=$(echo "${files[1]}")
file_info=$(identify ${non_raw_file})
size=$(echo ${file_info} | cut -d " " -f 3)
echo "Size = $size"
for f in *.raw; do
	ffmpeg -f rawvideo -r 1 -s $size -pix_fmt gray -i "$f" "${f%.raw}.pgm"
done