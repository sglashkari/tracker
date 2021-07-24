#!/bin/bash
if [[ "$1" != "" ]];
then
	DIR="$1"
else
	echo "No arguments provided, current directory is selected."
	DIR=$(pwd)
fi
files="$DIR*.{jpg,pgm}"
file_info=$(identify ${files[1]})
size=$(echo ${file_info} | cut -d " " -f 3)
echo "Size = $size"
for f in "$DIR*.raw"; do
	ffmpeg -f rawvideo -r 1 -s $size -pix_fmt gray -i "$f" "${f%.raw}.pgm"
done