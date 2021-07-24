#!/bin/bash
if [[ "$1" != "" ]];
then
	DIR="$1"
	echo $DIR
else
	echo "No arguments provided, current directory is selected."
	DIR=$(pwd)
fi
non_raw_files=(${DIR}/*.{jpg,pgm})
raw_files=(${DIR}/*.raw)
file_info=$(identify ${non_raw_files[1]})
size=$(echo ${file_info} | cut -d " " -f 3)
echo "Size = $size"
for f in $raw_files; do
	ffmpeg -f rawvideo -r 1 -s $size -pix_fmt gray -i "$f" "${f%.raw}.pgm"
done