#!/bin/bash
if [[ "$1" != "" ]];
then
	DIR="$1"
else
	echo "No arguments provided, current directory is selected."
	DIR=$(pwd)
fi
echo $DIR
non_raw_file=$(ls $DIR | grep "jpg")
raw_files=$(ls $DIR | grep "raw")
file_info=$(identify "${DIR}/${non_raw_file}")
size=$(echo ${file_info} | cut -d " " -f 3)
echo "Size = $size"
for f in $raw_files; do
	ffmpeg -f rawvideo -r 1 -s $size -pix_fmt gray -i "${DIR}/$f" "${DIR}/${f%.raw}.pgm"
done