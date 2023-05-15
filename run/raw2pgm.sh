#!/bin/bash
if [[ "$1" != "" ]];
then
	DIR="$1"
else
	echo "No arguments provided, current directory is selected."
	DIR=$(pwd)
fi

echo $DIR
if [[ "$2" == "" ]];
then
	non_raw_file=$(ls $DIR | grep "jpg")
	file_info=$(identify "${DIR}/${non_raw_file}")
	size=$(echo ${file_info} | cut -d " " -f 3)
else
	size="$2"
fi
echo "Size = $size"
raw_files=$(ls $DIR | grep "raw")
for f in $raw_files; do
	ffmpeg -f rawvideo -r 1 -s $size -pix_fmt gray -i "${DIR}/$f" "${DIR}/${f%.raw}.pgm"
done