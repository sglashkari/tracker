#!/bin/bash

# Check that two arguments have been provided
if [ $# -ne 2 ]; then
  echo "Usage: $0 video_file audio_file"
  exit 1
fi

# Get video length in seconds using ffprobe
video_length=$(ffprobe -v error -select_streams v:0 -show_entries stream=duration -of default=noprint_wrappers=1:nokey=1 "$1")

# Get audio length in seconds using ffprobe
audio_length=$(ffprobe -v error -select_streams a:0 -show_entries stream=duration -of default=noprint_wrappers=1:nokey=1 "$2")

# Calculate the length difference
length_diff=$(echo "$video_length - $audio_length" | bc)

# Print the results
echo "Video length: $video_length seconds"
echo "Audio length: $audio_length seconds"
echo "Length difference: $length_diff seconds"
