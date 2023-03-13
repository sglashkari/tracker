#!/bin/bash

# Check if at least two arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <video_file> <audio_file> [output_file]"
    exit 1
fi

# Set variables for the input files
VIDEO_FILE="$1"
AUDIO_FILE="$2"

# Get the full path of the video file
AUDIO_PATH="$(realpath "${AUDIO_FILE%.*}")"

# Get the base name of the audio file (i.e. remove the path and extension)
VIDEO_BASENAME="$(basename "${VIDEO_FILE%.*}")"

# Set the output filename based on the video and audio file names, or use the third argument if provided
if [ -z "$3" ]; then
    OUTPUT_FILE="$AUDIO_PATH"_"$VIDEO_BASENAME".mp4
else
    OUTPUT_FILE="$3"
fi

# Combine the video and audio files using ffmpeg
ffmpeg -i "$VIDEO_FILE" -i "$AUDIO_FILE" -c:v copy -c:a aac -strict experimental "$OUTPUT_FILE"

# Check if ffmpeg was successful
if [ $? -eq 0 ]; then
    echo "Successfully combined $VIDEO_FILE and $AUDIO_FILE into $OUTPUT_FILE"
else
    echo "Failed to combine $VIDEO_FILE and $AUDIO_FILE"
fi
