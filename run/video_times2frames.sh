#!/bin/bash

# Check if the required argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_video_file> [<input_times_file>] [<output_file>]"
    exit 1
fi

# Set the input video file name
input_file="$1"

# Get the video duration
duration=$(ffmpeg -i "$input_file" 2>&1 | grep "Duration" | cut -d ' ' -f 4 | sed s/,//)

# Extract the frame rate
frame_rate=$(ffmpeg -i "$input_file" 2>&1 | grep -oP "(?<=, )\d+(?=\s+fps)")

# Check if the frame rate is found
if [ -z "$frame_rate" ]; then
    echo "Frame rate not found. Exiting."
    exit 1
fi

# Set the input times file name
if [ $# -ge 2 ]; then
    input_times_file="$2"
else
    input_times_file="video_times.csv"
fi

# Set the output file name
if [ $# -ge 3 ]; then
    output_file="$3"
else
    output_file="video_frames.csv"
fi

# Read the input file and extract frame numbers
while IFS=',' read -r start_time end_time || [ -n "$start_time" ]; do
    # Convert the start time to seconds
    start_seconds=$(echo "$start_time" | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }')
    # Convert the end time to seconds
    end_seconds=$(echo "$end_time" | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }')
    # Calculate the start frame number corresponding to the start time
    start_frame_number=$(echo "($start_seconds)*$frame_rate + 1" | bc)
    # Calculate the end frame number corresponding to the end time
    end_frame_number=$(echo "($end_seconds)*$frame_rate + 1" | bc)
    # Add the start time, end time, start frame number and end frame number to the output file
    echo "$start_frame_number, $end_frame_number" >> "$output_file"
done < "$input_times_file"

# Print the output file name
echo "Frame numbers saved to $output_file."
