#!/bin/bash
# check if directory is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 directory [-p|-r]"
    exit 1
fi
# set directory
directory="$1"
# default values
pattern="raw"
# check if pattern is provided
if [ $# -gt 1 ]; then
    if [ "$2" == "-p" ]; then
        pattern="pgm"
    elif [ "$2" != "-r" ]; then
        echo "Invalid pattern specified. Please specify -r for 'raw' or -p for 'pgm'."
        exit 1
    fi
fi
# display warning message
echo "WARNING: This script will delete all files containing the word '$pattern' in the directory '$directory'."
echo "Are you sure you want to continue? (y/n)"
# read user input
read confirm
# if user confirms, run the command
if [ "$confirm" == "y" ]; then
    if [ "$directory" != "." ] && [ "$directory" != "./" ]; then
        cd "$directory" || exit 1 # exit if the directory does not exist or cannot be accessed
    fi
    ls -h | sed -e "/$pattern/!d" | xargs rm
else
    echo "Aborted."
fi