#!/bin/bash

shopt -s nullglob

JOBID=$1

FILES=(${JOBID}_*.root)    #All of the first output files for any different jobs
NFILES=${#FILES[@]}
CURRENTDIR=$(pwd)/

if [ $NFILES -eq 0 ]; then
    echo "No files found."

else
    for fname in "${FILES[@]}" ; do
	FILENAME=${fname}
	echo "${CURRENTDIR}${FILENAME}"
    done > file.list
    
    echo "File list saved to file.list"
fi
