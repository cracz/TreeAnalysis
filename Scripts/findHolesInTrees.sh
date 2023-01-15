#!/bin/bash

if [ $# -ne 2 ]; then
    echo
    echo "Usage: ./findHoles.sh [job ID] [expected # of files]"
    #echo "Be sure to first remove all job IDs with the following command:"
    #echo "rename [job ID]_ \"\" *"
    exit
fi

JOBID=$1
NFILES=$2

if [ "$NFILES" -le 0 ]; then
    echo "Argument can't be less than or equal to zero."
    echo
    exit
fi

FILENUM="0"
MISSINGINDICES=()
HOLECOUNTER=0

while [ $FILENUM -lt $NFILES ]; do
    FILE="${JOBID}_${FILENUM}.root"

    if [ ! -e "$FILE" ]; then
	echo "File with index $FILENUM is missing!"
	MISSINGINDICES+=("${FILENUM}")
	HOLECOUNTER=$[$HOLECOUNTER+1]
    fi

    FILENUM=$[$FILENUM+1]
done

if [ "${#MISSINGINDICES[@]}" -eq 0 ]; then
    echo "No files missing."
    exit
else
    COMMAND="star-submit -r "
    COUNTER=1

    for i in "${MISSINGINDICES[@]}"
    do
	if [ "$COUNTER" -eq "${#MISSINGINDICES[@]}" ]; then
	    COMMAND+="${i}"
	else
	    COMMAND+="${i},"
	fi
	
	COUNTER=$[$COUNTER+1]
    done

    #echo
    #echo "Command for resubmission:"
    #echo
    #echo "${COMMAND} ${JOBID}.session.xml"

    echo
    echo ${HOLECOUNTER}" holes found in file list."
    echo
    read -p "Resubmit these jobs? [y||n] " choice

    if [ "$choice" = "y" ]; then
	cd ~/work/flow/
	${COMMAND} sched${JOBID}.session.xml
    fi
fi
