#!/bin/bash

if [ $# -ne 3 ]; then
    echo
    echo "Usage: ./findHoles.sh [job ID] [expected # of files] [sched directory]"
    exit
fi

JOBID=$1
NFILES=$2
SCHEDDIRECTORY=$3

if [ "$NFILES" -le 0 ]; then
    echo "Argument can't be less than or equal to zero."
    echo
    exit
fi

FILENUM="0"
MISSINGINDICES=()
HOLECOUNTER=0

while [ $FILENUM -lt $NFILES ]; do
    FILE="${JOBID}_${FILENUM}.picoDst.result.root"

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
    
    cd ${SCHEDDIRECTORY}
    ${COMMAND} sched${JOBID}.session.xml
fi
