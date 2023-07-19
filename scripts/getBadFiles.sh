#!/bin/bash

shopt -s nullglob

#FILES=()
#while IFS=  read -r -d $'\0'; do
#    FILES+=("$REPLY")
#done < <(find ../ -type f -size 284k -print0)

JOBID=$1

FILES=(${JOBID}_*.picoDst.result.root)    #All of the first output files for any different jobs
NFILES=${#FILES[@]}

if [ $NFILES -eq 0 ]; then
    echo "No bad files found."
fi

for fname in "${FILES[@]}" ; do
    #FILENAME=${fname}    

    #tmp=${FULLPATH#*/}       # remove prefix ending in "/"
    IDwithINDEX=${fname%.picoDst*}   # remove suffix starting with ".picoDst"
    
    OUTFILE=out/${IDwithINDEX}.out

    #mapfile -t BADFILES < <(cat ${OUTFILE} | grep "=root://")   # fill BADFILES with the picoDsts that failed to copy
    mapfile -t BADFILES < <(grep "=root://" ${OUTFILE})   # fill BADFILES with the picoDsts that failed to copy
    NBADFILES=${#BADFILES[@]}

    for entry in "${BADFILES[@]}" ; do
	BADFILEPATH=${entry#*=}
	echo ${BADFILEPATH}
    done
done > allBadFiles.txt

