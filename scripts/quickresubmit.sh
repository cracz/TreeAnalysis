#!/bin/bash

shopt -s nullglob

echo
echo "What energy is this?"
echo "[0] Cancel"
echo "[1] 3.0 GeV"
echo "[2] 3.2 GeV"
echo "[3] 3.5 GeV"
echo "[4] 3.9 GeV"
echo "[5] 4.5 GeV"
echo

read -p "Enter an integer: " CHOICE

if [[ -n ${CHOICE//[0-9]/} ]]; then
    echo "Please input only an integer."
    exit 1
elif [[ ${CHOICE} -eq 0 ]]; then
    echo "Cancelling."
    exit 0
fi

echo
read -p "Enter the job ID to clean up: " JOBID

CLEANOLD="rm "${JOBID}"_* correctionInfo_OUTPUT_"${JOBID}"_* out/"${JOBID}"_*"

MOVE="cd /star/u/cracz/work/flow/TreeAnalysis/submit/"
CLEANSCHED="rm -r sched"${JOBID}"*"

SUBMIT1="star-submit treeAnalyze"
SUBMIT2=".xml"

if [[ ${CHOICE} -eq 1 ]]; then
    SUBMIT2=".xml"
elif [[ ${CHOICE} -eq 2 ]]; then
    SUBMIT2="_3p2GeV.xml"
elif [[ ${CHOICE} -eq 3 ]]; then
    SUBMIT2="_3p5GeV.xml"
elif [[ ${CHOICE} -eq 4 ]]; then
    SUBMIT2="_3p9GeV.xml"
elif [[ ${CHOICE} -eq 5 ]]; then
    SUBMIT2="_4p5GeV.xml"
fi

echo
echo "Cleaning old results..."
echo "${CLEANOLD}"
${CLEANOLD}

echo
echo "Changing to submit directory and cleaning old scheduler files..."
echo "${MOVE}"
echo "${CLEANSCHED}"
${MOVE}
${CLEANSCHED}

echo
echo "Resubmitting jobs..."
echo "${SUBMIT1}${SUBMIT2}"
${SUBMIT1}${SUBMIT2}
