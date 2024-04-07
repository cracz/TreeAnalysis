#!/bin/bash

shopt -s nullglob

FILES=(correctionInfo_OUTPUT_*.root)    #All of the output correction files for any different jobs
NFILES=${#FILES[@]}

if [ $NFILES -eq 0 ]; then
    echo "No correction files found."
    exit 0
fi

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

CORRECTIONFILEDIR="~/work/flow/TreeAnalysis/CorrectionFiles/"
SUBDIR="fxt_3p0GeV/"
NAMEPART1="correctionInfo_INPUT"
NAMEPART2=".root"

#echo ${CHOICE//[0-9]/}

if [[ -n ${CHOICE//[0-9]/} ]]; then
    echo "Please input only an integer."
    exit 1
elif [[ ${CHOICE} -eq 0 ]]; then
    echo "Cancelling."
    exit 0
elif [[ ${CHOICE} -eq 1 ]]; then
    NAMEPART2=".root"
    SUBDIR="fxt_3p0GeV/"
elif [[ ${CHOICE} -eq 2 ]]; then
    NAMEPART2="_3p2GeV.root"
    SUBDIR="fxt_3p2GeV/"
elif [[ ${CHOICE} -eq 3 ]]; then
    NAMEPART2="_3p5GeV.root"
    SUBDIR="fxt_3p5GeV/"
elif [[ ${CHOICE} -eq 4 ]]; then
    NAMEPART2="_3p9GeV.root"
    SUBDIR="fxt_3p9GeV/"
elif [[ ${CHOICE} -eq 5 ]]; then
    NAMEPART2="_4p5GeV.root"
    SUBDIR="fxt_4p5GeV/"
fi

hadd -f ${CORRECTIONFILEDIR}${SUBDIR}${NAMEPART1}${NAMEPART2} correctionInfo_OUTPUT_*.root

echo
echo "Correction information saved to this destination path and file:"
echo "${CORRECTIONFILEDIR}${SUBDIR}${NAMEPART1}${NAMEPART2}"
echo
