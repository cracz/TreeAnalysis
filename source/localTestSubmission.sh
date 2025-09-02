#!/bin/bash

starver SL23d
#TreeAnalyzer /star/data01/pwg/cracz/Data_3p0GeV_FXT/FXT_3p0GeV_SL20d_2018_212.root TESTJOB ../Configs/fxt_3p0GeV/config_3p0GeV.txt ../CorrectionFiles/fxt_3p0GeV/correctionInfo_INPUT.root ../CorrectionFiles/fxt_3p0GeV/resolutionInfo_INPUT_3p0GeV_averagedRes.root

#TreeAnalyzer /star/data01/pwg/cracz/Data_3p5GeV_FXT/3EF0FC8F160EA89C7EC83BA7CF144481_0.root TESTJOB ../Configs/fxt_3p5GeV/config_3p5GeV.txt ../CorrectionFiles/fxt_3p5GeV/correctionInfo_INPUT_3p5GeV.root ../CorrectionFiles/fxt_3p5GeV/resolutionInfo_INPUT_3p5GeV_1to6_8to11_NewResolutions.root

valgrind --tool=memcheck --leak-check=yes --log-file=Valgrind.log --num-callers=30 --suppressions=$ROOTSYS/etc/valgrind-root.supp TreeAnalyzer /star/data01/pwg/cracz/Data_3p5GeV_FXT/3EF0FC8F160EA89C7EC83BA7CF144481_0.root TESTJOB ../Configs/fxt_3p5GeV/config_3p5GeV.txt ../CorrectionFiles/fxt_3p5GeV/correctionInfo_INPUT_3p5GeV.root ../CorrectionFiles/fxt_3p5GeV/resolutionInfo_INPUT_3p5GeV_1to6_8to11_NewResolutions.root
