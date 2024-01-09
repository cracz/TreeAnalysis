#!/bin/bash

starver SL23d
#TreeAnalyzer /star/data01/pwg/cracz/Data_3p0GeV_FXT/FXT_3p0GeV_SL20d_2018_212.root TESTJOB ../Configs/fxt_3p0GeV/config_3p0GeV.txt ../CorrectionFiles/fxt_3p0GeV/correctionInfo_INPUT.root ../CorrectionFiles/fxt_3p0GeV/resolutionInfo_INPUT_3p0GeV_averagedRes.root

TreeAnalyzer /star/data01/pwg/cracz/Data_3p2GeV_FXT/D12984BD22C3810D140327578CFA7483_0.root TESTJOB ../Configs/fxt_3p2GeV/config_3p2GeV.txt ../CorrectionFiles/fxt_3p2GeV/correctionInfo_INPUT_3p2GeV.root ../CorrectionFiles/fxt_3p2GeV/resolutionInfo_INPUT_3p2GeV_1to6_max13_averagedRes_SL23d.root

