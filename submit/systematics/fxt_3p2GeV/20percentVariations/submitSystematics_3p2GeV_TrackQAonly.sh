#!/bin/bash

nohup star-submit ../../../treeAnalyze_3p2GeV.xml >& submitOut_3p2GeV_Normal.txt &

nohup star-submit treeAnalyze_3p2GeV_dca_low.xml >& submitOut_3p2GeV_dca_low.txt &
nohup star-submit treeAnalyze_3p2GeV_dca_high.xml >& submitOut_3p2GeV_dca_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nhits_low.xml >& submitOut_3p2GeV_nhits_low.txt &
nohup star-submit treeAnalyze_3p2GeV_nhits_high.xml >& submitOut_3p2GeV_nhits_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nhitsdedx_high.xml >& submitOut_3p2GeV_nhitsdedx_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nhitsratio_low.xml >& submitOut_3p2GeV_nhitsratio_low.txt &
nohup star-submit treeAnalyze_3p2GeV_nhitsratio_high.xml >& submitOut_3p2GeV_nhitsratio_high.txt &
