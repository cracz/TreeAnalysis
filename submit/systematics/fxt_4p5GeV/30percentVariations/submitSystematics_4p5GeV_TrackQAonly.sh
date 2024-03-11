#!/bin/bash

#nohup star-submit ../../../treeAnalyze_4p5GeV.xml >& submitOut_4p5GeV_Normal.txt &

nohup star-submit treeAnalyze_4p5GeV_dca_low.xml >& submitOut_4p5GeV_dca_low.txt &
nohup star-submit treeAnalyze_4p5GeV_dca_high.xml >& submitOut_4p5GeV_dca_high.txt &

nohup star-submit treeAnalyze_4p5GeV_nhits_low.xml >& submitOut_4p5GeV_nhits_low.txt &
nohup star-submit treeAnalyze_4p5GeV_nhits_high.xml >& submitOut_4p5GeV_nhits_high.txt &

nohup star-submit treeAnalyze_4p5GeV_nhitsdedx_high.xml >& submitOut_4p5GeV_nhitsdedx_high.txt &

nohup star-submit treeAnalyze_4p5GeV_nhitsratio_low.xml >& submitOut_4p5GeV_nhitsratio_low.txt &
nohup star-submit treeAnalyze_4p5GeV_nhitsratio_high.xml >& submitOut_4p5GeV_nhitsratio_high.txt &
