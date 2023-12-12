#!/bin/bash

#nohup star-submit ../../../treeAnalyze_3p2GeV.xml >& submitOut_3p2GeV_Normal.txt &

nohup star-submit treeAnalyze_3p2GeV_dca_low.xml >& submitOut_3p2GeV_dca_low.txt &
nohup star-submit treeAnalyze_3p2GeV_dca_high.xml >& submitOut_3p2GeV_dca_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nhits_low.xml >& submitOut_3p2GeV_nhits_low.txt &
nohup star-submit treeAnalyze_3p2GeV_nhits_high.xml >& submitOut_3p2GeV_nhits_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nhitsdedx_high.xml >& submitOut_3p2GeV_nhitsdedx_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nhitsratio_low.xml >& submitOut_3p2GeV_nhitsratio_low.txt &
nohup star-submit treeAnalyze_3p2GeV_nhitsratio_high.xml >& submitOut_3p2GeV_nhitsratio_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nSigPi_low.xml >& submitOut_3p2GeV_nSigPi_low.txt &
nohup star-submit treeAnalyze_3p2GeV_nSigPi_high.xml >& submitOut_3p2GeV_nSigPi_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nSigKa_low.xml >& submitOut_3p2GeV_nSigKa_low.txt &
nohup star-submit treeAnalyze_3p2GeV_nSigKa_high.xml >& submitOut_3p2GeV_nSigKa_high.txt &

nohup star-submit treeAnalyze_3p2GeV_nSigPr_low.xml >& submitOut_3p2GeV_nSigPr_low.txt &
nohup star-submit treeAnalyze_3p2GeV_nSigPr_high.xml >& submitOut_3p2GeV_nSigPr_high.txt &

nohup star-submit treeAnalyze_3p2GeV_m2Pi_low.xml >& submitOut_3p2GeV_m2Pi_low.txt &
nohup star-submit treeAnalyze_3p2GeV_m2Pi_high.xml >& submitOut_3p2GeV_m2Pi_high.txt &

nohup star-submit treeAnalyze_3p2GeV_m2Ka_low.xml >& submitOut_3p2GeV_m2Ka_low.txt &
nohup star-submit treeAnalyze_3p2GeV_m2Ka_high.xml >& submitOut_3p2GeV_m2Ka_high.txt &

nohup star-submit treeAnalyze_3p2GeV_epd_low.xml >& submitOut_3p2GeV_epd_low.txt &
nohup star-submit treeAnalyze_3p2GeV_epd_high.xml >& submitOut_3p2GeV_epd_high.txt &
