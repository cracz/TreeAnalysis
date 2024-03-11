#!/bin/bash

#nohup star-submit ../../../treeAnalyze.xml >& submitOut_Normal.txt &

nohup star-submit treeAnalyze_rvtx_low.xml >& submitOut_rvtx_low.txt &
nohup star-submit treeAnalyze_rvtx_high.xml >& submitOut_rvtx_high.txt &

nohup star-submit treeAnalyze_zvtx_low.xml >& submitOut_zvtx_low.txt &
nohup star-submit treeAnalyze_zvtx_high.xml >& submitOut_zvtx_high.txt &

nohup star-submit treeAnalyze_dca_low.xml >& submitOut_dca_low.txt &
nohup star-submit treeAnalyze_dca_high.xml >& submitOut_dca_high.txt &

nohup star-submit treeAnalyze_nhits_low.xml >& submitOut_nhits_low.txt &
nohup star-submit treeAnalyze_nhits_high.xml >& submitOut_nhits_high.txt &

nohup star-submit treeAnalyze_nhitsdedx_high.xml >& submitOut_nhitsdedx_high.txt &

nohup star-submit treeAnalyze_nhitsratio_low.xml >& submitOut_nhitsratio_low.txt &
nohup star-submit treeAnalyze_nhitsratio_high.xml >& submitOut_nhitsratio_high.txt &

nohup star-submit treeAnalyze_nSigPi_low.xml >& submitOut_nSigPi_low.txt &
nohup star-submit treeAnalyze_nSigPi_high.xml >& submitOut_nSigPi_high.txt &

nohup star-submit treeAnalyze_nSigKa_low.xml >& submitOut_nSigKa_low.txt &
nohup star-submit treeAnalyze_nSigKa_high.xml >& submitOut_nSigKa_high.txt &

nohup star-submit treeAnalyze_nSigPr_low.xml >& submitOut_nSigPr_low.txt &
nohup star-submit treeAnalyze_nSigPr_high.xml >& submitOut_nSigPr_high.txt &

nohup star-submit treeAnalyze_m2Pi_low.xml >& submitOut_m2Pi_low.txt &
nohup star-submit treeAnalyze_m2Pi_high.xml >& submitOut_m2Pi_high.txt &

nohup star-submit treeAnalyze_m2Ka_low.xml >& submitOut_m2Ka_low.txt &
nohup star-submit treeAnalyze_m2Ka_high.xml >& submitOut_m2Ka_high.txt &

nohup star-submit treeAnalyze_epd_low.xml >& submitOut_epd_low.txt &
nohup star-submit treeAnalyze_epd_high.xml >& submitOut_epd_high.txt &
