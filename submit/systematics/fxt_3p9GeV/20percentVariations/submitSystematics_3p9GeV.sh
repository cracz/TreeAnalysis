#!/bin/bash

nohup star-submit ../../../treeAnalyze_3p9GeV.xml >& submitOut_3p9GeV_Normal.txt &

nohup star-submit treeAnalyze_3p9GeV_dca_low.xml >& submitOut_3p9GeV_dca_low.txt &
nohup star-submit treeAnalyze_3p9GeV_dca_high.xml >& submitOut_3p9GeV_dca_high.txt &

nohup star-submit treeAnalyze_3p9GeV_nhits_low.xml >& submitOut_3p9GeV_nhits_low.txt &
nohup star-submit treeAnalyze_3p9GeV_nhits_high.xml >& submitOut_3p9GeV_nhits_high.txt &

nohup star-submit treeAnalyze_3p9GeV_nhitsdedx_high.xml >& submitOut_3p9GeV_nhitsdedx_high.txt &

nohup star-submit treeAnalyze_3p9GeV_nhitsratio_low.xml >& submitOut_3p9GeV_nhitsratio_low.txt &
nohup star-submit treeAnalyze_3p9GeV_nhitsratio_high.xml >& submitOut_3p9GeV_nhitsratio_high.txt &

nohup star-submit treeAnalyze_3p9GeV_nSigPi_low.xml >& submitOut_3p9GeV_nSigPi_low.txt &
nohup star-submit treeAnalyze_3p9GeV_nSigPi_high.xml >& submitOut_3p9GeV_nSigPi_high.txt &

nohup star-submit treeAnalyze_3p9GeV_nSigKa_low.xml >& submitOut_3p9GeV_nSigKa_low.txt &
nohup star-submit treeAnalyze_3p9GeV_nSigKa_high.xml >& submitOut_3p9GeV_nSigKa_high.txt &

nohup star-submit treeAnalyze_3p9GeV_nSigPr_low.xml >& submitOut_3p9GeV_nSigPr_low.txt &
nohup star-submit treeAnalyze_3p9GeV_nSigPr_high.xml >& submitOut_3p9GeV_nSigPr_high.txt &

nohup star-submit treeAnalyze_3p9GeV_zDe_low.xml >& submitOut_3p9GeV_zDe_low.txt &
nohup star-submit treeAnalyze_3p9GeV_zDe_high.xml >& submitOut_3p9GeV_zDe_high.txt &

nohup star-submit treeAnalyze_3p9GeV_zTr_low.xml >& submitOut_3p9GeV_zTr_low.txt &
nohup star-submit treeAnalyze_3p9GeV_zTr_high.xml >& submitOut_3p9GeV_zTr_high.txt &

nohup star-submit treeAnalyze_3p9GeV_m2Pi_low.xml >& submitOut_3p9GeV_m2Pi_low.txt &
nohup star-submit treeAnalyze_3p9GeV_m2Pi_high.xml >& submitOut_3p9GeV_m2Pi_high.txt &

nohup star-submit treeAnalyze_3p9GeV_m2Ka_low.xml >& submitOut_3p9GeV_m2Ka_low.txt &
nohup star-submit treeAnalyze_3p9GeV_m2Ka_high.xml >& submitOut_3p9GeV_m2Ka_high.txt &

nohup star-submit treeAnalyze_3p9GeV_m2De_low.xml >& submitOut_3p9GeV_m2De_low.txt &
nohup star-submit treeAnalyze_3p9GeV_m2De_high.xml >& submitOut_3p9GeV_m2De_high.txt &

nohup star-submit treeAnalyze_3p9GeV_m2Tr_low.xml >& submitOut_3p9GeV_m2Tr_low.txt &
nohup star-submit treeAnalyze_3p9GeV_m2Tr_high.xml >& submitOut_3p9GeV_m2Tr_high.txt &

#nohup star-submit treeAnalyze_3p9GeV_epd_low.xml >& submitOut_3p9GeV_epd_low.txt &
#nohup star-submit treeAnalyze_3p9GeV_epd_high.xml >& submitOut_3p9GeV_epd_high.txt &
