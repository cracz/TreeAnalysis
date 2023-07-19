#ifndef SETUPATTRIBUTES_H
#define SETUPATTRIBUTES_H

#include <iostream>
#include "TROOT.h"
#include "TObject.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"
#include "TKey.h"


class SetupAttributes
{
 private:
  Int_t runIteration = 0;
  Bool_t resolutionsFound = false;
  Bool_t tpcEfficienciesFound = false;
  Bool_t tofEfficienciesFound = false;
  Bool_t v1WeightsFound = false;

 public:
  TFile* correctionFile;
  TFile* resolutionFile;
  TFile* pikpEfficiencyFile;
  TFile* tofEfficiencyFile;
  TFile* v1WeightsInputFile;
  TH2D* h2_tracking_pp;
  TH2D* h2_tracking_pm;
  TH2D* h2_tracking_kp;
  TH2D* h2_tracking_km;
  TH2D* h2_tracking_pr;
  TH2D* h2_ratio_tof;
  TProfile2D* p2_TPCv1Weights;
  TProfile2D* p2_EPDv1Weights;
  Int_t getRunIteration();
  Bool_t resolutionsWereFound();
  Bool_t tpcEfficienciesWereFound();
  Bool_t tofEfficienciesWereFound();
  Bool_t v1WeightsWereFound();
  void setCorrectionFileAndRunIteration(TString fileName);
  void setResolutionFile(TString fileName);
  void setTPCEfficiencyFile(TString fileName);
  void setTOFEfficiencyFile(TString fileName);
  void setv1WeightsFile(TString fileName);
};// End namespace Setup


#endif
