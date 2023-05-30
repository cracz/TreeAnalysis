#ifndef SETUPATTRIBUTES_H
#define SETUPATTRIBUTES_H

#include <iostream>
#include <string>
#include "TROOT.h"
#include "TObject.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TString.h"
#include "TKey.h"


class SetupAttributes
{
 private:

 public:
  Int_t getRunIteration(TFile* correctionFile);
  Bool_t setResolutionFile(TFile*& filePointerAddress, TString fileName);
};// End namespace Setup


#endif
