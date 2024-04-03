// Preprocessor to avoid loading header file multiple times
#ifndef TreeMaker_def
#define TreeMaker_def

// Load C/C++ header files
#include <stdio.h>
#include <iostream>
#include <vector>

// Load ROOT header files
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

// Load STARLibrary header files
#include "StMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoBbcHit.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoETofHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoETofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"

// Configuration file reader
#include "../ConfigReader/ConfigReader.h"

// IEvent file headers
//#include "StRoot/IClasses/IEventPlane.h"
//#include "StRoot/IClasses/IEvent.h"


class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

const Double_t PI = TMath::Pi();
const Int_t Ncentralities =  16;

const Int_t N_track    = 2048;       // Max tracks in an event
//This number may be unnecessary

const Int_t CENT_BINS  = 16;
const Int_t FIRST_CENT = 16 - CENT_BINS;
const Double_t D_M0_PI = 0.139571;   //Rest masses
const Double_t D_M0_KA = 0.493677;
const Double_t D_M0_PR = 0.938272;
const Double_t D_M0_DE = 1.875613;   // Deuteron
const Double_t D_M0_TR = 2.808921;   // Triton


class TreeMaker : public StMaker
{
private:
  StPicoDstMaker* mPicoDstMaker;
  ConfigReader    configs;
  TString         JobIdName;
  Double_t        cutTest;
  TFile*          outputFile;
  Double_t        Y_MID;
    
  TTree*              FxtTree;
  Int_t               tree_runId;
  Int_t               tree_eventId;
  //Int_t               tree_refmult;
  Float_t             tree_bField;
  Float_t             tree_Vx;
  Float_t             tree_Vy;
  Float_t             tree_Vz;
  UShort_t            tree_centrality;
  UShort_t            tree_tracknumber;
  UShort_t            tree_nEPDhits;
  Short_t             tree_EPDid[744];
  Float_t             tree_EPDnMip[744];
//  IEvent	      tree_tracks;  
  
  //This information contained in IEvent (IEventPlane)
  Short_t             tree_Charge[N_track];
  Float_t             tree_Px[N_track];
  Float_t             tree_Py[N_track];
  Float_t             tree_Pz[N_track];
  Float_t             tree_DCA[N_track];
  Float_t             tree_nSigmaPi[N_track];
  Float_t             tree_nSigmaKa[N_track];
  Float_t             tree_nSigmaPr[N_track];
  Float_t             tree_Beta[N_track];
  Float_t             tree_EtofBeta[N_track];
  Float_t             tree_dEdx[N_track];
  Float_t             tree_dEdxError[N_track];
  Int_t               tree_nHits[N_track];
  Int_t               tree_nHitsFit[N_track];
  Int_t               tree_nHitsPoss[N_track];
  Int_t               tree_nHitsDedx[N_track];
  
  
  TH1D*                h_eventCheck;
  TH1D*                h_zvtx;
  TH1D*                h_trackmult;
  TH1D*                h_refmult;
  TH1D*                h_tofmult;
  TH2D*                h2_trans_vtx;
  TH2D*                h2_refmult_vs_trackmult;
  TH2D*                h2_tofmult_vs_trackmult;
  TH2D*                h2_tofmult_vs_refmult;

  TF1* pileupBoundHigh_3p2GeV;
  TF1* pileupBoundLow_3p2GeV;


public:
  TreeMaker(StPicoDstMaker* Maker, std::string configFileName, TString JobId, Int_t EventsNumber, Double_t inputParameter);
  virtual ~TreeMaker();
  Int_t    Init();
  Int_t    Make();
  Int_t    Finish();
  Bool_t   IsGoodRun(Int_t runNumber, Double_t sqrt_s_NN);
  bool     foundInArray(int number, int array[], const int size);
  Double_t getRapidity(Double_t px, Double_t py, Double_t pz, Double_t mass);
  ClassDef(TreeMaker,1) // Class title
};

// End of preprocessor
#endif
