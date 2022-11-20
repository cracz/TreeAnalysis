//////
// THIS IS A VERSION OF FlowAnalyzer.cxx THAT ANALYZES TREES RATHER THAN PICODSTS.
//
// This program reads trees made from picoDsts that have already passed 
// through event cuts. It performs an analysis of anisotropic
// flow coefficients. Many of the controls are, or should be, located
// within the config file that is supplied to the ConfigReader.
// The config files are used to maximize the modularity of the program
// and minimize hardcoding of settings/parameters. Event/track QA cuts, 
// bad runs, and centrality definitions should have been handled when 
// the trees were produced.
//
// Author: Cameron Racz
// Date: 2021
//////


//////
// Do not use vector::erase() to change any vectors. 
// This will invalidate any for() loops iterating over the vectors
// and make things much more complicated. For bad events after 
// creation of the "Event" vector, use the badEvent flag.
//////


// C++ headers
#include <iostream>
#include <vector>
#include <sys/resource.h>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TChain.h"
#include "TF1.h"
#include "TSystem.h"
#include "TKey.h"
#include "TStopwatch.h"

// EPD Util headers
#include "StRoot/StEpdUtil/StEpdGeom.h"

// Bichsel header
#include "StRoot/StBichsel/Bichsel.h"

// Configuration file reader
#include "StRoot/ConfigReader/ConfigReader.h"

// My Util Header
#include "FlowUtils.h"

// Bichsel Function
Double_t bichselZ(Double_t *x,Double_t *par) 
{
  Double_t pove   = x[0];
  Double_t poverm = pove/par[0];
  return TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(poverm),par[1]));
}

Double_t bichselZcharge2(Double_t *x,Double_t *par) 
{
  Double_t pove   = 2.0 * x[0];
  Double_t poverm = pove/par[0];
  return TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(poverm),par[1]));
}


//=========================================================
//          SOME CONTROLS
//=========================================================
const Int_t CENT_BINS  = 16;             // Number of centrality bins to show (max 16)  LEAVE AT 16 FOR NOW, BEST FOR RESOLUTION STUFF
const Int_t FIRST_CENT = 16 - CENT_BINS;            // Starting point for centrality dependent plots

const Double_t D_M0_PI = 0.139571;   //Rest masses
const Double_t D_M0_KA = 0.493677;
const Double_t D_M0_PR = 0.938272;
const Double_t D_M0_DE = 1.875613;   // Deuteron
const Double_t D_M0_TR = 2.808921;   // Triton
const Double_t D_M0_HE3 = 2.809414;   // Helium-3
const Double_t D_M0_AL = 3.727379;   // Alpha

const Double_t PI = TMath::Pi();

Int_t RUN_ITERATION = 0;
// 0 = No correction info yet; save raw (Xn,Yn) distributions
// 1 = Correction file found, but only <Xn> and <Yn> for re-centering.
//     Also save <sin> <cos> at this step for shifting in the next step.
// 2 = Correction file found, and <sin> <cos> values found so that shifting can be performed.
//=========================================================
//          
//=========================================================

int main(int argc, char *argv[])
{
  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();

  std::cout << "Initializing..." << std::endl;

  TString inFile = argv[1];
  TString jobID  = argv[2];
  std::string configFileName = argv[3];
  TString correctionFileName = argv[4];
  TString resolutionFileName = argv[5];

  if (gSystem->AccessPathName(inFile)) { std::cout << "Error reading input file!" << std::endl; return 1;}

  //=========================================================
  //          Set up various files
  //=========================================================
  ConfigReader configs;
  configs.read(configFileName);
  if (configs.errorFound()) { std::cout << "There was an error reading the configurations! Aborting analysis!" << std::endl; return 1; }
  configs.printAll();

  const Double_t ORDER_N = configs.order_n;   // Order of anisotropic flow (v_n)
  const Double_t ORDER_M = configs.order_m;   // Order of event plane angle (psi_m)
  const Double_t Y_MID   = configs.y_mid;     // Mid rapidity for the current energy
  const Double_t Y_BEAM  = configs.y_beam;    // Beam rapidity for the current energy
  TString ORDER_N_STR;
  TString ORDER_M_STR;
  ORDER_N_STR.Form("%.0f", ORDER_N);
  ORDER_M_STR.Form("%.0f", ORDER_M);
  const Double_t PSI_BOUNDS = TMath::Pi()/ORDER_M + 1;  // Boundaries for many histograms
  const Double_t Q_BOUNDS = 100;
  const Bool_t ODD_PLANE = ((int)ORDER_M % 2 == 1) ? true : false;

  //=== INITIALIZE TTREE
  Int_t N_track = 0;  // Max number of tracks in an event. Depends on energy and centrality definition!
  if      (configs.sqrt_s_NN == 3.0)  { N_track = 195;  }
  else if (configs.sqrt_s_NN == 3.22) { N_track = 287;  }
  else if (configs.sqrt_s_NN == 3.9)  { N_track = 287;  }  // UPDATE THIS WHEN CENTRALITY IS OFFICIAL!
  else if (configs.sqrt_s_NN == 4.49) { N_track = 195;  }  // UPDATE THIS WHEN CENTRALITY IS OFFICIAL!
  else if (configs.sqrt_s_NN == 7.2)  { N_track = 240;  }
  //else if (configs.sqrt_s_NN == 14.5) { N_track = 2048; } // UPDATE THIS WHEN THE CENTRALITY IS OFFICIAL
  else if (configs.sqrt_s_NN == 19.6) { N_track = 2048; } // UPDATE THIS WHEN THE CENTRALITY IS OFFICIAL
  else
    {
      std::cout << "Unknown energy! N_track is not set!" << std::endl;
      return 1;
    }

  Int_t i_runID;
  Int_t i_eventID;
  Float_t f_bField;
  Float_t f_xvtx;
  Float_t f_yvtx;
  Float_t f_zvtx;
  UShort_t i_centrality;
  UShort_t i_trackNumber;
  Short_t charge[N_track];
  Float_t Px[N_track];
  Float_t Py[N_track];
  Float_t Pz[N_track];
  Float_t DCA[N_track];
  Float_t nSigmaPi[N_track];
  Float_t nSigmaKa[N_track];
  Float_t nSigmaPr[N_track];
  Float_t tofBeta[N_track];
  Int_t nHits[N_track];
  Int_t nHitsFit[N_track];
  Int_t nHitsPoss[N_track];
  Int_t nHitsDedx[N_track];
  Float_t dEdx[N_track];
  UShort_t i_nEPDhits;
  Short_t EPDids[744];
  Float_t EPDnMip[744];


  TFile *inputFile = TFile::Open(inFile);
  if (!inputFile) { std::cout << "Input file could not be opened properly!" << std::endl; return 1; }

  TTree *tree = (TTree*)inputFile->Get("Autree");  
  tree->SetBranchAddress("runId", &i_runID);
  tree->SetBranchAddress("eventId",&i_eventID);
  tree->SetBranchAddress("bField",&f_bField);
  tree->SetBranchAddress("Vx",&f_xvtx);
  tree->SetBranchAddress("Vy",&f_yvtx);
  tree->SetBranchAddress("Vz",&f_zvtx);
  tree->SetBranchAddress("centrality",&i_centrality);
  tree->SetBranchAddress("tracknumber",&i_trackNumber);
  tree->SetBranchAddress("Charge",&charge);
  tree->SetBranchAddress("Px",&Px);
  tree->SetBranchAddress("Py",&Py);
  tree->SetBranchAddress("Pz",&Pz);
  tree->SetBranchAddress("DCA",&DCA);
  tree->SetBranchAddress("nSigmaPi",&nSigmaPi);
  tree->SetBranchAddress("nSigmaKa",&nSigmaKa);
  tree->SetBranchAddress("nSigmaPr",&nSigmaPr);
  tree->SetBranchAddress("tofBeta",&tofBeta);
  tree->SetBranchAddress("dEdx",&dEdx);
  tree->SetBranchAddress("nHits",&nHits);
  tree->SetBranchAddress("nHitsFit",&nHitsFit);
  tree->SetBranchAddress("nHitsPoss",&nHitsPoss);
  tree->SetBranchAddress("nHitsDedx",&nHitsDedx);
  tree->SetBranchAddress("nEPDhits",&i_nEPDhits);
  tree->SetBranchAddress("EPDid",&EPDids);
  tree->SetBranchAddress("EPDnMip",&EPDnMip);
  //=== END TTREE SETUP


  // INPUT FILE FOR CORRECTION INFORMATION
  TFile *correctionInputFile = TFile::Open(correctionFileName, "READ");
  if (!correctionInputFile)
    {
      RUN_ITERATION = 0;
      std::cout << "No correction file found." << std::endl
		<< "Re-centering and shifting will not be performed." << std::endl;
    }
  else
    {
      TKey *key;
      TIter next(correctionInputFile->GetListOfKeys());
      TProfile *profile;

      while( (key = (TKey*)next()) )
	{
	  TClass *cl = gROOT->GetClass(key->GetClassName());
	  
	  if (cl->InheritsFrom("TProfile"))
	    {
	      profile = (TProfile*)key->ReadObj();
	      if (profile->GetEntries() == 0)
		{
		  std::cout << "TProfiles are empty!" << std::endl
			    << "Re-centering will be performed and TProfiles will be filled." << std::endl;
		  RUN_ITERATION = 1;
		  break;
		}
	      else if (profile->GetEntries() != 0)
		{
		  std::cout << "Non-empty TProfiles found!" << std::endl
			    << "Re-centering and event plane shifting will be performed." << std::endl;
		  RUN_ITERATION = 2;
		  break;
		}
	    }
	}
    }
  ////

  // INPUT FILE FOR EVENT PLANE RESOLUTION INFORMATION
  Bool_t resolutionsFound = false;
  TFile *resolutionInputFile;
  if (RUN_ITERATION == 2) 
    { 
      resolutionInputFile = TFile::Open(resolutionFileName, "READ"); 
      if (!resolutionInputFile) { std::cout << "No resolution file was found!" << std::endl; }
      else 
	{ 
	  resolutionsFound = true;
	  std::cout << "Resolution file found!" << std::endl; 
	}
    }
  ////

  // INPUT FILE FOR TPC EFFICIENCY CORRECTIONS
  TString pdtEfficiencyFileName = "pdt_efficiency.root";
  TString pikEfficiencyFileName = "pik_efficiency.root";
  TFile *pdtEfficiencyFile;
  TFile *pikEfficiencyFile;
  Bool_t efficienciesFound = false;
  TH2D *h2_tracking_pp;
  TH2D *h2_tracking_pm;
  TH2D *h2_tracking_kp;
  TH2D *h2_tracking_km;
  TH2D *h2_tracking_pr;
  TH2D *h2_tracking_de;
  TH2D *h2_tracking_tr;
  if (RUN_ITERATION == 2)
    {
      pdtEfficiencyFile = TFile::Open(pdtEfficiencyFileName, "READ");
      pikEfficiencyFile = TFile::Open(pikEfficiencyFileName, "READ");
      if (!pdtEfficiencyFile || !pikEfficiencyFile) { std::cout << "One or both efficiency files missing! All efficiencies will default to 1!" << std::endl; }
      else 
	{ 
	  efficienciesFound = true;
	  std::cout << "TPC pik and pdt efficiency files were found!" << std::endl; 

	  h2_tracking_pp = (TH2D*)pikEfficiencyFile->Get("h2_ratio_pp");
	  h2_tracking_pm = (TH2D*)pikEfficiencyFile->Get("h2_ratio_pm");
	  h2_tracking_kp = (TH2D*)pikEfficiencyFile->Get("h2_ratio_kp");
	  h2_tracking_km = (TH2D*)pikEfficiencyFile->Get("h2_ratio_km");

	  h2_tracking_pr = (TH2D*)pdtEfficiencyFile->Get("tracking_p");
	  h2_tracking_de = (TH2D*)pdtEfficiencyFile->Get("tracking_d");
	  h2_tracking_tr = (TH2D*)pdtEfficiencyFile->Get("tracking_t");
	}
    }
  ////

  // OUTPUT FILE FOR CORRECTION INFORMATION
  TString correctionOutputName = "correctionInfo_OUTPUT_"+jobID+".root";
  TFile *correctionOutputFile;
  if (RUN_ITERATION == 0 || RUN_ITERATION == 1) { correctionOutputFile = new TFile(correctionOutputName, "RECREATE"); }
  ////

  // MAIN OUTPUT FILE
  TString outFile = jobID+".picoDst.result.root";
  TFile *outputFile = new TFile(outFile,"RECREATE");
  outputFile->cd();
  ////
  //=========================================================
  //          END file setup
  //=========================================================


  //=========================================================
  //          Bichsel Function Setup
  //=========================================================
  Double_t log2dx = 1.0;
  Double_t xStart = 0.01;
  Double_t xStop  = 3.0;
  Int_t npx = 10000;
  //                      Mass  log2(dx)
  Double_t params[2] = {  1.0,   log2dx  };

  params[0] = D_M0_DE;
  TF1 *bichselZ_de = new TF1(Form("BichselZ_de_log2dx_%i",(int)log2dx),bichselZ,xStart,xStop,2);
  if (!bichselZ_de) { std::cout << "De function error" << std::endl; return 1; }
  bichselZ_de->SetParameters(params); 
  bichselZ_de->SetNpx(npx);

  params[0] = D_M0_TR;
  TF1 *bichselZ_tr = new TF1(Form("BichselZ_tr_log2dx_%i",(int)log2dx),bichselZ,xStart,xStop,2);
  if (!bichselZ_tr) { std::cout << "Tr function error" << std::endl; return 1; }
  bichselZ_tr->SetParameters(params); 
  bichselZ_tr->SetNpx(npx);

  params[0] = D_M0_HE3;
  TF1 *bichselZ_he3 = new TF1(Form("BichselZ_he3_log2dx_%i",(int)log2dx),bichselZcharge2,xStart,xStop,2);
  if (!bichselZ_he3) { std::cout << "He3 function error" << std::endl; return 1; }
  bichselZ_he3->SetParameters(params); 
  bichselZ_he3->SetNpx(npx);

  params[0] = D_M0_AL;
  TF1 *bichselZ_al = new TF1(Form("BichselZ_al_log2dx_%i",(int)log2dx),bichselZcharge2,xStart,xStop,2);
  if (!bichselZ_al) { std::cout << "Al function error" << std::endl; return 1; }
  bichselZ_al->SetParameters(params); 
  bichselZ_al->SetNpx(npx);
  //=========================================================
  //          END Bichsel Function Setup
  //=========================================================


  // HISTOGRAMS

  // temporary variables when histogram bins/bounds depend on the energy
  int tempBins1 = 0;
  double tempLowBound1 = 0;
  double tempHighBound1 = 0;
  int tempBins2 = 0;
  double tempLowBound2 = 0;
  double tempHighBound2 = 0;


  TH1D *h_eventCheck = (TH1D*)inputFile->Get("h_eventCheck");
  h_eventCheck->SetStats(0);

  TH1D *h_trackCheck = new TH1D("h_trackCheck","Track number after each cut;;Tracks", 3, 0, 3);
  h_trackCheck->GetXaxis()->SetBinLabel(1,"Event Cuts Only");
  h_trackCheck->GetXaxis()->SetBinLabel(2,"QA Cuts");
  h_trackCheck->GetXaxis()->SetBinLabel(3,"PID Cuts");
  h_trackCheck->SetStats(0);


  TH1D *h_simulationCheck_pp = new TH1D ("h_simulationCheck_pp", "N_{#pi+} with no TPC efficiency", 3, 0, 3);
  TH1D *h_simulationCheck_pm = new TH1D ("h_simulationCheck_pm", "N_{#pi-} with no TPC efficiency", 3, 0, 3);
  TH1D *h_simulationCheck_kp = new TH1D ("h_simulationCheck_kp", "N_{K+} with no TPC efficiency", 3, 0, 3);
  TH1D *h_simulationCheck_km = new TH1D ("h_simulationCheck_km", "N_{K-} with no TPC efficiency", 3, 0, 3);
  TH1D *h_simulationCheck_pr = new TH1D ("h_simulationCheck_pr", "N_{pr} with no TPC efficiency", 3, 0, 3);
  TH1D *h_simulationCheck_de = new TH1D ("h_simulationCheck_de", "N_{de} with no TPC efficiency", 3, 0, 3);
  TH1D *h_simulationCheck_tr = new TH1D ("h_simulationCheck_tr", "N_{tr} with no TPC efficiency", 3, 0, 3);
  TH1D *h_simulationCheck_total = new TH1D ("h_simulationCheck_total", "Total N_{trk}", 3, 0, 3);

  TH1D *h_nhits       = new TH1D("h_nhits", "nHits;Number of hits;Tracks", 50, 0, 50);
  TH1D *h_nhits_dEdx  = new TH1D("h_nhits_dEdx","nHitsdEdx;Number of hits;Tracks", 50, 0, 50);
  TH1D *h_nhitsFit    = new TH1D("h_nhitsFit", "nHitsFit;Number of hits for fitting;Tracks", 50, 0, 50);
  TH1D *h_nhitsPoss   = new TH1D("h_nhitsPoss", "nHitsFit;Number of hits possible;Tracks", 50, 0, 50);
  TH1D *h_nhits_ratio = new TH1D("h_nhits_ratio","nhitsFit/nhitsPoss;nhitsFit/nhitsPoss;Tracks",200,0.0,2.0);    
  TH1D *h_DCA         = new TH1D("h_DCA","DCA (cm);DCA (cm);Tracks",100,0.0,10.0);
  TH1D *h_DCA_pr_alt  = new TH1D("h_DCA_pr_alt","Alt Proton DCA (cm);DCA (cm);Tracks",100,0.0,10.0);
  TH1D *h_DCA_de      = new TH1D("h_DCA_de","Deuteron DCA (cm);DCA (cm);Tracks",100,0.0,10.0);
  TH1D *h_DCA_tr      = new TH1D("h_DCA_tr","Triton DCA (cm);DCA (cm);Tracks",100,0.0,10.0);

  TH1D *h_primTracks = new TH1D("h_primTracks","Raw Number of Primary Tracks;Tracks;Events", 200, 0, 200);

  tempBins1      = (configs.fixed_target) ? 200 : 500;
  tempLowBound1  = (configs.fixed_target) ? 190.0 : -210.0;
  tempHighBound1 = 210.0;
  TH1D* h_zvtx = new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", tempBins1, tempLowBound1, tempHighBound1);

  TH1D *h_tileWeights = new TH1D("h_tileWeights", "EPD Tile Weights;nMIP Weights;Hits", 5, -1, 4);
  TH1D *h_centralities = new TH1D("h_centralities", "Centralities;Centrality ID;Events", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TH1D *h_centralities_final = new TH1D("h_centralities_final", "Final Good Centralities;Centrality ID;Events", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TH1D *h_trackmult = (TH1D*)inputFile->Get("h_trackmult");
  TH1D *h_refmult   = (TH1D*)inputFile->Get("h_refmult");
  TH1D *h_tofmult   = (TH1D*)inputFile->Get("h_tofmult");

  TH1D *h_pT = new TH1D("h_pT","p_{T};p_{T};Tracks",1000,0.0,10.0);
  TH1D *h_eta = new TH1D("h_eta","#eta;#eta;Tracks",600,-6.0,6.0);
  TH1D *h_phi = new TH1D("h_phi","#phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  TH2D *h2_dEdx_vs_qp = new TH2D("h2_dEdx_vs_qp", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 800, -2, 6, 1000, 0, 20);
  TH2D *h2_dEdx_vs_qp_charge2 = new TH2D("h2_dEdx_vs_qp_charge2", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 800, -2, 6, 1000, 0, 20);
  TH2D *h2_dEdx_vs_qp_half = new TH2D("h2_dEdx_vs_qp_half", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 600, 0, 6, 1000, 0, 20);
  TH2D *h2_beta_vs_qp = new TH2D("h2_beta_vs_qp","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_m2_vs_qp = new TH2D("h2_m2_vs_qp", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);

  TH1D *h_tofBeta = new TH1D("h_tofBeta", "TOF #beta;#beta;Tracks", 150, 0, 1.5);
  TH1D *h_m2 = new TH1D("h_m2", "m^{2};m^{2} (GeV^{2}/c^{4});Tracks", 1000, 0, 15);
  TH1D *h_m2_alpha_he3 = new TH1D("h_m2_alpha_he3", "m^{2};m^{2} (GeV^{2}/c^{4});Tracks", 1000, 0, 15);

  TH1D *h_mom_pp = new TH1D("h_mom_pp", "#pi^{+} Total Momentum;|p| (GeV);", 100, 0, 5);
  TH1D *h_mom_pm = new TH1D("h_mom_pm", "#pi^{-} Total Momentum;|p| (GeV);", 100, 0, 5);
  TH1D *h_mom_kp = new TH1D("h_mom_kp", "K^{+} Total Momentum;|p| (GeV);",   100, 0, 5);
  TH1D *h_mom_km = new TH1D("h_mom_km", "K^{-} Total Momentum;|p| (GeV);",   100, 0, 5);
  TH1D *h_mom_pr = new TH1D("h_mom_pr", "Proton Total Momentum;|p| (GeV);",  100, 0, 5);
  TH1D *h_mom_de = new TH1D("h_mom_de", "Deuteron Total Momentum;|p| (GeV);",100, 0, 5);
  TH1D *h_mom_tr = new TH1D("h_mom_tr", "Triton Total Momentum;|p| (GeV);",  100, 0, 5);

  TH1D *h_pT_pp = new TH1D("h_pT_pp","#pi^{+} p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  TH1D *h_pT_pm = new TH1D("h_pT_pm","#pi^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  TH1D *h_pT_kp = new TH1D("h_pT_kp","K^{+} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  TH1D *h_pT_km = new TH1D("h_pT_km","K^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  TH1D *h_pT_pr = new TH1D("h_pT_pr","Proton p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  TH1D *h_pT_de = new TH1D("h_pT_de","Deuteron p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  TH1D *h_pT_tr = new TH1D("h_pT_tr","Triton p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);

  TH2D *h2_pT_vs_cent_pr = new TH2D("h2_pT_vs_cent_pr", "Proton;Centrality;p_{T} (GeV/c)", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 1000, 0.0, 5.0);
  TH2D *h2_pT_vs_cent_de = new TH2D("h2_pT_vs_cent_de", "Deuteron;Centrality;p_{T} (GeV/c)", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 1000, 0.0, 5.0);
  TH2D *h2_pT_vs_cent_tr = new TH2D("h2_pT_vs_cent_tr", "Triton;Centrality;p_{T} (GeV/c)", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 1000, 0.0, 5.0);

  TH2D *h2_dndDeltaPhi_vs_cent_pr_alt = new TH2D("h2_dndDeltaPhi_pr_alt", "Alt Proton;Centrality;dN/d(#Delta#phi)", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 100, -1.0*PI, PI);
  TH2D *h2_dndDeltaPhi_vs_cent_de     = new TH2D("h2_dndDeltaPhi_de", "Deuteron;Centrality;dN/d(#Delta#phi)", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 100, -1.0*PI, PI);
  TH2D *h2_dndDeltaPhi_vs_cent_tr     = new TH2D("h2_dndDeltaPhi_tr", "Triton;Centrality;dN/d(#Delta#phi)", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 100, -1.0*PI, PI);

  TH1D *h_eta_pp = new TH1D("h_eta_pp","#pi^{+} #eta;#eta;Tracks",500,-5.0,5.0);
  TH1D *h_eta_pm = new TH1D("h_eta_pm","#pi^{-} #eta;#eta;Tracks",500,-5.0,5.0);
  TH1D *h_eta_kp = new TH1D("h_eta_kp","K^{+} #eta;#eta;Tracks",500,-5.0,5.0);
  TH1D *h_eta_km = new TH1D("h_eta_km","K^{-} #eta;#eta;Tracks",500,-5.0,5.0);
  TH1D *h_eta_pr = new TH1D("h_eta_pr","Proton #eta;#eta;Tracks",500,-5.0,5.0);
  TH1D *h_eta_de = new TH1D("h_eta_de","Deuteron #eta;#eta;Tracks",500,-5.0,5.0);
  TH1D *h_eta_tr = new TH1D("h_eta_tr","Triton #eta;#eta;Tracks",500,-5.0,5.0);

  TH1D *h_dndy_pp = new TH1D("h_dndy_pp", "#pi^{+} Raw Rapidity Spectrum;y;dN/dy", 80, -2, 2);
  TH1D *h_dndy_pm = new TH1D("h_dndy_pm", "#pi^{-} Raw Rapidity Spectrum;y;dN/dy", 80, -2, 2);
  TH1D *h_dndy_kp = new TH1D("h_dndy_kp", "K^{+} Raw Rapidity Spectrum;y;dN/dy",   80, -2, 2);
  TH1D *h_dndy_km = new TH1D("h_dndy_km", "K^{-} Raw Rapidity Spectrum;y;dN/dy",   80, -2, 2);
  TH1D *h_dndy_pr = new TH1D("h_dndy_pr", "Proton Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);
  TH1D *h_dndy_de = new TH1D("h_dndy_de", "Deuteron Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);
  TH1D *h_dndy_tr = new TH1D("h_dndy_tr", "Triton Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);
  
  TH1D *h_phi_pp = new TH1D("h_phi_pp","#pi^{+} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  TH1D *h_phi_pm = new TH1D("h_phi_pm","#pi^{-} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  TH1D *h_phi_kp = new TH1D("h_phi_kp","K^{+} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  TH1D *h_phi_km = new TH1D("h_phi_km","K^{-} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  TH1D *h_phi_pr = new TH1D("h_phi_pr","Proton #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  TH1D *h_phi_de = new TH1D("h_phi_de","Deuteron #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  TH1D *h_phi_tr = new TH1D("h_phi_tr","Triton #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);

  TH1D *h_dndm_pp = new TH1D("h_dndm_pp", "#pi^{+} Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}", 60, 0, 3);
  TH1D *h_dndm_pm = new TH1D("h_dndm_pm", "#pi^{-} Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}", 60, 0, 3);
  TH1D *h_dndm_kp = new TH1D("h_dndm_kp", "K^{+} Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}",   60, 0, 3);
  TH1D *h_dndm_km = new TH1D("h_dndm_km", "K^{-} Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}",   60, 0, 3);
  TH1D *h_dndm_pr = new TH1D("h_dndm_pr", "Proton Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}",  60, 0, 3);
  TH1D *h_dndm_de = new TH1D("h_dndm_de", "Deuteron Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}",60, 0, 3);
  TH1D *h_dndm_tr = new TH1D("h_dndm_tr", "Triton Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}",  60, 0, 3);
  
  TH1D *h_mult_pp = new TH1D("h_mult_pp","#pi^{+} track multiplicity;#pi^{+} Mult;Events",1001,-0.5,1000.5);
  TH1D *h_mult_pm = new TH1D("h_mult_pm","#pi^{-} track multiplicity;#pi^{-} Mult;Events",1001,-0.5,1000.5);
  TH1D *h_mult_kp = new TH1D("h_mult_kp","K^{#plus} track multiplicity;K^{+} Mult;Events",1001,-0.5,1000.5);
  TH1D *h_mult_km = new TH1D("h_mult_km","K^{-} track multiplicity;K^{-} Mult;Events",1001,-0.5,1000.5);
  TH1D *h_mult_pr = new TH1D("h_mult_pr","Proton track multiplicity;Proton Mult;Events",1001,-0.5,1000.5);
  TH1D *h_mult_de = new TH1D("h_mult_de","Deuteron track multiplicity;Deuteron Mult;Events",1001,-0.5,1000.5);
  TH1D *h_mult_tr = new TH1D("h_mult_tr","Triton track multiplicity;Triton Mult;Events",1001,-0.5,1000.5);
  
  TH2D *h2_dEdx_vs_qp_pp = new TH2D("h2_dEdx_vs_qp_pp", "#pi^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  TH2D *h2_dEdx_vs_qp_pm = new TH2D("h2_dEdx_vs_qp_pm", "#pi^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  TH2D *h2_dEdx_vs_qp_kp = new TH2D("h2_dEdx_vs_qp_kp", "K^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  TH2D *h2_dEdx_vs_qp_km = new TH2D("h2_dEdx_vs_qp_km", "K^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  TH2D *h2_dEdx_vs_qp_pr = new TH2D("h2_dEdx_vs_qp_pr", "Proton dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  TH2D *h2_dEdx_vs_qp_de = new TH2D("h2_dEdx_vs_qp_de", "Deuteron dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  TH2D *h2_dEdx_vs_qp_tr = new TH2D("h2_dEdx_vs_qp_tr", "Triton dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  
  TH2D *h2_beta_vs_qp_pp = new TH2D("h2_beta_vs_qp_pp","#pi^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_beta_vs_qp_pm = new TH2D("h2_beta_vs_qp_pm","#pi^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_beta_vs_qp_kp = new TH2D("h2_beta_vs_qp_kp","K^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_beta_vs_qp_km = new TH2D("h2_beta_vs_qp_km","K^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  //TH2D *h2_beta_vs_qp_pr = new TH2D("h2_beta_vs_qp_pr","Proton 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  //TH2D *h2_beta_vs_qp_de = new TH2D("h2_beta_vs_qp_de","Deuteron 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  //TH2D *h2_beta_vs_qp_tr = new TH2D("h2_beta_vs_qp_tr","Triton 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
 
  TH2D *h2_m2_vs_qp_pp = new TH2D("h2_m2_vs_qp_pp", "#pi^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  TH2D *h2_m2_vs_qp_pm = new TH2D("h2_m2_vs_qp_pm", "#pi^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  TH2D *h2_m2_vs_qp_kp = new TH2D("h2_m2_vs_qp_kp", "K^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
  TH2D *h2_m2_vs_qp_km = new TH2D("h2_m2_vs_qp_km", "K^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
  //TH2D *h2_m2_vs_qp_pr = new TH2D("h2_m2_vs_qp_pr", "Proton m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  //TH2D *h2_m2_vs_qp_de = new TH2D("h2_m2_vs_qp_de", "Deuteron m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  //TH2D *h2_m2_vs_qp_tr = new TH2D("h2_m2_vs_qp_tr", "Triton m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);

  if (configs.sqrt_s_NN == 3.0)
    {
      tempBins1 = 300;
      tempLowBound1 = -1.2;
      tempHighBound1 = 1.2;
      tempBins2 = 300;
      tempLowBound2  = 0.0;
      tempHighBound2 = 2.5;
    }
  else if (configs.sqrt_s_NN == 3.22 || configs.sqrt_s_NN == 3.9)
    {
      tempBins1 = 300;
      tempLowBound1 = -1.2;
      tempHighBound1 = 1.2;
      tempBins2 = 300;
      tempLowBound2  = 0.0;
      tempHighBound2 = 2.5;
    }
  else if (configs.sqrt_s_NN == 4.49)
    {
      tempBins1 = 300;
      tempLowBound1 = -1.0;
      tempHighBound1 = 1.7;
      tempBins2 = 300;
      tempLowBound2  = 0.0;
      tempHighBound2 = 2.5;
    }
  else if (configs.sqrt_s_NN == 7.2)
    {
      tempBins1 = 300;
      tempLowBound1 = -0.2;
      tempHighBound1 = 2.2;
      tempBins2 = 300;
      tempLowBound2  = 0.0;
      tempHighBound2 = 2.5;
    }
  else if (configs.sqrt_s_NN == 19.6)
    {
      tempBins1 = 400;
      tempLowBound1 = -2.0;
      tempHighBound1 = 2.0;
      tempBins2 = 500;
      tempLowBound2  = 0.0;
      tempHighBound2 = 5.0;
    }
  else
    {
      tempBins1 = 0;
      tempLowBound1 = 0.0;
      tempHighBound1 = 0.0;
      tempBins2 = 0;
      tempLowBound2  = 0.0;
      tempHighBound2 = 0.0;
    }
  TH2D *h2_pT_vs_yCM_pp = new TH2D("h2_pT_vs_yCM_pp", "#pi^{+};y-y_{mid};p_{T} (GeV/c)", tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_pT_vs_yCM_pm = new TH2D("h2_pT_vs_yCM_pm", "#pi^{-};y-y_{mid};p_{T} (GeV/c)", tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_pT_vs_yCM_kp = new TH2D("h2_pT_vs_yCM_kp", "K^{+};y-y_{mid};p_{T} (GeV/c)",   tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_pT_vs_yCM_km = new TH2D("h2_pT_vs_yCM_km", "K^{-};y-y_{mid};p_{T} (GeV/c)",   tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_pT_vs_yCM_pr = new TH2D("h2_pT_vs_yCM_pr", "Proton;y-y_{mid};p_{T} (GeV/c)",  tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_pT_vs_yCM_pr_alt = new TH2D("h2_pT_vs_yCM_pr_alt", "Proton;y-y_{mid};p_{T} (GeV/c)",  tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_pT_vs_yCM_de = new TH2D("h2_pT_vs_yCM_de", "Deuteron;y-y_{mid};p_{T} (GeV/c)",tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_pT_vs_yCM_tr = new TH2D("h2_pT_vs_yCM_tr", "Triton;y-y_{mid};p_{T} (GeV/c)",  tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);

  TH2D *h2_pT_vs_yCM_pp_noEff = new TH2D("h2_pT_vs_yCM_pp_noEff", "#pi^{+} with No TPC Efficiency;y-y_{mid};p_{T} (GeV/c)",  
					 tempBins1, tempLowBound1, tempHighBound1, 500, 0.0, 5.0);
  TH2D *h2_pT_vs_yCM_pm_noEff = new TH2D("h2_pT_vs_yCM_pm_noEff", "#pi^{-} with No TPC Efficiency;y-y_{mid};p_{T} (GeV/c)",  
					 tempBins1, tempLowBound1, tempHighBound1, 500, 0.0, 5.0);
  TH2D *h2_pT_vs_yCM_kp_noEff = new TH2D("h2_pT_vs_yCM_kp_noEff", "K^{+} with No TPC Efficiency;y-y_{mid};p_{T} (GeV/c)",  
					 tempBins1, tempLowBound1, tempHighBound1, 500, 0.0, 5.0);
  TH2D *h2_pT_vs_yCM_km_noEff = new TH2D("h2_pT_vs_yCM_km_noEff", "K^{-} with No TPC Efficiency;y-y_{mid};p_{T} (GeV/c)",  
					 tempBins1, tempLowBound1, tempHighBound1, 500, 0.0, 5.0);
  TH2D *h2_pT_vs_yCM_pr_noEff = new TH2D("h2_pT_vs_yCM_pr_noEff", "Protons with No TPC Efficiency;y-y_{mid};p_{T} (GeV/c)",  
					 tempBins1, tempLowBound1, tempHighBound1, 500, 0.0, 5.0);
  TH2D *h2_pT_vs_yCM_de_noEff = new TH2D("h2_pT_vs_yCM_de_noEff", "Deuterons with No TPC Efficiency;y-y_{mid};p_{T} (GeV/c)",  
					 tempBins1, tempLowBound1, tempHighBound1, 500, 0.0, 5.0);
  TH2D *h2_pT_vs_yCM_tr_noEff = new TH2D("h2_pT_vs_yCM_tr_noEff", "Tritons with No TPC Efficiency;y-y_{mid};p_{T} (GeV/c)",  
					 tempBins1, tempLowBound1, tempHighBound1, 500, 0.0, 5.0);

  TH1D *h_psiTpc_RAW  = new TH1D("h_psiTpc_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", TPC);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcA_RAW = new TH1D("h_psiTpcA_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", TPC A);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcB_RAW = new TH1D("h_psiTpcB_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", TPC B);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpd_RAW  = new TH1D("h_psiEpd_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", EPD);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdA_RAW = new TH1D("h_psiEpdA_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", EPD A);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdB_RAW = new TH1D("h_psiEpdB_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", EPD B);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  TProfile *p_meanpT_vs_yCM_pp = new TProfile("p_meanpT_vs_yCM_pp","#pi^{+} <p_{T}>;y-y_{mid};<p_{T}>", 20, -1.0, 1.0);
  TProfile *p_meanpT_vs_yCM_pm = new TProfile("p_meanpT_vs_yCM_pm","#pi^{-} <p_{T}>;y-y_{mid};<p_{T}>", 20, -1.0, 1.0);
  TProfile *p_meanpT_vs_yCM_kp = new TProfile("p_meanpT_vs_yCM_kp","K^{+} <p_{T}>;y-y_{mid};<p_{T}>", 20, -1.0, 1.0);
  TProfile *p_meanpT_vs_yCM_km = new TProfile("p_meanpT_vs_yCM_km","K^{-} <p_{T}>;y-y_{mid};<p_{T}>", 20, -1.0, 1.0);
  TProfile *p_meanpT_vs_yCM_pr = new TProfile("p_meanpT_vs_yCM_pr","Proton <p_{T}>;y-y_{mid};<p_{T}>", 20, -1.0, 1.0);
  TProfile *p_meanpT_vs_yCM_pr_alt = new TProfile("p_meanpT_vs_yCM_pr_alt","Proton <p_{T}>;y-y_{mid};<p_{T}>", 20, -1.0, 1.0);
  TProfile *p_meanpT_vs_yCM_de = new TProfile("p_meanpT_vs_yCM_de","Deuteron <p_{T}>;y-y_{mid};<p_{T}>", 20, -1.0, 1.0);
  TProfile *p_meanpT_vs_yCM_tr = new TProfile("p_meanpT_vs_yCM_tr","Triton <p_{T}>;y-y_{mid};<p_{T}>", 20, -1.0, 1.0);


  TProfile *p_vn_pp = new TProfile("p_vn_pp", "#pi^{+} v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pm = new TProfile("p_vn_pm", "#pi^{-} v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_kp = new TProfile("p_vn_kp", "K^{+} v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_km = new TProfile("p_vn_km", "K^{-} v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pr = new TProfile("p_vn_pr", "Proton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pr_alt = new TProfile("p_vn_pr_alt", "Proton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pr_pTlt1 = new TProfile("p_vn_pr_pTlt1", "Proton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
					 CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pr_pTgt1 = new TProfile("p_vn_pr_pTgt1", "Proton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
					 CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_de = new TProfile("p_vn_de", "Deuteron v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_tr = new TProfile("p_vn_tr", "Triton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);

  TProfile *p_vn_pr_alt_y0to0p6 = new TProfile("p_vn_pr_alt_y0to0p6", 
					       "Proton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
					       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_de_y0to0p6 = new TProfile("p_vn_de_y0to0p6", 
					   "Deuteron v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}",
					   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_tr_y0to0p6 = new TProfile("p_vn_tr_y0to0p6", 
					   "Triton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
					   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);

  TProfile *p_vn_pr_alt_y0p6to1p0 = new TProfile("p_vn_pr_alt_y0p6to1p0", 
						 "Proton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
						 CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_de_y0p6to1p0 = new TProfile("p_vn_de_y0p6to1p0", 
					     "Deuteron v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
					     CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_tr_y0p6to1p0 = new TProfile("p_vn_tr_y0p6to1p0", 
					     "Triton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
					     CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);

  // vn profiles at "extended" rapidity range 0.5 < y_CM < 1.0
  TProfile *p_vn_pp_ext = new TProfile("p_vn_pp_ext", "#pi^{+} v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pm_ext = new TProfile("p_vn_pm_ext", "#pi^{-} v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_kp_ext = new TProfile("p_vn_kp_ext", "K^{+} v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_km_ext = new TProfile("p_vn_km_ext", "K^{-} v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pr_ext = new TProfile("p_vn_pr_ext", "Proton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);

  // vn profiles at the "forward" raidity range y_CM < 0
  TProfile *p_vn_pr_for = new TProfile("p_vn_pr_for", "Proton v_{"+ORDER_N_STR+"};Centrality;v_{"+ORDER_N_STR+"}{#psi_{"+ORDER_M_STR+"}}/R_{"+ORDER_N_STR+ORDER_M_STR+"}", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);


  // "Observed" flow values (uncorrected by event plane resolution)
  TProfile *p_vn_pp_obs = new TProfile("p_vn_pp_obs", "#pi^{+} v_{"+ORDER_N_STR+"}^{obs};Centrality;v_{"+ORDER_N_STR+"}^{obs}{#psi_{"+ORDER_M_STR+"}}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pm_obs = new TProfile("p_vn_pm_obs", "#pi^{-} v_{"+ORDER_N_STR+"}^{obs};Centrality;v_{"+ORDER_N_STR+"}^{obs}{#psi_{"+ORDER_M_STR+"}}", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_kp_obs = new TProfile("p_vn_kp_obs", "K^{+} v_{"+ORDER_N_STR+"}^{obs};Centrality;v_{"+ORDER_N_STR+"}^{obs}{#psi_{"+ORDER_M_STR+"}}", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_km_obs = new TProfile("p_vn_km_obs", "K^{-} v_{"+ORDER_N_STR+"}^{obs};Centrality;v_{"+ORDER_N_STR+"}^{obs}{#psi_{"+ORDER_M_STR+"}}", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pr_obs = new TProfile("p_vn_pr_obs", "Proton v_{"+ORDER_N_STR+"}^{obs};Centrality;v_{"+ORDER_N_STR+"}^{obs}{#psi_{"+ORDER_M_STR+"}}", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pr_alt_obs = new TProfile("p_vn_pr_alt_obs", "Proton v_{"+ORDER_N_STR+"}^{obs};Centrality;v_{"+ORDER_N_STR+"}^{obs}{#psi_{"+ORDER_M_STR+"}}", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_de_obs = new TProfile("p_vn_de_obs", "Deuteron v_{"+ORDER_N_STR+"}^{obs};Centrality;v_{"+ORDER_N_STR+"}^{obs}{#psi_{"+ORDER_M_STR+"}}", 
				   CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_tr_obs = new TProfile("p_vn_tr_obs", "Triton v_{"+ORDER_N_STR+"}^{obs};Centrality;v_{"+ORDER_N_STR+"}^{obs}{#psi_{"+ORDER_M_STR+"}}", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);


  // Differential Flow Profiles
  //TProfile2D *p2_v1_pT_eta_TPCB_pr = new TProfile2D("p2_v1_pT_eta_TPCB_pr", "TPC B Proton v_{1};#eta;p_{T}", 25, -1.0, 0.0, 25, 0, 2.5);
  //TProfile *p_v1_EPD_ring = new TProfile("p_v1_EPD_ring", "EPD v_{1} by Ring;Ring;v_{1}", 16, 0.5, 16.5);

  TProfile2D *p2_vn_yCM_cent_pp = new TProfile2D("p2_vn_yCM_cent_pp", "#pi^{+} v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);
  TProfile2D *p2_vn_yCM_cent_pm = new TProfile2D("p2_vn_yCM_cent_pm", "#pi^{-} v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);
  TProfile2D *p2_vn_yCM_cent_kp = new TProfile2D("p2_vn_yCM_cent_kp", "K^{+} v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);
  TProfile2D *p2_vn_yCM_cent_km = new TProfile2D("p2_vn_yCM_cent_km", "K^{-} v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);
  TProfile2D *p2_vn_yCM_cent_pr = new TProfile2D("p2_vn_yCM_cent_pr", "Proton v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);
  TProfile2D *p2_vn_yCM_cent_pr_alt = 
    new TProfile2D("p2_vn_yCM_cent_pr_alt", "Proton v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);
  TProfile2D *p2_vn_yCM_cent_pr_symmetry = 
    new TProfile2D("p2_vn_yCM_cent_pr_symmetry", "Proton v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);
  TProfile2D *p2_vn_yCM_cent_de = new TProfile2D("p2_vn_yCM_cent_de", "Deuteron v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);
  TProfile2D *p2_vn_yCM_cent_tr = new TProfile2D("p2_vn_yCM_cent_tr", "Triton v_{"+ORDER_N_STR+"};Centrality;y-y_{mid}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 20, -1, 1);

  TProfile2D *p2_vn_yOverYbeam_cent_pr_alt = 
    new TProfile2D("p2_vn_yOverYbeam_cent_pr_alt", "Proton v_{"+ORDER_N_STR+"};Centrality;y/|y_{beam}|", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 25, -0.5, 0.0);
  TProfile2D *p2_vn_yOverYbeam_cent_de = 
    new TProfile2D("p2_vn_yOverYbeam_cent_de", "Deuteron v_{"+ORDER_N_STR+"};Centrality;y/|y_{beam}|", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 25, -0.5, 0.0);
  TProfile2D *p2_vn_yOverYbeam_cent_tr = 
    new TProfile2D("p2_vn_yOverYbeam_cent_tr", "Triton v_{"+ORDER_N_STR+"};Centrality;y/|y_{beam}|", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 25, -0.5, 0.0);

  
  TProfile2D *p2_vn_pT_cent_pp = new TProfile2D("p2_vn_pT_cent_pp", "#pi^{+} v_{"+ORDER_N_STR+"};Centrality;p_{T}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_pT_cent_pm = new TProfile2D("p2_vn_pT_cent_pm", "#pi^{-} v_{"+ORDER_N_STR+"};Centrality;p_{T}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_pT_cent_kp = new TProfile2D("p2_vn_pT_cent_kp", "K^{+} v_{"+ORDER_N_STR+"};Centrality;p_{T}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_pT_cent_km = new TProfile2D("p2_vn_pT_cent_km", "K^{-} v_{"+ORDER_N_STR+"};Centrality;p_{T}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_pT_cent_pr = new TProfile2D("p2_vn_pT_cent_pr", "Proton v_{"+ORDER_N_STR+"};Centrality;p_{T}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_pT_cent_pr_alt = new TProfile2D("p2_vn_pT_cent_pr_alt", "Proton v_{"+ORDER_N_STR+"};Centrality;p_{T}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 15, 0, 2.5);
  TProfile2D *p2_vn_pT_cent_de = new TProfile2D("p2_vn_pT_cent_de", "Deuteron v_{"+ORDER_N_STR+"};Centrality;p_{T}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 15, 0, 2.5);
  TProfile2D *p2_vn_pT_cent_tr = new TProfile2D("p2_vn_pT_cent_tr", "Triton v_{"+ORDER_N_STR+"};Centrality;p_{T}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 15, 0, 2.5);

  // yCM stratified by pT
  TProfile2D *p2_vn_pT_vs_yCM_pp = new TProfile2D("p2_vn_pT_vs_yCM_pp", "#pi^{+} v_{3};y-y_{mid};p_{T} (GeV/c)", 20, -1.0, 1.0, 10, 0.0, 2.0);
  TProfile2D *p2_vn_pT_vs_yCM_pm = new TProfile2D("p2_vn_pT_vs_yCM_pm", "#pi^{-} v_{3};y-y_{mid};p_{T} (GeV/c)", 20, -1.0, 1.0, 10, 0.0, 2.0);
  TProfile2D *p2_vn_pT_vs_yCM_pr = new TProfile2D("p2_vn_pT_vs_yCM_pr", "Proton v_{3};y-y_{mid};p_{T} (GeV/c)",  20, -1.0, 1.0, 10, 0.0, 2.5);
  ////

  TProfile2D *p2_vn_KT_cent_pp = new TProfile2D("p2_vn_KT_cent_pp", "#pi^{+} v_{"+ORDER_N_STR+"};Centrality;m_{T}-m_{0}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_KT_cent_pm = new TProfile2D("p2_vn_KT_cent_pm", "#pi^{-} v_{"+ORDER_N_STR+"};Centrality;m_{T}-m_{0}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_KT_cent_kp = new TProfile2D("p2_vn_KT_cent_kp", "K^{+} v_{"+ORDER_N_STR+"};Centrality;m_{T}-m_{0}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_KT_cent_km = new TProfile2D("p2_vn_KT_cent_km", "K^{-} v_{"+ORDER_N_STR+"};Centrality;m_{T}-m_{0}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_KT_cent_pr = new TProfile2D("p2_vn_KT_cent_pr", "Proton v_{"+ORDER_N_STR+"};Centrality;m_{T}-m_{0}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_KT_cent_pr_alt = new TProfile2D("p2_vn_KT_cent_pr_alt", "Proton v_{"+ORDER_N_STR+"};Centrality;m_{T}-m_{0}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_KT_cent_de = new TProfile2D("p2_vn_KT_cent_de", "Deuteron v_{"+ORDER_N_STR+"};Centrality;m_{T}-m_{0}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);
  TProfile2D *p2_vn_KT_cent_tr = new TProfile2D("p2_vn_KT_cent_tr", "Triton v_{"+ORDER_N_STR+"};Centrality;m_{T}-m_{0}", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS, 10, 0, 2);

  // Profiles for resolution terms
  TProfile *p_TpcAB = new TProfile("p_TpcAB","TPC A-B Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,A}_{"+ORDER_M_STR+"}-#psi^{TPC,B}_{"+ORDER_M_STR+"}))>",
				   CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_TpcAEpdA = new TProfile("p_TpcAEpdA","TPC A EPD A Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,A}_{"+ORDER_M_STR+"}-#psi^{EPD,A}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdB = new TProfile("p_TpcAEpdB","TPC A EPD B Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,A}_{"+ORDER_M_STR+"}-#psi^{EPD,B}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_TpcBEpdA = new TProfile("p_TpcBEpdA","TPC B EPD A Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,B}_{"+ORDER_M_STR+"}-#psi^{EPD,A}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdB = new TProfile("p_TpcBEpdB","TPC B EPD B Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,B}_{"+ORDER_M_STR+"}-#psi^{EPD,B}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_EpdAEpdB = new TProfile("p_EpdAEpdB","EPD A EPD B Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{EPD,A}_{"+ORDER_M_STR+"}-#psi^{EPD,B}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  ////

  tempBins1      = (configs.fixed_target) ? 400 : 800;
  tempLowBound1  = (configs.fixed_target) ? -6.0 : -6.0;
  tempHighBound1 = (configs.fixed_target) ? -2.0 : 6.0;
  TProfile2D *p2_pp_vs_eta = new TProfile2D("p2_pp_vs_eta","<TnMIP> for Supersectors vs #eta;#eta;Supersector", tempBins1, tempLowBound1, tempHighBound1, 12, 0.5, 12.5);

  TH2D *h2_ring_vs_eta = new TH2D("h2_ring_vs_eta","EPD East Ring vs #eta;#eta;Ring", 500, -6.0, -1.0, 16, 0.5, 16.5);

  TH2D *h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  TH2D *h2_trans_vtx_cut = new TH2D("h2_trans_vtx_cut","Final Primary Vertices;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);

  TH2D *h2_refmult_vs_trackmult = (TH2D*)inputFile->Get("h2_refmult_vs_trackmult");
  TH2D *h2_tofmult_vs_trackmult = (TH2D*)inputFile->Get("h2_tofmult_vs_trackmult");
  TH2D *h2_tofmult_vs_refmult   = (TH2D*)inputFile->Get("h2_tofmult_vs_refmult");

  TH2D *h2_hits_vs_cent_EpdA = new TH2D("h2_nHits_vs_cent_EpdA", "EPD A;Centrality;Hits", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS, 50, 0, 50);
  TH2D *h2_hits_vs_cent_EpdB = new TH2D("h2_nHits_vs_cent_EpdB", "EPD B;Centrality;Hits", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS, 50, 0, 50);
  TH2D *h2_hits_vs_cent_TpcB = new TH2D("h2_nHits_vs_cent_TpcB", "TPC B;Centrality;Hits", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS, 50, 0, 50);

  TH2D *h2_dEdx_vs_qp_id_pr = new TH2D("h2_dEdx_vs_qp_id_pr", ";|p| (GeV/c);dE/dx (keV/cm)", 25, 0.0, 2.5, 500, 0.0, 20.0);
  TH2D *h2_dEdx_vs_qp_id_pr_alt = new TH2D("h2_dEdx_vs_qp_id_pr_alt", ";|p| (GeV/c);dE/dx (keV/cm)", 25, 0.0, 2.5, 500, 0.0, 20.0);
  TH2D *h2_dEdx_vs_qp_id_de = new TH2D("h2_dEdx_vs_qp_id_de", ";|p| (GeV/c);dE/dx (keV/cm)", 25, 0.0, 2.5, 500, 0.0, 20.0);
  TH2D *h2_dEdx_vs_qp_id_tr = new TH2D("h2_dEdx_vs_qp_id_tr", ";|p| (GeV/c);dE/dx (keV/cm)", 25, 0.0, 2.5, 500, 0.0, 20.0);

  TH2D *h2_nSigp_vs_mom = new TH2D("h2_nSigp_vs_mom", ";|p| (GeV/c);n#sigma_{p}", 40, 0.0, 4.0, 600, -3.0, 3.0);
  TH2D *h2_zd_vs_mom = new TH2D("h2_zd_vs_mom", ";|p| (GeV/c);z_{d}", 40, 0.0, 4.0, 140, -0.7, 0.7);
  TH2D *h2_zt_vs_mom = new TH2D("h2_zt_vs_mom", ";|p| (GeV/c);z_{t}", 40, 0.0, 4.0, 140, -0.7, 0.7);
  

  tempBins1      = (configs.fixed_target) ?  300 :  600;
  tempLowBound1  = (configs.fixed_target) ? -2.2 : -2.5;
  tempHighBound1 = (configs.fixed_target) ?  0.2 :  2.5;
  TH2D *h2_phi_vs_eta_TPC = new TH2D("h2_phi_vs_eta_TPC", "TPC;#eta;#phi", tempBins1, tempLowBound1, tempHighBound1, 300, -4, 4);

  tempBins1      = (configs.fixed_target) ?  300 :  600;
  tempLowBound1  = (configs.fixed_target) ? -6.0 : -6.0;
  tempHighBound1 = (configs.fixed_target) ? -2.5 :  6.0;
  TH2D *h2_phi_vs_eta_EPD = new TH2D("h2_phi_vs_eta_EPD", "EPD;#eta;#phi", tempBins1, tempLowBound1, tempHighBound1, 300, -4, 4);


  if (configs.sqrt_s_NN == 3.0 || configs.sqrt_s_NN == 3.22 || configs.sqrt_s_NN == 3.9)
    {
      tempBins1 = 300;
      tempLowBound1 = -1.2;
      tempHighBound1 = 1.2;
      tempBins2 = 300;
      tempLowBound2  = 0.0;
      tempHighBound2 = 2.5;
    }
  else if (configs.sqrt_s_NN == 7.2)
    {
      tempBins1 = 300;
      tempLowBound1 = -0.2;
      tempHighBound1 = 2.2;
      tempBins2 = 300;
      tempLowBound2  = 0.0;
      tempHighBound2 = 2.5;
    }
  else if (configs.sqrt_s_NN == 19.6)
    {
      tempBins1 = 400;
      tempLowBound1 = -2.0;
      tempHighBound1 = 2.0;
      tempBins2 = 500;
      tempLowBound2  = 0.0;
      tempHighBound2 = 5.0;
    }
  else
    {
      tempBins1 = 0;
      tempLowBound1 = 0.0;
      tempHighBound1 = 0.0;
      tempBins2 = 0;
      tempLowBound2  = 0.0;
      tempHighBound2 = 0.0;
    }
  TH2D *h2_pToverA_vs_yCM_de = 
    new TH2D("h2_pToverA_vs_yCM_de","Deuteron;y-y_{mid};p_{T}/A (GeV/c)",tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_pToverA_vs_yCM_tr = 
    new TH2D("h2_pToverA_vs_yCM_tr", "Triton;y-y_{mid};p_{T}/A (GeV/c)", tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);

  TH2D *h2_KToverA_vs_yCM_pr = 
    new TH2D("h2_KToverA_vs_yCM_pr", "Proton;y-y_{mid};m_{T}-m_{0} (GeV)",  tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_KToverA_vs_yCM_de = 
    new TH2D("h2_KToverA_vs_yCM_de", "Deuteron;y-y_{mid};m_{T}-m_{0} (GeV)",tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  TH2D *h2_KToverA_vs_yCM_tr = 
    new TH2D("h2_KToverA_vs_yCM_tr", "Triton;y-y_{mid};m_{T}-m_{0} (GeV)",  tempBins1, tempLowBound1, tempHighBound1, tempBins2, tempLowBound2, tempHighBound2);
  /////

  // Here the name refers to the eta region that will be displayed/searched using the event plane angle from the opposite region
  tempBins1      = (configs.fixed_target) ?  50  :  100;
  tempLowBound1  = (configs.fixed_target) ? -2.0 : -2.3;
  tempHighBound1 = (configs.fixed_target) ?  0.0 :  2.3;
  TProfile2D *h2_vnScanTpc = new TProfile2D("h2_vnScanTpc", "<cos("+ORDER_N_STR+"(#phi^{TPC} - #psi^{EPD}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
					    tempBins1, tempLowBound1, tempHighBound1, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_vnScanTpcEpdA = new TProfile2D("h2_vnScanTpcEpdA", "<cos("+ORDER_N_STR+"(#phi^{TPC} - #psi^{EPD,A}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
						tempBins1, tempLowBound1, tempHighBound1, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_vnScanTpcEpdB = new TProfile2D("h2_vnScanTpcEpdB", "<cos("+ORDER_N_STR+"(#phi^{TPC} - #psi^{EPD,B}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
						tempBins1, tempLowBound1, tempHighBound1, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_vnScanEpd = new TProfile2D("h2_vnScanEpd", "<cos("+ORDER_N_STR+"(#phi^{EPD} - #psi^{TPC}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
					    50, -5.2, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_vnScanEpdTpcA = new TProfile2D("h2_vnScanEpdTpcA", "<cos("+ORDER_N_STR+"(#phi^{EPD} - #psi^{TPC,A}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
						50, -5.2, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_vnScanEpdTpcB = new TProfile2D("h2_vnScanEpdTpcB", "<cos("+ORDER_N_STR+"(#phi^{EPD} - #psi^{TPC,B}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
						50, -5.2, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  h2_vnScanTpc->SetStats(0);
  h2_vnScanTpcEpdA->SetStats(0);
  h2_vnScanTpcEpdB->SetStats(0);
  h2_vnScanEpd->SetStats(0);
  h2_vnScanEpdTpcA->SetStats(0);
  h2_vnScanEpdTpcB->SetStats(0);
  ////

  // Raw correlations between event planes
  TH2D *h2_psiEpdATpcA = new TH2D("h2_psiEpdATpcA", "#psi^{EPD,A} vs #psi^{TPC,A} (Order "+ORDER_M_STR+");#psi^{TPC}_{A};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBTpcA = new TH2D("h2_psiEpdBTpcA", "#psi^{EPD,B} vs #psi^{TPC,A} (Order "+ORDER_M_STR+");#psi^{TPC}_{A};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiEpdATpcB = new TH2D("h2_psiEpdATpcB", "#psi^{EPD,A} vs #psi^{TPC,B} (Order "+ORDER_M_STR+");#psi^{TPC}_{B};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBTpcB = new TH2D("h2_psiEpdBTpcB", "#psi^{EPD,B} vs #psi^{TPC,B} (Order "+ORDER_M_STR+");#psi^{TPC}_{B};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiTpcATpcB = new TH2D("h2_psiTpcATpcB", "#psi^{TPC,A} vs #psi^{TPC,B} (Order "+ORDER_M_STR+");#psi^{TPC}_{B};#psi^{TPC}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiEpdAEpdB = new TH2D("h2_psiEpdAEpdB", "#psi^{EPD,A} vs #psi^{EPD,B} (Order "+ORDER_M_STR+");#psi^{EPD}_{B};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  ////


  TH1D *h_XnTpc  = new TH1D("h_XnTpc", "X_n Distribution (TPC);X_n;Events",    250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpc  = new TH1D("h_YnTpc", "Y_n Distribution (TPC);Y_n;Events",    250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcA = new TH1D("h_XnTpcA", "X_n Distribution (TPC A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcA = new TH1D("h_YnTpcA", "Y_n Distribution (TPC A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcB = new TH1D("h_XnTpcB", "X_n Distribution (TPC B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcB = new TH1D("h_YnTpcB", "Y_n Distribution (TPC B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpd  = new TH1D("h_XnEpd", "X_n Distribution (EPD);X_n;Events",    250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpd  = new TH1D("h_YnEpd", "Y_n Distribution (EPD);Y_n;Events",    250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdA = new TH1D("h_XnEpdA", "X_n Distribution (EPD A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdA = new TH1D("h_YnEpdA", "Y_n Distribution (EPD A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdB = new TH1D("h_XnEpdB", "X_n Distribution (EPD B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdB = new TH1D("h_YnEpdB", "Y_n Distribution (EPD B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);

  // CORRECTION HISTOGRAMS
  TProfile *p_sinAvgsTpc  = new TProfile("p_sinAvgsTpc", "Sin Averages (TPC);j (Correction term);<sin(jn#psi^{TPC}_{n})>",      configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_cosAvgsTpc  = new TProfile("p_cosAvgsTpc", "Cos Averages (TPC);j (Correction term);<sin(jn#psi^{TPC}_{n})>",      configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_sinAvgsTpcA = new TProfile("p_sinAvgsTpcA", "Sin Averages (TPC A);j (Correction term);<sin(jn#psi^{TPC,A}_{n})>", configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_cosAvgsTpcA = new TProfile("p_cosAvgsTpcA", "Cos Averages (TPC A);j (Correction term);<sin(jn#psi^{TPC,A}_{n})>", configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_sinAvgsTpcB = new TProfile("p_sinAvgsTpcB", "Sin Averages (TPC B);j (Correction term);<sin(jn#psi^{TPC,B}_{n})>", configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_cosAvgsTpcB = new TProfile("p_cosAvgsTpcB", "Cos Averages (TPC B);j (Correction term);<sin(jn#psi^{TPC,B}_{n})>", configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_sinAvgsEpd  = new TProfile("p_sinAvgsEpd", "Sin Averages (EPD);j (Correction term);<sin(jn#psi^{EPD}_{n})>",      configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_cosAvgsEpd  = new TProfile("p_cosAvgsEpd", "Cos Averages (EPD);j (Correction term);<sin(jn#psi^{EPD}_{n})>",      configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_sinAvgsEpdA = new TProfile("p_sinAvgsEpdA", "Sin Averages (EPD A);j (Correction term);<sin(jn#psi^{EPD,A}_{n})>", configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_cosAvgsEpdA = new TProfile("p_cosAvgsEpdA", "Cos Averages (EPD A);j (Correction term);<sin(jn#psi^{EPD,A}_{n})>", configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_sinAvgsEpdB = new TProfile("p_sinAvgsEpdB", "Sin Averages (EPD B);j (Correction term);<sin(jn#psi^{EPD,B}_{n})>", configs.shift_terms, 1, configs.shift_terms + 1);
  TProfile *p_cosAvgsEpdB = new TProfile("p_cosAvgsEpdB", "Cos Averages (EPD B);j (Correction term);<sin(jn#psi^{EPD,B}_{n})>", configs.shift_terms, 1, configs.shift_terms + 1);

  // RECENTERED (RC) HISTOGRAMS
  TH1D *h_XnTpc_RC  = new TH1D("h_XnTpc_RC", "Re-centered X_n Distribution (TPC);X_n;Events",    200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpc_RC  = new TH1D("h_YnTpc_RC", "Re-centered Y_n Distribution (TPC);Y_n;Events",    200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcA_RC = new TH1D("h_XnTpcA_RC", "Re-centered X_n Distribution (TPC A);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcA_RC = new TH1D("h_YnTpcA_RC", "Re-centered Y_n Distribution (TPC A);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcB_RC = new TH1D("h_XnTpcB_RC", "Re-centered X_n Distribution (TPC B);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcB_RC = new TH1D("h_YnTpcB_RC", "Re-centered Y_n Distribution (TPC B);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpd_RC  = new TH1D("h_XnEpd_RC", "Re-centered X_n Distribution (EPD);X_n;Events",    200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpd_RC  = new TH1D("h_YnEpd_RC", "Re-centered Y_n Distribution (EPD);Y_n;Events",    200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdA_RC = new TH1D("h_XnEpdA_RC", "Re-centered X_n Distribution (EPD A);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdA_RC = new TH1D("h_YnEpdA_RC", "Re-centered Y_n Distribution (EPD A);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdB_RC = new TH1D("h_XnEpdB_RC", "Re-centered X_n Distribution (EPD B);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdB_RC = new TH1D("h_YnEpdB_RC", "Re-centered Y_n Distribution (EPD B);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);

  TH1D *h_psiTpc_RC  = new TH1D("h_psiTpc_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", TPC);#psi_{"+ORDER_M_STR+"};Events",    400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcA_RC = new TH1D("h_psiTpcA_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", TPC A);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcB_RC = new TH1D("h_psiTpcB_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", TPC B);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpd_RC  = new TH1D("h_psiEpd_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", EPD);#psi_{"+ORDER_M_STR+"};Events",    400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdA_RC = new TH1D("h_psiEpdA_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", EPD A);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdB_RC = new TH1D("h_psiEpdB_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", EPD B);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  // RECENTERED AND SHIFTED HISTOGRAMS
  TH1D *h_psiTpc_FLAT  = new TH1D("h_psiTpc_FLAT", "Flattened Event Plane Angle (TPC, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events",    400, -PSI_BOUNDS, PSI_BOUNDS);      
  TH1D *h_psiTpcA_FLAT = new TH1D("h_psiTpcA_FLAT", "Flattened Event Plane Angle (TPC A, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcB_FLAT = new TH1D("h_psiTpcB_FLAT", "Flattened Event Plane Angle (TPC B, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpd_FLAT  = new TH1D("h_psiEpd_FLAT", "Flattened Event Plane Angle (EPD, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events",    400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdA_FLAT = new TH1D("h_psiEpdA_FLAT", "Flattened Event Plane Angle (EPD A, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdB_FLAT = new TH1D("h_psiEpdB_FLAT", "Flattened Event Plane Angle (EPD B, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  
  FlowUtils::Event eventInfo;        // These hold info for each event
  FlowUtils::Particle particleInfo;  // These hold info for each track/hit
  StEpdGeom *epdGeom = new StEpdGeom();

  // EVENT LOOP
  Int_t events2read = tree->GetEntries();
  std::cout << "Setup complete, beginning analysis on " << events2read << " events..." << std::endl;
  for (Long64_t ievent = 0; ievent < events2read; ievent++)
    {
      eventInfo.reset();

      tree->GetEntry(ievent);

      // At this point, all bad runs and bad trigger events are removed by TreeMaker.
      // Now check event vertex
      TVector3 pVtx(f_xvtx, f_yvtx, f_zvtx);
      Double_t d_rvtx = TMath::Sqrt(f_xvtx * f_xvtx + (f_yvtx + 2) * (f_yvtx + 2));
      
      h_zvtx->Fill(f_zvtx); // All events that pass minbias and wide TreeMaker cuts. Possibly move this past Vz cut.

      Bool_t goodVertexZ = (f_zvtx > configs.z_vtx_low && f_zvtx < configs.z_vtx_high);
      if (!goodVertexZ) continue;
      h_eventCheck->Fill(3);    // Count # of events after Vz cut
      h2_trans_vtx->Fill(pVtx.X(),pVtx.Y()); // transverse vtx position after z cut

      Bool_t goodVertexR = (d_rvtx < configs.r_vtx);
      if(!goodVertexR) continue;
      h_eventCheck->Fill(4);    // Count # of events after Vr cut
      h2_trans_vtx_cut->Fill(pVtx.X(),pVtx.Y()); // transverse vtx position after z cut and R cut

      eventInfo.centID = i_centrality;
      if (i_centrality == -99) continue;  // Remove undefined events
      h_eventCheck->Fill(5); // Count # of events after centrality cut
      h_centralities->Fill(i_centrality);

      Double_t d_px;
      Double_t d_py;
      Double_t d_pz;
      Double_t d_pT;
      Double_t d_mom;
      Double_t d_eta;
      Double_t d_phi;
      Double_t d_DCA;
      Double_t d_nSigmaPi;
      Double_t d_nSigmaKa;
      Double_t d_nSigmaPr;
      Double_t d_tofBeta;
      Double_t d_dEdx;
      Int_t i_nHits;
      Int_t i_nHitsFit;
      Int_t i_nHitsPoss;
      Int_t i_nHitsDedx;
      Short_t s_charge;

      Int_t N_pp = 0;
      Int_t N_pm = 0;
      Int_t N_kp = 0;
      Int_t N_km = 0;
      Int_t N_pr = 0;
      Int_t N_de = 0;
      Int_t N_tr = 0;


      // TRACK LOOP OVER PRIMARY TRACKS
      Int_t nTracks = (Int_t)i_trackNumber;
      for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
	{
	  particleInfo.reset();

	  eventInfo.primTracks++;
	  h_trackCheck->Fill(0);  // Only event cuts

	  d_px  = Px[iTrk];
	  d_py  = Py[iTrk];
	  d_pz  = Pz[iTrk];
	  d_phi = FlowUtils::phi(Px[iTrk], Py[iTrk]);
	  d_eta = FlowUtils::pseudorapidity(Px[iTrk], Py[iTrk], Pz[iTrk]);
	  d_pT  = FlowUtils::transMomentum(Px[iTrk], Py[iTrk]);
	  d_mom = FlowUtils::totalMomentum(Px[iTrk], Py[iTrk], Pz[iTrk]);
	  s_charge = charge[iTrk];
	  d_DCA = DCA[iTrk];
	  d_nSigmaPi = nSigmaPi[iTrk];
	  d_nSigmaKa = nSigmaKa[iTrk];
	  d_nSigmaPr = nSigmaPr[iTrk];
	  d_tofBeta = tofBeta[iTrk];
	  d_dEdx = dEdx[iTrk];
	  i_nHits = nHits[iTrk];
	  i_nHitsFit = nHitsFit[iTrk];
	  i_nHitsPoss = nHitsPoss[iTrk];
	  i_nHitsDedx = nHitsDedx[iTrk];


	  h_nhits->Fill(i_nHits);
	  h_nhitsFit->Fill(i_nHitsFit);
	  h_nhitsPoss->Fill(i_nHitsPoss);
	  h_nhits_ratio->Fill((double)i_nHitsFit/(double)i_nHitsPoss);
	  h_nhits_dEdx->Fill(i_nHitsDedx);
	  h_DCA->Fill(d_DCA);

	  //=========================================================
	  //          Track QA Cuts
	  //=========================================================
	  bool b_bad_hits  = ( i_nHits < configs.nHits );
	  bool b_bad_dEdx  = ( i_nHitsDedx <= configs.nHits_dEdx );
	  bool b_bad_ratio = ( ((double)i_nHitsFit / (double)i_nHitsPoss) <= configs.nHits_ratio );
	  bool b_bad_DCA   = ( d_DCA >= configs.dca );

	  if (b_bad_hits || b_bad_dEdx || b_bad_ratio || b_bad_DCA) continue;
	  //=========================================================
	  //          End Track QA Cuts
	  //=========================================================
	  h_trackCheck->Fill(1); // After QA Cuts

	  h_phi->Fill(d_phi);
	  h_eta->Fill(d_eta);
	  h_pT->Fill(d_pT);
	  if (s_charge <= 1 && s_charge >= -1) h2_dEdx_vs_qp->Fill(s_charge * d_mom, d_dEdx);
	  if (s_charge == 2 || s_charge == -2) h2_dEdx_vs_qp_charge2->Fill(s_charge * d_mom, d_dEdx);
	  if (s_charge == 1) h2_dEdx_vs_qp_half->Fill(s_charge * d_mom, d_dEdx);

	  // Get event planes from the TPC here before the TOF cut
	  if (s_charge != 0)
	    {
	      eventInfo.nTracksTpc++;

	      particleInfo.phi = d_phi;
	      particleInfo.eta = d_eta;
	      particleInfo.pT  = d_pT;
	      particleInfo.weight = (configs.sqrt_s_NN == 7.2) ? 1.0/*TMath::Abs(d_eta-Y_MID)*/ : d_pT;

	      h2_phi_vs_eta_TPC->Fill(d_eta, d_phi);

	      eventInfo.incrementQvectorTPC(ODD_PLANE, ORDER_M, Y_MID, particleInfo.eta, particleInfo.phi, particleInfo.weight);

	      // TPC A
	      if (d_eta > configs.tpc_A_low_eta && d_eta < configs.tpc_A_high_eta)
		{
		  eventInfo.nTracksTpcA++;
		  particleInfo.isInTpcA = true;
		  eventInfo.incrementQvectorTPCA(ODD_PLANE, ORDER_M, Y_MID, particleInfo.eta, particleInfo.phi, particleInfo.weight);
		}
	      // TPC B
	      if (d_eta > configs.tpc_B_low_eta && d_eta < configs.tpc_B_high_eta)
		{
		  eventInfo.nTracksTpcB++;
		  particleInfo.isInTpcB = true;
		  eventInfo.incrementQvectorTPCB(ODD_PLANE, ORDER_M, Y_MID, particleInfo.eta, particleInfo.phi, particleInfo.weight);
		}
	      
	      
	      // TOF information here before full PID next
	      // d_tofBeta = -999.0 if it's not a TOF hit
	      //=========================================================
	      //          TOF Beta Cuts
	      //=========================================================
	      Double_t d_m2 = -999.0;
	      Bool_t tofTrack = (d_tofBeta != -999.0);

	      if (tofTrack)
		{
		  d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0);		    
		  h_m2->Fill(d_m2);
		  h_tofBeta->Fill(d_tofBeta);		  
		  h2_beta_vs_qp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
		  h2_m2_vs_qp->Fill(s_charge*d_mom, d_m2);
		}
	      //=========================================================
	      //          End TOF Beta Cuts
	      //=========================================================
	      


	      Double_t d_zDeuteron = (s_charge == 1) ? TMath::Log(d_dEdx / bichselZ_de->Eval(d_mom)) : -999.0;
	      Double_t d_zTriton   = (s_charge == 1) ? TMath::Log(d_dEdx / bichselZ_tr->Eval(d_mom)) : -999.0;
	      Double_t d_zHelium3  = (s_charge == 2) ? TMath::Log(d_dEdx / bichselZ_he3->Eval(d_mom)) : -999.0;
	      Double_t d_zAlpha    = (s_charge == 2) ? TMath::Log(d_dEdx / bichselZ_al->Eval(d_mom)) : -999.0;

	      if (tofTrack && d_zHelium3 < 0.2 && d_zHelium3 > -0.2) { h_m2_alpha_he3->Fill(d_m2); }
	      if (tofTrack && d_zAlpha < 0.2 && d_zAlpha > -0.2) { h_m2_alpha_he3->Fill(d_m2); }

	      //=========================================================
	      //          PID Cuts
	      //=========================================================
	      Bool_t pion   = false;
	      Bool_t kaon   = false;
	      Bool_t proton = (d_nSigmaPr > configs.nSig_pr_low) && (d_nSigmaPr < configs.nSig_pr_high) && (s_charge == 1);
	      //Bool_t proton = false;
	      Bool_t deuteron = false;
	      Bool_t triton   = false;
	      //Bool_t deuteron = (d_zDeuteron > configs.z_de_low) && (d_zDeuteron < configs.z_de_high);
	      //Bool_t triton   = (d_zTriton > configs.z_tr_low) && (d_zTriton < configs.z_tr_high);

	      if (tofTrack)
		{
		  pion = (d_nSigmaPi > configs.nSig_pi_low) &&
		    (d_nSigmaPi < configs.nSig_pi_high) &&
		    (d_m2 > configs.m2_pi_low) &&
		    (d_m2 < configs.m2_pi_high);

		  kaon = (d_nSigmaKa > configs.nSig_ka_low) &&
		    (d_nSigmaKa < configs.nSig_ka_high) &&
		    (d_m2 > configs.m2_ka_low) &&
		    (d_m2 < configs.m2_ka_high);
		  /*
		  deuteron = (d_zDeuteron > configs.z_de_low) &&
		    (d_zDeuteron < configs.z_de_high) &&
		    (d_m2 > configs.m2_de_low) &&
		    (d_m2 < configs.m2_de_high);

		  triton   = (d_zTriton > configs.z_tr_low) &&
		  (d_zTriton < configs.z_tr_high) &&
		    (d_m2 > configs.m2_tr_low) &&
		    (d_m2 < configs.m2_tr_high);
		  */
		}

	      // 3.0 GeV d and t PID
	      if (!pion && !kaon && configs.sqrt_s_NN == 3.0)
		{
		  deuteron = FlowUtils::momDepDeuteronID(configs.sqrt_s_NN, d_mom, d_zDeuteron, tofTrack, d_m2, 
		  					 configs.z_de_low, configs.z_de_high, configs.m2_de_low, configs.m2_de_high);
		  //deuteron = FlowUtils::momDepDeuteronID_lowSystematics(configs.sqrt_s_NN, d_mom, d_zDeuteron, tofTrack, d_m2, 
		  //							configs.z_de_low, configs.z_de_high, configs.m2_de_low, configs.m2_de_high);
		  //deuteron = FlowUtils::momDepDeuteronID_highSystematics(configs.sqrt_s_NN, d_mom, d_zDeuteron, tofTrack, d_m2, 
		  //							 configs.z_de_low, configs.z_de_high, configs.m2_de_low, configs.m2_de_high);

		  triton = FlowUtils::momDepTritonID(configs.sqrt_s_NN, d_mom, d_zTriton, tofTrack, d_m2, 
		  				     configs.z_tr_low, configs.z_tr_high, configs.m2_tr_low, configs.m2_tr_high);
		  //triton = FlowUtils::momDepTritonID_lowSystematics(configs.sqrt_s_NN, d_mom, d_zTriton, tofTrack, d_m2, 
		  //						    configs.z_tr_low, configs.z_tr_high, configs.m2_tr_low, configs.m2_tr_high);
		  //triton = FlowUtils::momDepTritonID_highSystematics(configs.sqrt_s_NN, d_mom, d_zTriton, tofTrack, d_m2, 
		  //						     configs.z_tr_low, configs.z_tr_high, configs.m2_tr_low, configs.m2_tr_high);

		}
	      // Basic d and t PID
	      else if (!pion && !kaon && (configs.sqrt_s_NN == 7.2 || configs.sqrt_s_NN == 3.22 || configs.sqrt_s_NN == 3.9))
		{
		  //  DEUTERON
		  if (tofTrack)
		    {
		      if (d_zDeuteron > configs.z_de_low &&
			  d_zDeuteron < configs.z_de_high &&
			  d_m2 > configs.m2_de_low &&
			  d_m2 < configs.m2_de_high)
			deuteron = true;
		    }

		  //  TRITON
		  if (tofTrack)
		    {
		      if (d_zTriton > configs.z_tr_low &&
			  d_zTriton < configs.z_tr_high &&
			  d_m2 > configs.m2_tr_low &&
			  d_m2 < configs.m2_tr_high)
			triton = true;
		    }
		}

	      /*
		if (deuteron && proton) 
		{ 
		if (TMath::Abs(d_zDeuteron) < TMath::Abs(d_nSigmaPr)) { proton = false; }
		else if (TMath::Abs(d_zDeuteron) == TMath::Abs(d_nSigmaPr)) { proton = false; deuteron = false; }
		else { deuteron = false; }
		}
		if (triton && proton)
		{
		if (TMath::Abs(d_zTriton) < TMath::Abs(d_nSigmaPr)) { proton = false; }
		else if (TMath::Abs(d_zTriton) == TMath::Abs(d_nSigmaPr)) { proton = false; triton = false; }
		else { triton = false; }
		}
		if (deuteron && triton)
		{
		if (TMath::Abs(d_zDeuteron) < TMath::Abs(d_zTriton)) { triton = false; }
		else if (TMath::Abs(d_zDeuteron) == TMath::Abs(d_zTriton)) { triton = false; deuteron = false; }
		else { deuteron = false; }
		}
	      */

	      if (deuteron && triton) { deuteron = false; triton = false; } // Ignore tags of both d and t.
	      if (deuteron && proton) { proton = false; } // d and t will have some contamination from p, but that has been minimized
	      if (triton && proton)   { proton = false; }

	      if (pion && proton)   { proton = false; }
	      if (kaon && proton)   { proton = false; }
	      //=========================================================
	      //          END PID Cuts
	      //=========================================================

	      if (pion || kaon || proton || deuteron || triton) h_trackCheck->Fill(2);

	      if (!pion && !kaon) 
		{
		  h2_nSigp_vs_mom->Fill(d_mom, d_nSigmaPr);
		  h2_zd_vs_mom->Fill(d_mom, d_zDeuteron);
		  h2_zt_vs_mom->Fill(d_mom, d_zTriton);
		}


	      Double_t d_rapidity = 999.0;
	      Double_t d_mT = -999.0;
	      Double_t massNumber;

	      if(pion) // PID Pions
		{ 
		  d_rapidity = FlowUtils::rapidity(d_px, d_py, d_pz, D_M0_PI);
		  d_mT = FlowUtils::transMass(d_px, d_py, D_M0_PI);

		  particleInfo.rapidity = d_rapidity;
		  particleInfo.KT = d_mT - D_M0_PI;

		  if(s_charge > 0)
		    {
		      particleInfo.ppTag = true;
		      N_pp++;
		      h_eta_pp->Fill(d_eta);
		      h_phi_pp->Fill(d_phi);
		      h_pT_pp->Fill(d_pT);
		      h_mom_pp->Fill(d_mom);
		      h_dndy_pp->Fill(d_rapidity);
		      h_dndm_pp->Fill(d_mT - D_M0_PI);
		      h2_pT_vs_yCM_pp->Fill(d_rapidity - Y_MID, d_pT);
		      h2_dEdx_vs_qp_pp->Fill(s_charge*d_mom, d_dEdx);
		      h2_beta_vs_qp_pp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
		      h2_m2_vs_qp_pp->Fill(s_charge*d_mom, d_m2);
			  
		      if (d_rapidity - Y_MID > configs.yCM_norm_pi_low && d_rapidity - Y_MID < configs.yCM_norm_pi_high && 
			  d_pT >= configs.pt_norm_pi_low && d_pT <= configs.pt_norm_pi_high)
			{
			  //h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  //h2_y_vs_eta_pp->Fill(d_eta, d_rapidity);
			  p_meanpT_vs_yCM_pp->Fill(d_rapidity - Y_MID, particleInfo.pT);
			}
		    }
		  else if(s_charge < 0)
		    {
		      particleInfo.pmTag = true;
		      N_pm++;
		      h_eta_pm->Fill(d_eta);
		      h_phi_pm->Fill(d_phi);
		      h_pT_pm->Fill(d_pT);
		      h_mom_pm->Fill(d_mom);
		      h_dndy_pm->Fill(d_rapidity);
		      h_dndm_pm->Fill(d_mT - D_M0_PI);
		      h2_pT_vs_yCM_pm->Fill(d_rapidity - Y_MID, d_pT);
		      h2_dEdx_vs_qp_pm->Fill(s_charge*d_mom, d_dEdx);
		      h2_beta_vs_qp_pm->Fill(s_charge*d_mom, 1.0/d_tofBeta);
		      h2_m2_vs_qp_pm->Fill(s_charge*d_mom, d_m2);

		      if (d_rapidity - Y_MID > configs.yCM_norm_pi_low && d_rapidity - Y_MID < configs.yCM_norm_pi_high && 
			  d_pT >= configs.pt_norm_pi_low && d_pT <= configs.pt_norm_pi_high)
			{
			  //h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  //h2_y_vs_eta_pm->Fill(d_eta, d_rapidity);
			  p_meanpT_vs_yCM_pm->Fill(d_rapidity - Y_MID, particleInfo.pT);
			}
		    }
		}
	      else if(kaon) // PID Kaons
		{
		  d_rapidity = FlowUtils::rapidity(d_px, d_py, d_pz, D_M0_KA);
		  d_mT = FlowUtils::transMass(d_px, d_py, D_M0_KA);

		  particleInfo.rapidity = d_rapidity;
		  particleInfo.KT = d_mT - D_M0_KA;
				  
		  if(s_charge > 0)
		    {
		      particleInfo.kpTag = true;
		      N_kp++;
		      h_eta_kp->Fill(d_eta);
		      h_phi_kp->Fill(d_phi);
		      h_pT_kp->Fill(d_pT);
		      h_mom_kp->Fill(d_mom);
		      h_dndy_kp->Fill(d_rapidity);
		      h_dndm_kp->Fill(d_mT - D_M0_KA);
		      h2_pT_vs_yCM_kp->Fill(d_rapidity - Y_MID, d_pT);
		      h2_dEdx_vs_qp_kp->Fill(s_charge*d_mom, d_dEdx);
		      h2_beta_vs_qp_kp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
		      h2_m2_vs_qp_kp->Fill(s_charge*d_mom, d_m2);

		      if (d_rapidity - Y_MID > configs.yCM_norm_ka_low && d_rapidity - Y_MID < configs.yCM_norm_ka_high && 
			  d_pT >= configs.pt_norm_ka_low && d_pT <= configs.pt_norm_ka_high)
			{
			  //h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  //h2_y_vs_eta_kp->Fill(d_eta, d_rapidity);
			  p_meanpT_vs_yCM_kp->Fill(d_rapidity - Y_MID, particleInfo.pT);
			}

		    }
		  else if(s_charge < 0)
		    {
		      particleInfo.kmTag = true;
		      N_km++;
		      h_eta_km->Fill(d_eta);
		      h_phi_km->Fill(d_phi);
		      h_pT_km->Fill(d_pT);
		      h_mom_km->Fill(d_mom);
		      h_dndy_km->Fill(d_rapidity);
		      h_dndm_km->Fill(d_mT - D_M0_KA);
		      h2_pT_vs_yCM_km->Fill(d_rapidity - Y_MID, d_pT);
		      h2_dEdx_vs_qp_km->Fill(s_charge*d_mom, d_dEdx);
		      h2_beta_vs_qp_km->Fill(s_charge*d_mom, 1.0/d_tofBeta);
		      h2_m2_vs_qp_km->Fill(s_charge*d_mom, d_m2);

		      if (d_rapidity - Y_MID > configs.yCM_norm_ka_low && d_rapidity - Y_MID < configs.yCM_norm_ka_high && 
			  d_pT >= configs.pt_norm_ka_low && d_pT <= configs.pt_norm_ka_high)
			{
			  //h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  //h2_y_vs_eta_km->Fill(d_eta, d_rapidity);
			  p_meanpT_vs_yCM_km->Fill(d_rapidity - Y_MID, particleInfo.pT);
			}
		    }
		}
	      else if(proton) // PID Proton
		{
		  particleInfo.prTag = true;
		  N_pr++;
		  d_rapidity = FlowUtils::rapidity(d_px, d_py, d_pz, D_M0_PR);
		  d_mT = FlowUtils::transMass(d_px, d_py, D_M0_PR);
		  massNumber = 1.0;

		  particleInfo.rapidity = d_rapidity;
		  particleInfo.KT = d_mT - D_M0_PR;
				  
		  h_eta_pr->Fill(d_eta);
		  h_phi_pr->Fill(d_phi);
		  h_pT_pr->Fill(d_pT);
		  h_mom_pr->Fill(d_mom);
		  h_dndy_pr->Fill(d_rapidity);
		  h_dndm_pr->Fill(d_mT - D_M0_PR);
		  h2_pT_vs_yCM_pr->Fill(d_rapidity - Y_MID, d_pT);
		  h2_KToverA_vs_yCM_pr->Fill(d_rapidity - Y_MID, (d_mT - D_M0_PR)/1.0);
		  h2_dEdx_vs_qp_pr->Fill(s_charge*d_mom, d_dEdx);
		  //h2_beta_vs_qp_pr->Fill(s_charge*d_mom, 1.0/d_tofBeta);
		  //h2_m2_vs_qp_pr->Fill(s_charge*d_mom, d_m2);

		  // Normal acceptance region
		  if (d_rapidity - Y_MID > configs.yCM_norm_pr_low && d_rapidity - Y_MID < configs.yCM_norm_pr_high && 
		      d_pT >= configs.pt_norm_pr_low && d_pT <= configs.pt_norm_pr_high)
		    {
		      p_meanpT_vs_yCM_pr->Fill(d_rapidity - Y_MID, d_pT);
		      h2_dEdx_vs_qp_id_pr->Fill(d_mom, d_dEdx);
		      //h2_y_vs_eta->Fill(d_eta, d_rapidity);
		      //h2_y_vs_eta_pr->Fill(d_eta, d_rapidity);
		    }
		  // Alternate acceptance region
		  if (d_rapidity - Y_MID > configs.yCM_alt_pr_low && d_rapidity - Y_MID < configs.yCM_alt_pr_high && 
		      (d_mT-D_M0_PR)/massNumber >= configs.KT_pdt_low && (d_mT-D_M0_PR)/massNumber <= configs.KT_pdt_high)
		    {
		      h2_dEdx_vs_qp_id_pr_alt->Fill(d_mom, d_dEdx);
		      p_meanpT_vs_yCM_pr_alt->Fill(d_rapidity - Y_MID, d_pT);
		      h2_pT_vs_cent_pr->Fill(i_centrality, d_pT);
		      h2_pT_vs_yCM_pr_alt->Fill(d_rapidity - Y_MID, d_pT);
		      h_DCA_pr_alt->Fill(d_DCA);
		    }
		}		    
	      else if(deuteron) // PID Deuteron
		{
		  particleInfo.deTag = true;
		  N_de++;
		  d_rapidity = FlowUtils::rapidity(d_px, d_py, d_pz, D_M0_DE);
		  d_mT = FlowUtils::transMass(d_px, d_py, D_M0_DE);
		  massNumber = 2.0;

		  particleInfo.rapidity = d_rapidity;
		  particleInfo.KT = d_mT - D_M0_DE;
		  
		  h_eta_de->Fill(d_eta);
		  h_phi_de->Fill(d_phi);
		  h_pT_de->Fill(d_pT);
		  h_mom_de->Fill(d_mom);
		  h_dndy_de->Fill(d_rapidity);
		  h_dndm_de->Fill(d_mT - D_M0_DE);
		  h2_pT_vs_yCM_de->Fill(d_rapidity - Y_MID, d_pT);
		  h2_dEdx_vs_qp_de->Fill(s_charge*d_mom, d_dEdx);
		  //h2_beta_vs_qp_de->Fill(s_charge*d_mom, 1.0/d_tofBeta);
		  //h2_m2_vs_qp_de->Fill(s_charge*d_mom, d_m2);
		  h2_pToverA_vs_yCM_de->Fill(d_rapidity - Y_MID, d_pT/2.0);
		  h2_KToverA_vs_yCM_de->Fill(d_rapidity - Y_MID, (d_mT - D_M0_DE)/2.0);
		  //h2_dEdx_vs_qp_half_postZdCut->Fill(s_charge * d_mom, d_dEdx);

		  if (d_rapidity - Y_MID > configs.yCM_norm_de_low && d_rapidity - Y_MID < configs.yCM_norm_de_high && 
		      (d_mT-D_M0_DE)/massNumber >= configs.KT_pdt_low && (d_mT-D_M0_DE)/massNumber <= configs.KT_pdt_high)
		    {
		      //h2_y_vs_eta->Fill(d_eta, d_rapidity);
		      //h2_y_vs_eta_de->Fill(d_eta, d_rapidity);
		      p_meanpT_vs_yCM_de->Fill(d_rapidity - Y_MID, d_pT);
		      h2_dEdx_vs_qp_id_de->Fill(d_mom, d_dEdx);
		      h2_pT_vs_cent_de->Fill(i_centrality, d_pT);
		      h_DCA_de->Fill(d_DCA);
		    }
		}		    
	      else if(triton) // PID Triton
		{
		  particleInfo.trTag = true;
		  N_tr++;
		  d_rapidity = FlowUtils::rapidity(d_px, d_py, d_pz, D_M0_TR);
		  d_mT = FlowUtils::transMass(d_px, d_py, D_M0_TR);
		  massNumber = 3.0;

		  particleInfo.rapidity = d_rapidity;
		  particleInfo.KT = d_mT - D_M0_TR;
				  
		  h_eta_tr->Fill(d_eta);
		  h_phi_tr->Fill(d_phi);
		  h_pT_tr->Fill(d_pT);
		  h_mom_tr->Fill(d_mom);
		  h_dndy_tr->Fill(d_rapidity);
		  h_dndm_tr->Fill(d_mT - D_M0_TR);
		  h2_pT_vs_yCM_tr->Fill(d_rapidity - Y_MID, d_pT);
		  h2_dEdx_vs_qp_tr->Fill(s_charge*d_mom, d_dEdx);
		  //h2_beta_vs_qp_tr->Fill(s_charge*d_mom, 1.0/d_tofBeta);
		  //h2_m2_vs_qp_tr->Fill(s_charge*d_mom, d_m2);
		  h2_pToverA_vs_yCM_tr->Fill(d_rapidity - Y_MID, d_pT/3.0);
		  h2_KToverA_vs_yCM_tr->Fill(d_rapidity - Y_MID, (d_mT - D_M0_TR)/3.0);
		  //h2_dEdx_vs_qp_half_postZtCut->Fill(s_charge * d_mom, d_dEdx);

		  if (d_rapidity - Y_MID > configs.yCM_norm_tr_low && d_rapidity - Y_MID < configs.yCM_norm_tr_high && 
		      (d_mT-D_M0_TR)/massNumber >= configs.KT_pdt_low && (d_mT-D_M0_TR)/massNumber <= configs.KT_pdt_high)
		    {
		      //h2_y_vs_eta->Fill(d_eta, d_rapidity);
		      //h2_y_vs_eta_tr->Fill(d_eta, d_rapidity);
		      p_meanpT_vs_yCM_tr->Fill(d_rapidity - Y_MID, d_pT);
		      h2_dEdx_vs_qp_id_tr->Fill(d_mom, d_dEdx);
		      h2_pT_vs_cent_tr->Fill(i_centrality, d_pT);
		      h_DCA_tr->Fill(d_DCA);
		    }
		}		    

	      eventInfo.tpcParticles.push_back(particleInfo);
	    }// End if(s_charge != 0)
	}//End TPC track loop
      particleInfo.reset();

      h_mult_pp->Fill(N_pp);
      h_mult_pm->Fill(N_pm);
      h_mult_kp->Fill(N_kp);
      h_mult_km->Fill(N_km);
      h_mult_pr->Fill(N_pr);
      h_mult_de->Fill(N_de);
      h_mult_tr->Fill(N_tr);

      //=========================================================
      //                EPD STUFF
      //=========================================================
      Short_t tileID;
      TVector3 tileVector;     // Vector from vertex to center of tile that was hit
      Int_t tileSector;
      Int_t tileRow;
      Double_t tileWeight;
      Double_t tileEta;
      Double_t tilePhi;
      Double_t tilenMip;

      FlowUtils::Particle epdParticleInfo;
      for (UShort_t iEpdHit = 0; iEpdHit < i_nEPDhits; iEpdHit++)
	{
	  epdParticleInfo.reset();

	  tileID = EPDids[iEpdHit];
	  Bool_t epdAside = (tileID < 0);  // East side
	  Bool_t epdBside = (configs.fixed_target) ? (tileID < 0) : (tileID > 0); // EPD B is also on the East side in FXT
	  
	  tileVector = epdGeom->TileCenter(tileID) - pVtx;
	  tileSector = FlowUtils::epdSector(tileID);
	  tileRow = FlowUtils::epdRow(tileID);
	  tileEta = tileVector.Eta();
	  tilePhi = tileVector.Phi();
	  tilenMip = EPDnMip[iEpdHit];
	  tileWeight = (tilenMip > configs.epd_threshold) ? ( (tilenMip > configs.epd_max_weight)?configs.epd_max_weight:tilenMip ) : 0;

	  h_eta->Fill(tileEta);
	  //if (configs.sqrt_s_NN == 7.2) tileWeight *= TMath::Abs(tileEta - Y_MID);

	  if (epdAside)
	    {
	      eventInfo.nHitsEpd++;
	      h2_ring_vs_eta->Fill(tileEta, tileRow);

	      eventInfo.incrementQvectorEPD(ODD_PLANE, ORDER_M, Y_MID, tileEta, tilePhi, tileWeight);
	    }

	  if (epdAside && tileRow >= configs.epdA_inner_row && tileRow <= configs.epdA_outer_row)
	    {
	      eventInfo.nHitsEpdA++;
	      epdParticleInfo.isInEpdA = true;
	      epdParticleInfo.phi    = tilePhi;
	      epdParticleInfo.eta    = tileEta;
	      epdParticleInfo.weight = tileWeight;

	      h_tileWeights->Fill(tileWeight);
	      h2_phi_vs_eta_EPD->Fill(tileEta, tilePhi);
	      p2_pp_vs_eta->Fill(tileEta, tileSector, tileWeight);

	      eventInfo.incrementQvectorEPDA(ODD_PLANE, ORDER_M, Y_MID, tileEta, tilePhi, tileWeight);
	      eventInfo.epdParticles.push_back(epdParticleInfo);
	    }
	  else if (epdBside && tileRow >= configs.epdB_inner_row && tileRow <= configs.epdB_outer_row)
	    {
	      eventInfo.nHitsEpdB++;
	      epdParticleInfo.isInEpdB = true;
	      epdParticleInfo.phi    = tilePhi;
	      epdParticleInfo.eta    = tileEta;
	      epdParticleInfo.weight = tileWeight;

	      h_tileWeights->Fill(tileWeight);
	      h2_phi_vs_eta_EPD->Fill(tileEta, tilePhi);
	      p2_pp_vs_eta->Fill(tileEta, tileSector, tileWeight);

	      if (configs.fixed_target) h2_ring_vs_eta->Fill(tileEta, tileRow);  // Only need the East side for this plot

	      eventInfo.incrementQvectorEPDB(ODD_PLANE, ORDER_M, Y_MID, tileEta, tilePhi, tileWeight);
	      eventInfo.epdParticles.push_back(epdParticleInfo);
	    }
	}// End EPD hit loop
      epdParticleInfo.reset();
      //=========================================================
      //            END EPD STUFF
      //=========================================================

      //if (eventInfo.nTracksTpcA < configs.min_tracks) continue;
      if (eventInfo.nTracksTpcB < configs.min_tracks) continue;
      //if (eventInfo.nHitsEpd    < configs.min_tracks) continue;
      if (eventInfo.nHitsEpdA   < configs.min_tracks) continue;
      //if (eventInfo.nHitsEpdB   < configs.min_tracks) continue;
      if (configs.fixed_target && configs.sqrt_s_NN == 3.0 && eventInfo.nHitsEpdB < configs.min_tracks+4) continue;
      else if (configs.fixed_target && configs.sqrt_s_NN == 3.22 && eventInfo.nHitsEpdB < configs.min_tracks+4) continue;
      else if (configs.fixed_target && configs.sqrt_s_NN == 3.9  && eventInfo.nHitsEpdB < configs.min_tracks+4) continue;

      h_centralities_final->Fill(i_centrality);

      FlowUtils::checkZeroQ(eventInfo);  // Remove events with no flow
      if (eventInfo.badEvent) continue;

      FlowUtils::getAllPsi(eventInfo, ORDER_M);

      h2_hits_vs_cent_EpdA->Fill(eventInfo.centID, eventInfo.nHitsEpdA);
      h2_hits_vs_cent_EpdB->Fill(eventInfo.centID, eventInfo.nHitsEpdB);
      h2_hits_vs_cent_TpcB->Fill(eventInfo.centID, eventInfo.nTracksTpcB);

      h_primTracks->Fill(eventInfo.primTracks);      

      h_XnTpc->Fill(eventInfo.XnTpc);
      h_YnTpc->Fill(eventInfo.YnTpc);
      h_XnTpcA->Fill(eventInfo.XnTpcA);
      h_YnTpcA->Fill(eventInfo.YnTpcA);
      h_XnTpcB->Fill(eventInfo.XnTpcB);
      h_YnTpcB->Fill(eventInfo.YnTpcB);
      h_XnEpd->Fill(eventInfo.XnEpd);
      h_YnEpd->Fill(eventInfo.YnEpd);
      h_XnEpdA->Fill(eventInfo.XnEpdA);
      h_YnEpdA->Fill(eventInfo.YnEpdA);
      h_XnEpdB->Fill(eventInfo.XnEpdB);
      h_YnEpdB->Fill(eventInfo.YnEpdB);

      h_psiTpc_RAW->Fill(eventInfo.psiTpc);
      h_psiTpcA_RAW->Fill(eventInfo.psiTpcA);
      h_psiTpcB_RAW->Fill(eventInfo.psiTpcB);
      h_psiEpd_RAW->Fill(eventInfo.psiEpd);
      h_psiEpdA_RAW->Fill(eventInfo.psiEpdA);
      h_psiEpdB_RAW->Fill(eventInfo.psiEpdB);



      //=========================================================
      //          Re-centering (Xn, Yn) Distributions
      //=========================================================

      if (RUN_ITERATION == 1 || RUN_ITERATION == 2)
	{
	  FlowUtils::recenterQ(eventInfo, correctionInputFile, ORDER_M);

	  if (eventInfo.badEvent) continue;

	  h_XnTpc_RC->Fill(eventInfo.XnTpc);
	  h_XnTpcA_RC->Fill(eventInfo.XnTpcA);
	  h_XnTpcB_RC->Fill(eventInfo.XnTpcB);
	  h_XnEpd_RC->Fill(eventInfo.XnEpd);
	  h_XnEpdA_RC->Fill(eventInfo.XnEpdA);
	  h_XnEpdB_RC->Fill(eventInfo.XnEpdB);

	  h_YnTpc_RC->Fill(eventInfo.YnTpc);
	  h_YnTpcA_RC->Fill(eventInfo.YnTpcA);
	  h_YnTpcB_RC->Fill(eventInfo.YnTpcB);
	  h_YnEpd_RC->Fill(eventInfo.YnEpd);
	  h_YnEpdA_RC->Fill(eventInfo.YnEpdA);
	  h_YnEpdB_RC->Fill(eventInfo.YnEpdB);

	  h_psiTpc_RC->Fill(eventInfo.psiTpc);
	  h_psiTpcA_RC->Fill(eventInfo.psiTpcA);
	  h_psiTpcB_RC->Fill(eventInfo.psiTpcB);
	  h_psiEpd_RC->Fill(eventInfo.psiEpd);
	  h_psiEpdA_RC->Fill(eventInfo.psiEpdA);
	  h_psiEpdB_RC->Fill(eventInfo.psiEpdB);

	  // Accumulate terms for averages over the re-centered angles for event plane angle shifting
	  for (int j = 1; j <= configs.shift_terms; j++)
	    {
	      p_sinAvgsTpc->Fill(j,  TMath::Sin((Double_t)j * ORDER_M * eventInfo.psiTpc));
	      p_cosAvgsTpc->Fill(j,  TMath::Cos((Double_t)j * ORDER_M * eventInfo.psiTpc));
	      p_sinAvgsTpcA->Fill(j, TMath::Sin((Double_t)j * ORDER_M * eventInfo.psiTpcA));
	      p_cosAvgsTpcA->Fill(j, TMath::Cos((Double_t)j * ORDER_M * eventInfo.psiTpcA));
	      p_sinAvgsTpcB->Fill(j, TMath::Sin((Double_t)j * ORDER_M * eventInfo.psiTpcB));
	      p_cosAvgsTpcB->Fill(j, TMath::Cos((Double_t)j * ORDER_M * eventInfo.psiTpcB));
	      p_sinAvgsEpd->Fill(j,  TMath::Sin((Double_t)j * ORDER_M * eventInfo.psiEpd));
	      p_cosAvgsEpd->Fill(j,  TMath::Cos((Double_t)j * ORDER_M * eventInfo.psiEpd));
	      p_sinAvgsEpdA->Fill(j, TMath::Sin((Double_t)j * ORDER_M * eventInfo.psiEpdA));
	      p_cosAvgsEpdA->Fill(j, TMath::Cos((Double_t)j * ORDER_M * eventInfo.psiEpdA));
	      p_sinAvgsEpdB->Fill(j, TMath::Sin((Double_t)j * ORDER_M * eventInfo.psiEpdB));
	      p_cosAvgsEpdB->Fill(j, TMath::Cos((Double_t)j * ORDER_M * eventInfo.psiEpdB));
	    }
	}
      //=========================================================
      //          End Re-centering
      //=========================================================



      //=========================================================
      //          Event Plane Angle Shifting and Flow
      //=========================================================

      if (RUN_ITERATION == 2)
	{
	  FlowUtils::shiftPsi(eventInfo, correctionInputFile, ORDER_M, configs.shift_terms);

	  h_psiTpc_FLAT->Fill(eventInfo.psiTpc);
	  h_psiTpcA_FLAT->Fill(eventInfo.psiTpcA);
	  h_psiTpcB_FLAT->Fill(eventInfo.psiTpcB);
	  h_psiEpd_FLAT->Fill(eventInfo.psiEpd);
	  h_psiEpdA_FLAT->Fill(eventInfo.psiEpdA);
	  h_psiEpdB_FLAT->Fill(eventInfo.psiEpdB);
	  //=========================================================
	  //          End Event Plane Angle Shifting
	  //=========================================================


	  // 2D Correlations between event planes
	  h2_psiEpdATpcA->Fill(eventInfo.psiTpcA,eventInfo.psiEpdA);
	  h2_psiEpdBTpcA->Fill(eventInfo.psiTpcA,eventInfo.psiEpdB);

	  h2_psiEpdATpcB->Fill(eventInfo.psiTpcB,eventInfo.psiEpdA);
	  h2_psiEpdBTpcB->Fill(eventInfo.psiTpcB,eventInfo.psiEpdB);

	  h2_psiTpcATpcB->Fill(eventInfo.psiTpcB,eventInfo.psiTpcA);

	  h2_psiEpdAEpdB->Fill(eventInfo.psiEpdB,eventInfo.psiEpdA);
	  //


	  // 1D correlation averages used in calculating resolution using the 3 sub-event method
	  p_TpcAB->Fill(eventInfo.centID,    TMath::Cos(ORDER_N * (eventInfo.psiTpcA - eventInfo.psiTpcB)));

	  p_TpcAEpdA->Fill(eventInfo.centID, TMath::Cos(ORDER_N * (eventInfo.psiTpcA - eventInfo.psiEpdA)));
	  p_TpcAEpdB->Fill(eventInfo.centID, TMath::Cos(ORDER_N * (eventInfo.psiTpcA - eventInfo.psiEpdB)));
	  p_TpcBEpdA->Fill(eventInfo.centID, TMath::Cos(ORDER_N * (eventInfo.psiTpcB - eventInfo.psiEpdA)));
	  p_TpcBEpdB->Fill(eventInfo.centID, TMath::Cos(ORDER_N * (eventInfo.psiTpcB - eventInfo.psiEpdB)));

	  p_EpdAEpdB->Fill(eventInfo.centID, TMath::Cos(ORDER_N * (eventInfo.psiEpdA - eventInfo.psiEpdB)));
	  //


	  //=========================================================
	  //          v_n Scan Plots
	  //=========================================================
	  // 2D searches through eta and centrality for correlations between detectors
	  // SHOULD PROBABLY ONLY USE THIS SECTION IF THE EPD REGIONS COVER THE WHOLE EPD!! Otherwise there might be gaps or undercounting in some bins.
	  Int_t tpcHits = eventInfo.tpcParticles.size();
	  Int_t epdHits = eventInfo.epdParticles.size();
	  Double_t phiTpc;
	  Double_t etaTpc;
	  Double_t phiEpd;
	  Double_t etaEpd;
	  Double_t psiTpc  = eventInfo.psiTpc;
	  Double_t psiEpd  = eventInfo.psiEpd;
	  Double_t psiEpdA = eventInfo.psiEpdA;
	  Double_t psiEpdB = eventInfo.psiEpdB;
	  Double_t psiTpcA = eventInfo.psiTpcA;
	  Double_t psiTpcB = eventInfo.psiTpcB;
	  Int_t centralityID = eventInfo.centID;

	  for (int j = 0; j < epdHits; j++)
	    {
	      phiEpd = eventInfo.epdParticles.at(j).phi;
	      etaEpd = eventInfo.epdParticles.at(j).eta;

	      h2_vnScanEpd->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpc)));
	      h2_vnScanEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcA)));
	      h2_vnScanEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcB)));
	    }
	  for (int j = 0; j < tpcHits; j++)
	    {
	      phiTpc = eventInfo.tpcParticles.at(j).phi;
	      etaTpc = eventInfo.tpcParticles.at(j).eta;

	      h2_vnScanTpc->Fill(etaTpc, centralityID, TMath::Cos(ORDER_M * (phiTpc - psiEpd)));
	      h2_vnScanTpcEpdA->Fill(etaTpc, centralityID, TMath::Cos(ORDER_M * (phiTpc - psiEpdA)));
	      h2_vnScanTpcEpdB->Fill(etaTpc, centralityID, TMath::Cos(ORDER_M * (phiTpc - psiEpdB)));
	    }
	  //=========================================================
	  //          End v_n Scan Plots
	  //=========================================================


	  //=========================================================
	  //        Flow Calculations
	  //=========================================================
	  Double_t jthPhi;
	  //Double_t jthEta;
	  Double_t jthpT;
	  Double_t jthKT;
	  Double_t jthRapidity;
	  Double_t psi = eventInfo.psiEpdA;
	  //Double_t psi_epdA = eventInfo.psiEpdA;
	  //Double_t psi_tpcB = eventInfo.psiTpcB;
	  Double_t tpcEfficiency = 1.0;  // Default
	  Int_t centID = eventInfo.centID;

	  if (configs.sqrt_s_NN == 3.0 && centID < 4) continue;  // ONLY LOOKING AT CENTRALITY 60% AND LOWER FOR 3.0 GeV
	  else if (configs.sqrt_s_NN == 3.22 && centID < 4) continue;


	  // JUST v1 FOR WEIGHTING (NOT READY TO USE)
	  /*
	  // TPC B region
	  for (UInt_t j = 0; j < eventInfo.tpcParticles.size(); j++)
	    {
	      tpcEfficiency = 1.0; // Make sure efficiency is reset each time.

	      jthPhi = eventInfo.tpcParticles.at(j).phi;
	      jthEta = eventInfo.tpcParticles.at(j).eta;
	      jthpT  = eventInfo.tpcParticles.at(j).pT;
	      jthKT  = eventInfo.tpcParticles.at(j).KT;
	      jthRapidity = eventInfo.tpcParticles.at(j).rapidity;
	      if (jthPhi == FlowUtils::D_BAD_VALUE || jthpT == FlowUtils::D_BAD_VALUE || jthEta == FlowUtils::D_BAD_VALUE ||
		  jthKT == FlowUtils::D_BAD_VALUE  || jthRapidity == FlowUtils::D_BAD_VALUE) 
		continue;

	      if (efficienciesFound && configs.sqrt_s_NN == 3.0)
		{
		  if      (eventInfo.tpcParticles.at(j).prTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_pr);
		  else if (eventInfo.tpcParticles.at(j).deTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_de);
		  else if (eventInfo.tpcParticles.at(j).trTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_tr);
		}
	      if (tpcEfficiency == -1) { continue; }
		  
	      if (eventInfo.tpcParticles.at(j).isInTpcB)
		{ 
		  //p2_v1_eta_cent_TPCB->Fill(centID, jthEta, TMath::Cos(1.0 * (jthPhi - psi_epdA))); 

		  if (eventInfo.tpcParticles.at(j).prTag)
		    { p2_v1_pT_eta_TPCB_pr->Fill(jthEta, jthpT, TMath::Cos(1.0 * (jthPhi - psi_epdA)) / tpcEfficiency); }
		}
	    }// END TPC v1

	  // EPD regions

	  for (UInt_t j = 0; j < eventInfo.epdParticles.size(); j++)
	    {
	      jthPhi = eventInfo.epdParticles.at(j).phi;
	      jthEta = eventInfo.epdParticles.at(j).eta;
	      if (jthPhi == FlowUtils::D_BAD_VALUE || jthEta == FlowUtils::D_BAD_VALUE)
		continue;

	      if (eventInfo.epdParticles.at(j).isInEpdA)
		{ p2_v1_eta_cent_EPDA->Fill(centID, jthEta, TMath::Cos(1.0 * (jthPhi - psi_tpcB))); }
	      else if (eventInfo.epdParticles.at(j).isInEpdB)
		{ p2_v1_eta_cent_EPDB->Fill(centID, jthEta, TMath::Cos(1.0 * (jthPhi - psi_tpcB))); }	    
	    }// End EPD v1
	  */

	  // "OBSERVED" FLOW VALUES HERE, no resolutions or efficiencies
	  for (UInt_t j = 0; j < eventInfo.tpcParticles.size(); j++)
	    {
	      jthPhi = eventInfo.tpcParticles.at(j).phi;
	      jthpT  = eventInfo.tpcParticles.at(j).pT;
	      jthKT  = eventInfo.tpcParticles.at(j).KT;
	      jthRapidity = eventInfo.tpcParticles.at(j).rapidity;
	      if (jthPhi == FlowUtils::D_BAD_VALUE || jthpT == FlowUtils::D_BAD_VALUE || 
		  jthKT == FlowUtils::D_BAD_VALUE  || jthRapidity == FlowUtils::D_BAD_VALUE) 
		continue;

	      // PI+
	      // Normal acceptance 0 < yCM < 0.5
	      if (eventInfo.tpcParticles.at(j).ppTag && 
		  jthpT >= configs.pt_norm_pi_low && jthpT <= configs.pt_norm_pi_high &&
		  jthRapidity - Y_MID > configs.yCM_norm_pi_low && jthRapidity - Y_MID < configs.yCM_norm_pi_high)
		{ p_vn_pp_obs->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi))); }
	      // PI-
	      // Normal acceptance 0 < yCM < 0.5
	      else if (eventInfo.tpcParticles.at(j).pmTag && 
		       jthpT >= configs.pt_norm_pi_low && jthpT <= configs.pt_norm_pi_high &&
		       jthRapidity - Y_MID > configs.yCM_norm_pi_low && jthRapidity - Y_MID < configs.yCM_norm_pi_high)
		{ p_vn_pm_obs->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi))); }
	      // K+
	      // Normal acceptance 0 < yCM < 0.5
	      else if (eventInfo.tpcParticles.at(j).kpTag && 
		       jthpT >= configs.pt_norm_ka_low && jthpT <= configs.pt_norm_ka_high &&
		       jthRapidity - Y_MID > configs.yCM_norm_ka_low && jthRapidity - Y_MID < configs.yCM_norm_ka_high)
		{ p_vn_kp_obs->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi))); }
	      // K-
	      // Normal acceptance 0 < yCM < 0.5
	      else if (eventInfo.tpcParticles.at(j).kmTag && 
		       jthpT >= configs.pt_norm_ka_low && jthpT <= configs.pt_norm_ka_high &&
		       jthRapidity - Y_MID > configs.yCM_norm_ka_low && jthRapidity - Y_MID < configs.yCM_norm_ka_high)
		{ p_vn_km_obs->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi))); }
	      // PROTON
	      else if (eventInfo.tpcParticles.at(j).prTag)
		{
		  // NORMAL ACCEPTANCE 0 < y_cm < 0.5
		  if (jthRapidity - Y_MID > configs.yCM_norm_pr_low && jthRapidity - Y_MID < configs.yCM_norm_pr_high &&
		      jthpT > configs.pt_norm_pr_low && jthpT < configs.pt_norm_pr_high)
		    { p_vn_pr_obs->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi))); }
		  // ALTERNATE ACCEPTANCE REGION
		  if (jthKT/1.0 >= configs.KT_pdt_low && jthKT/1.0 <= configs.KT_pdt_high &&
		      jthRapidity - Y_MID > configs.yCM_alt_pr_low && jthRapidity - Y_MID < configs.yCM_alt_pr_high)
		    { p_vn_pr_alt_obs->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi))); }
		}
	      // DEUTERON
	      else if (eventInfo.tpcParticles.at(j).deTag && 
		       jthKT/2.0 >= configs.KT_pdt_low && jthKT/2.0 <= configs.KT_pdt_high &&
		       jthRapidity - Y_MID > configs.yCM_norm_de_low && jthRapidity - Y_MID < configs.yCM_norm_de_high)
		{ p_vn_de_obs->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi))); }

	      // TRITON
	      else if (eventInfo.tpcParticles.at(j).trTag && 
		       jthKT/3.0 >= configs.KT_pdt_low && jthKT/3.0 <= configs.KT_pdt_high &&
		       jthRapidity - Y_MID > configs.yCM_norm_tr_low && jthRapidity - Y_MID < configs.yCM_norm_tr_high)
		{ p_vn_tr_obs->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi))); }
	    }
	  //////


	  // RESOLUTION CORRECTED FLOW VALUES HERE
	  if (resolutionsFound)
	    {
	      TH1D *resolutionHistogram = (TH1D*)resolutionInputFile->Get("h_resolutions");
	      Double_t resolution = resolutionHistogram->GetBinContent(centID+1);	      
	      if (resolution == 0.0) continue;  // Skip centralities without resolutions.


	      for (UInt_t j = 0; j < eventInfo.tpcParticles.size(); j++)
		{
		  tpcEfficiency = 1.0; // Make sure efficiency is reset each time.

		  jthPhi = eventInfo.tpcParticles.at(j).phi;
		  jthpT  = eventInfo.tpcParticles.at(j).pT;
		  jthKT  = eventInfo.tpcParticles.at(j).KT;
		  jthRapidity = eventInfo.tpcParticles.at(j).rapidity;
		  if (jthPhi == FlowUtils::D_BAD_VALUE || jthpT == FlowUtils::D_BAD_VALUE || 
		      jthKT == FlowUtils::D_BAD_VALUE  || jthRapidity == FlowUtils::D_BAD_VALUE) 
		    continue;
		  
		  h_simulationCheck_total->Fill(1);
		  if (efficienciesFound && configs.sqrt_s_NN == 3.0)
		    {
		      if      (eventInfo.tpcParticles.at(j).ppTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_pp);
		      else if (eventInfo.tpcParticles.at(j).pmTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_pm);
		      else if (eventInfo.tpcParticles.at(j).kpTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_kp);
		      else if (eventInfo.tpcParticles.at(j).kmTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_km);
		      else if (eventInfo.tpcParticles.at(j).prTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_pr);
		      else if (eventInfo.tpcParticles.at(j).deTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_de);
		      else if (eventInfo.tpcParticles.at(j).trTag) tpcEfficiency = FlowUtils::getTpcEff(jthRapidity - Y_MID, jthpT, h2_tracking_tr);
		    }
		  if (tpcEfficiency == -1) // Checks here for tracks with no recorded efficiency values.
		    { 
		      if (eventInfo.tpcParticles.at(j).ppTag)
			{
			  h_simulationCheck_pp->Fill(1);
			  h2_pT_vs_yCM_pp_noEff->Fill(jthRapidity - Y_MID, jthpT);
			}
		      else if (eventInfo.tpcParticles.at(j).pmTag)
			{
			  h_simulationCheck_pm->Fill(1);
			  h2_pT_vs_yCM_pm_noEff->Fill(jthRapidity - Y_MID, jthpT);
			}
		      else if (eventInfo.tpcParticles.at(j).kpTag)
			{
			  h_simulationCheck_kp->Fill(1);
			  h2_pT_vs_yCM_kp_noEff->Fill(jthRapidity - Y_MID, jthpT);
			}
		      else if (eventInfo.tpcParticles.at(j).kmTag)
			{
			  h_simulationCheck_km->Fill(1);
			  h2_pT_vs_yCM_km_noEff->Fill(jthRapidity - Y_MID, jthpT);
			}
		      else if (eventInfo.tpcParticles.at(j).prTag)
			{
			  h_simulationCheck_pr->Fill(1);
			  h2_pT_vs_yCM_pr_noEff->Fill(jthRapidity - Y_MID, jthpT);
			}
		      else if (eventInfo.tpcParticles.at(j).deTag) 
			{
			  h_simulationCheck_de->Fill(1); 
			  h2_pT_vs_yCM_de_noEff->Fill(jthRapidity - Y_MID, jthpT);
			}
		      else if (eventInfo.tpcParticles.at(j).trTag) 
			{
			  h_simulationCheck_tr->Fill(1); 
			  h2_pT_vs_yCM_tr_noEff->Fill(jthRapidity - Y_MID, jthpT);
			}
		      
		      continue; 
		    }


		  // PI+
		  if (eventInfo.tpcParticles.at(j).ppTag)
		    {
		      // Normal acceptance 0 < yCM < 0.5
		      if (jthpT >= configs.pt_norm_pi_low && jthpT <= configs.pt_norm_pi_high &&
			  jthRapidity - Y_MID > configs.yCM_norm_pi_low && jthRapidity - Y_MID < configs.yCM_norm_pi_high)
			{
			  p2_vn_yCM_cent_pp->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p_vn_pp->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			  p2_vn_pT_cent_pp->Fill(centID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_KT_cent_pp->Fill(centID, jthKT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		      // Extended rapidity acceptance 0.5 <= yCM < 1.0
		      else if (jthpT >= configs.pt_yExt_pi_low && jthpT <= configs.pt_yExt_pi_high &&
			       jthRapidity - Y_MID >= configs.yCM_yExt_pi_low && jthRapidity - Y_MID < configs.yCM_yExt_pi_high)
			{
			  p2_vn_yCM_cent_pp->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p_vn_pp_ext->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		      // rapidity stratified by pT
		      if (jthRapidity - Y_MID > -1.0 && jthRapidity - Y_MID < 1.0 && 
			  jthpT > configs.pt_norm_pi_low && jthpT < configs.pt_norm_pi_high)
			{ p2_vn_pT_vs_yCM_pp->Fill(jthRapidity - Y_MID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }
		    }
		  
		  // PI-
		  else if (eventInfo.tpcParticles.at(j).pmTag)
		    {
		      // Normal acceptance 0 < yCM < 0.5
		      if (jthpT >= configs.pt_norm_pi_low && jthpT <= configs.pt_norm_pi_high &&
			  jthRapidity - Y_MID > configs.yCM_norm_pi_low && jthRapidity - Y_MID < configs.yCM_norm_pi_high)
			{
			  p2_vn_yCM_cent_pm->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p_vn_pm->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			  p2_vn_pT_cent_pm->Fill(centID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_KT_cent_pm->Fill(centID, jthKT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		      // Extended rapidity acceptance 0.5 <= yCM < 1.0
		      else if (jthpT >= configs.pt_yExt_pi_low && jthpT <= configs.pt_yExt_pi_high &&
			       jthRapidity - Y_MID >= configs.yCM_yExt_pi_low && jthRapidity - Y_MID < configs.yCM_yExt_pi_high)
			{
			  p2_vn_yCM_cent_pm->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p_vn_pm_ext->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		      // rapidity stratified by pT
		      if (jthRapidity - Y_MID > -1.0 && jthRapidity - Y_MID < 1.0 && 
			  jthpT > configs.pt_norm_pi_low && jthpT < configs.pt_norm_pi_high)
			{ p2_vn_pT_vs_yCM_pm->Fill(jthRapidity - Y_MID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }
		    }

		  // K+
		  else if (eventInfo.tpcParticles.at(j).kpTag)
		    {
		      // Normal acceptance 0 < yCM < 0.5
		      if (jthpT >= configs.pt_norm_ka_low && jthpT <= configs.pt_norm_ka_high &&
			  jthRapidity - Y_MID > configs.yCM_norm_ka_low && jthRapidity - Y_MID < configs.yCM_norm_ka_high)
			{
			  p2_vn_yCM_cent_kp->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p_vn_kp->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			  p2_vn_pT_cent_kp->Fill(centID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_KT_cent_kp->Fill(centID, jthKT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		      // Extended rapidity acceptance 0.5 <= yCM < 1.0
		      else if (jthpT >= configs.pt_yExt_ka_low && jthpT <= configs.pt_yExt_ka_high &&
			       jthRapidity - Y_MID >= configs.yCM_yExt_ka_low && jthRapidity - Y_MID < configs.yCM_yExt_ka_high)
			{
			  p2_vn_yCM_cent_kp->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p_vn_kp_ext->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		    }

		  // K-
		  else if (eventInfo.tpcParticles.at(j).kmTag)
		    {
		      // Normal acceptance 0 < yCM < 0.5
		      if (jthpT >= configs.pt_norm_ka_low && jthpT <= configs.pt_norm_ka_high &&
			  jthRapidity - Y_MID > configs.yCM_norm_ka_low && jthRapidity - Y_MID < configs.yCM_norm_ka_high)
			{
			  p2_vn_yCM_cent_km->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p_vn_km->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			  p2_vn_pT_cent_km->Fill(centID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_KT_cent_km->Fill(centID, jthKT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		      // Extended rapidity acceptance 0.5 <= yCM < 1.0
		      else if (jthpT >= configs.pt_yExt_ka_low && jthpT <= configs.pt_yExt_ka_high &&
			       jthRapidity - Y_MID >= configs.yCM_yExt_ka_low && jthRapidity - Y_MID < configs.yCM_yExt_ka_high)
			{
			  p2_vn_yCM_cent_km->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p_vn_km_ext->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		    }

		  // PROTON
		  else if (eventInfo.tpcParticles.at(j).prTag)
		    {
		      if (jthRapidity - Y_MID > -1.0 && jthRapidity - Y_MID < 1.0 && 
			  jthpT > 0.4 && jthpT < 2.5)
			{ p2_vn_pT_vs_yCM_pr->Fill(jthRapidity - Y_MID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }
		      
		      // RAPIDITY DEPENDENT PLOT
		      if (jthRapidity - Y_MID > configs.yCM_yDep_pr_low && jthRapidity - Y_MID < configs.yCM_yDep_pr_high && 
			  jthpT > configs.pt_yDep_pr_low && jthpT < configs.pt_yDep_pr_high)
			{ p2_vn_yCM_cent_pr->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }

		      // NORMAL ACCEPTANCE 0 < y_cm < 0.5
		      if (jthRapidity - Y_MID > configs.yCM_norm_pr_low && jthRapidity - Y_MID < configs.yCM_norm_pr_high &&
			  jthpT > configs.pt_norm_pr_low && jthpT < configs.pt_norm_pr_high)
			{ 
			  p_vn_pr->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			  p2_vn_pT_cent_pr->Fill(centID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_KT_cent_pr->Fill(centID, jthKT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));

			  if (jthpT < 1.0)
			    { p_vn_pr_pTlt1->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }
			  else if (jthpT > 1.0)
			    { p_vn_pr_pTgt1->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }
			}
		      // EXTENDED RAPIDITY 0.5 <= y_cm < 1.0
		      else if (jthRapidity - Y_MID >= configs.yCM_yExt_pr_low && jthRapidity - Y_MID < configs.yCM_yExt_pr_high &&
			       jthpT > configs.pt_yExt_pr_low && jthpT < configs.pt_yExt_pr_high)
			{ p_vn_pr_ext->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }

		      // ALTERNATE ACCEPTANCE REGION
		      if (jthKT/1.0 >= configs.KT_pdt_low && jthKT/1.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > configs.yCM_alt_pr_low && jthRapidity - Y_MID < configs.yCM_alt_pr_high)
			{ 
			  p_vn_pr_alt->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_yCM_cent_pr_alt->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_pT_cent_pr_alt->Fill(centID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_KT_cent_pr_alt->Fill(centID, jthKT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  h2_dndDeltaPhi_vs_cent_pr_alt->Fill(centID, jthPhi - psi);
			}
		      if (jthKT/1.0 >= configs.KT_pdt_low && jthKT/1.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > 0.0 && jthRapidity - Y_MID < 0.6)
			{
			  p_vn_pr_alt_y0to0p6->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			}

		      if (jthKT/1.0 >= configs.KT_pdt_low && jthKT/1.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > 0.6 && jthRapidity - Y_MID < 1.0)
			{
			  p_vn_pr_alt_y0p6to1p0->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			}
		      if (jthKT/1.0 >= configs.KT_pdt_low && jthKT/1.0 <= configs.KT_pdt_high &&
			  jthRapidity/Y_BEAM > -1.0 && jthRapidity/Y_BEAM < 0.0)
			{ 
			  p2_vn_yOverYbeam_cent_pr_alt->Fill(centID, jthRapidity/Y_BEAM, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		      // RAPIDITY SYMMETRIC ACCEPTANCE REGION
		      if (jthRapidity - Y_MID > configs.yCM_ySym_pr_low && jthRapidity - Y_MID < configs.yCM_ySym_pr_high && 
			  jthpT > configs.pt_ySym_pr_low && jthpT < configs.pt_ySym_pr_high)
			{ p2_vn_yCM_cent_pr_symmetry->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }

		      // ONLY FORWARD ACCEPTANCE REGION
		      if (jthRapidity - Y_MID > configs.yCM_yFor_pr_low && jthRapidity - Y_MID < configs.yCM_yFor_pr_high && 
			  jthpT > configs.pt_yFor_pr_low && jthpT < configs.pt_yFor_pr_high)
			{ p_vn_pr_for->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); }
		    }

		  // DEUTERON
		  else if (eventInfo.tpcParticles.at(j).deTag)
		    {
		      if (jthKT/2.0 >= configs.KT_pdt_low && jthKT/2.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > configs.yCM_norm_de_low && jthRapidity - Y_MID < configs.yCM_norm_de_high)
			{
			  p2_vn_yCM_cent_de->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));

			  p_vn_de->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			  p2_vn_pT_cent_de->Fill(centID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_KT_cent_de->Fill(centID, jthKT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  h2_dndDeltaPhi_vs_cent_de->Fill(centID, jthPhi - psi);
			}

		      if (jthKT/2.0 >= configs.KT_pdt_low && jthKT/2.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > 0.0 && jthRapidity - Y_MID < 0.6)
			{
			  p_vn_de_y0to0p6->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			}

		      if (jthKT/2.0 >= configs.KT_pdt_low && jthKT/2.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > 0.6 && jthRapidity - Y_MID < 1.0)
			{
			  p_vn_de_y0p6to1p0->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			}

		      if (jthKT/2.0 >= configs.KT_pdt_low && jthKT/2.0 <= configs.KT_pdt_high &&
			  jthRapidity/Y_BEAM > -1.0 && jthRapidity/Y_BEAM < 0.0)
			{
			  p2_vn_yOverYbeam_cent_de->Fill(centID, jthRapidity/Y_BEAM, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}
		    }

		  // TRITON
		  else if (eventInfo.tpcParticles.at(j).trTag)
		    {
		      if (jthKT/3.0 >= configs.KT_pdt_low && jthKT/3.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > configs.yCM_norm_tr_low && jthRapidity - Y_MID < configs.yCM_norm_tr_high)
			{
			  p2_vn_yCM_cent_tr->Fill(centID, jthRapidity - Y_MID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));

			  p_vn_tr->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			  p2_vn_pT_cent_tr->Fill(centID, jthpT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  p2_vn_KT_cent_tr->Fill(centID, jthKT, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			  h2_dndDeltaPhi_vs_cent_tr->Fill(centID, jthPhi - psi);
			}

		      if (jthKT/3.0 >= configs.KT_pdt_low && jthKT/3.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > 0.0 && jthRapidity - Y_MID < 0.6)
			{
			  p_vn_tr_y0to0p6->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			}

		      if (jthKT/3.0 >= configs.KT_pdt_low && jthKT/3.0 <= configs.KT_pdt_high &&
			  jthRapidity - Y_MID > 0.6 && jthRapidity - Y_MID < 1.0)
			{
			  p_vn_tr_y0p6to1p0->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency)); 
			}

		      if (jthKT/3.0 >= configs.KT_pdt_low && jthKT/3.0 <= configs.KT_pdt_high &&
			  jthRapidity/Y_BEAM > -1.0 && jthRapidity/Y_BEAM < 0.0)
			{
			  p2_vn_yOverYbeam_cent_tr->Fill(centID, jthRapidity/Y_BEAM, TMath::Cos(ORDER_N * (jthPhi - psi)) / (resolution * tpcEfficiency));
			}

		    }
		}// End tpc particles loop
	    }// End if(resolutionsFound)
	  //=========================================================
	  //            End Flow Calculations
	  //=========================================================
	}// End if(RUN_ITERATION == 2)
    }//END EVENT LOOP
  eventInfo.reset();

  // Switch axis labels on some plots, put in centrality percentages

  // The indices here are equivalent to the corresponding centrality ID
  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};

  Int_t labelIndex;
  for (int i = 1; i <= CENT_BINS; i++) 
    {
      labelIndex = FIRST_CENT + i - 1;
      h2_vnScanTpc->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]); 
      h2_vnScanTpcEpdA->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]); 
      h2_vnScanTpcEpdB->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]); 
      h2_vnScanEpd->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_vnScanEpdTpcA->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_vnScanEpdTpcB->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
    }
  ////

  // Manually write the few plots that were pulled from the trees
  h_eventCheck->Write();
  h_trackmult->Write();
  h_refmult->Write();
  h_tofmult->Write();
  h2_refmult_vs_trackmult->Write();
  h2_tofmult_vs_trackmult->Write();
  h2_tofmult_vs_refmult->Write();
  ////

  // Check ending memory usage to help spot memory leaks
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;  
  ret = getrusage(who, &usage);

  if (ret == 0) { std::cout << "Ending memory usage: " << usage.ru_maxrss / 1000 << " MB" << std::endl; }
  else { std::cout << "Could not retrieve memory usage!" << std::endl; }
  //

  // Write all main output
  outputFile->cd();
  outputFile->Write();
  //

  // Save re-centering and Fourier shifting information
  if (RUN_ITERATION == 0 || RUN_ITERATION == 1)
    {
      correctionOutputFile->cd();

      p_sinAvgsTpc   ->Write();
      p_cosAvgsTpc   ->Write();
      p_sinAvgsTpcA  ->Write();
      p_cosAvgsTpcA  ->Write();
      p_sinAvgsTpcB  ->Write();
      p_cosAvgsTpcB  ->Write();
      p_sinAvgsEpd   ->Write();
      p_cosAvgsEpd   ->Write();
      p_sinAvgsEpdA  ->Write();
      p_cosAvgsEpdA  ->Write();
      p_sinAvgsEpdB  ->Write();
      p_cosAvgsEpdB  ->Write();
      h_XnTpc        ->Write();
      h_YnTpc        ->Write();
      h_XnTpcA       ->Write();
      h_YnTpcA       ->Write();
      h_XnTpcB       ->Write();
      h_YnTpcB       ->Write();
      h_XnEpd        ->Write();
      h_YnEpd        ->Write();
      h_XnEpdA       ->Write();
      h_YnEpdA       ->Write();
      h_XnEpdB       ->Write();
      h_YnEpdB       ->Write();

      gROOT->GetListOfFiles()->Remove(correctionOutputFile);
      correctionOutputFile->Close();
    }
  //

  // Close files
  gROOT->GetListOfFiles()->Remove(outputFile);
  outputFile->Close();
  if (RUN_ITERATION == 1 || RUN_ITERATION == 2)
    {
      gROOT->GetListOfFiles()->Remove(correctionInputFile);
      correctionInputFile->Close();
    }
  //

  std::cout << "Done!" << std::endl;
}//End main()
