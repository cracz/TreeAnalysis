// C++ headers
#include <iostream>
#include <vector>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TChain.h"
#include "TSystem.h"
#include "TKey.h"
#include "TVector3.h"

// Configuration file reader
#include "StRoot/ConfigReader/ConfigReader.h"

// My Util Header
#include "../include/FlowUtils.h"


//const Double_t D_M0_PI = 0.139571;   //Rest masses
//const Double_t D_M0_KA = 0.493677;
const Double_t D_M0_PR = 0.938272;
//const Double_t D_M0_DE = 1.875613;   // Deuteron
//const Double_t D_M0_TR = 2.808921;   // Triton
//const Double_t D_M0_HE3 = 2.809414;   // Helium-3
//const Double_t D_M0_AL = 3.727379;   // Alpha

//void nSigmaInfo(TString inFile, TString jobID, std::string configFileName)
int main(int argc, char *argv[])
{
  std::cout << "Initializing..." << std::endl;

  TString inFile = argv[1];
  TString jobID  = argv[2];
  std::string configFileName = argv[3];

  if (gSystem->AccessPathName(inFile)) { std::cout << "Error reading input file!" << std::endl; return 1;}

  //=========================================================
  //          Set up various files
  //=========================================================
  ConfigReader configs;
  configs.read(configFileName);
  if (configs.errorFound()) { std::cout << "There was an error reading the configurations! Aborting analysis!" << std::endl; return 1; }
  configs.checkForNonSetKeys();

  //=== INITIALIZE TTREE
  Int_t N_track = 0;  // Max number of tracks in an event. Depends on energy and centrality definition!
  if      (configs.sqrt_s_NN == 3.0)  { N_track = 195;  }
  else if (configs.sqrt_s_NN == 3.2)  { N_track = 287;  }
  else if (configs.sqrt_s_NN == 3.5)  { N_track = 325;  }
  else if (configs.sqrt_s_NN == 3.9)  { N_track = 344;  }
  else if (configs.sqrt_s_NN == 4.5)  { N_track = 367;  }  // UPDATE THIS WHEN CENTRALITY IS OFFICIAL!
  else if (configs.sqrt_s_NN == 7.2)  { N_track = 240;  }
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
  Float_t dEdxError[N_track];
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
  tree->SetBranchAddress("dEdxError",&dEdxError);
  tree->SetBranchAddress("nHits",&nHits);
  tree->SetBranchAddress("nHitsFit",&nHitsFit);
  tree->SetBranchAddress("nHitsPoss",&nHitsPoss);
  tree->SetBranchAddress("nHitsDedx",&nHitsDedx);
  tree->SetBranchAddress("nEPDhits",&i_nEPDhits);
  tree->SetBranchAddress("EPDid",&EPDids);
  tree->SetBranchAddress("EPDnMip",&EPDnMip);
  //=== END TTREE SETUP


  // MAIN OUTPUT FILE
  TString outFile = jobID+".root";
  TFile *outputFile = new TFile(outFile,"RECREATE");
  outputFile->cd();
  ////

  //=========================================================
  //          END file setup
  //=========================================================

  //TH3D *h3_lndEdx_qp_y = new TH3D("h3_lndEdx_qp_y", ";q*|p|;ln(dE/dx)", 2000, -10, 10, 500, 0, 4.5, 40, 0, 4);
  TH2D* h2;
  std::vector<TH2D*> vectorOfHistos;
  Double_t lowerRapidityBound = -0.1;
  Double_t upperRapidityBound = 0.0;
  TString title;
  TString name;

  for (int ithRapidityBin = 0; ithRapidityBin < 40; ithRapidityBin++)
    {
      name.Form("hLndEdxVsMomVsRap_nhitsdedxcut_%d", ithRapidityBin);
      title.Form("Rapidity %.1f - %.1f;q*|p|;ln(dE/dx)", lowerRapidityBound, upperRapidityBound);
      h2 = new TH2D(name, title, 2000, -10, 10, 500, 0, 4.5);
      vectorOfHistos.push_back(h2);
      lowerRapidityBound -= 0.1;
      upperRapidityBound -= 0.1;
    }

  //TH2D *hLndEdxVsMomVsRap_nhitsdedxcut_0 = new TH2D("hLndEdxVsMomVsRap_nhitsdedxcut_0", "Rapidity 0-0.1;q*|p|;ln(dE/dx)", 2000, -10, 10, 500, 0, 4.5);
  

  // EVENT LOOP
  Int_t events2read = tree->GetEntries();
  std::cout << "Setup complete, beginning analysis on " << events2read << " events..." << std::endl;
  for (Long64_t ievent = 0; ievent < events2read; ievent++)
    {
      tree->GetEntry(ievent);

      // At this point, all bad runs and bad trigger events are removed by TreeMaker.
      // Now check event vertex
      TVector3 pVtx(f_xvtx, f_yvtx, f_zvtx);
      Double_t d_rvtx = TMath::Sqrt(f_xvtx * f_xvtx + (f_yvtx + 2) * (f_yvtx + 2));
      
      Bool_t goodVertexZ = (f_zvtx > configs.z_vtx_low && f_zvtx < configs.z_vtx_high);
      if (!goodVertexZ) continue;

      Bool_t goodVertexR = (d_rvtx < 1.5);
      if(!goodVertexR) continue;

      if (i_centrality == -99) continue;  // Remove undefined events

      Double_t d_px;
      Double_t d_py;
      Double_t d_pz;
      Double_t d_mom;
      Double_t d_DCA;
      Double_t d_dEdx;
      Double_t d_rapidity_assumingProton;
      Int_t i_nHits;
      Int_t i_nHitsFit;
      Int_t i_nHitsPoss;
      Int_t i_nHitsDedx;
      Short_t s_charge;

      // TRACK LOOP OVER PRIMARY TRACKS
      Int_t nTracks = (Int_t)i_trackNumber;
      for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
	{
	  d_px  = Px[iTrk];
	  d_py  = Py[iTrk];
	  d_pz  = Pz[iTrk];
	  d_mom = FlowUtils::totalMomentum(Px[iTrk], Py[iTrk], Pz[iTrk]);
	  s_charge = charge[iTrk];
	  d_DCA = DCA[iTrk];
	  d_dEdx = dEdx[iTrk];
	  d_rapidity_assumingProton = FlowUtils::rapidity(d_px, d_py, d_pz, D_M0_PR);
	  i_nHits = nHits[iTrk];
	  i_nHitsFit = nHitsFit[iTrk];
	  i_nHitsPoss = nHitsPoss[iTrk];
	  i_nHitsDedx = nHitsDedx[iTrk];

	  //=========================================================
	  //          Track QA Cuts
	  //=========================================================

	  bool b_bad_hits  = ( i_nHits < configs.nHits );
	  bool b_bad_dEdx  = ( i_nHitsDedx < 20 );
	  bool b_bad_ratio = ( ((double)i_nHitsFit / (double)i_nHitsPoss) <= configs.nHits_ratio );
	  bool b_bad_DCA   = ( d_DCA >= configs.dca );

	  if (b_bad_hits || b_bad_dEdx || b_bad_ratio || b_bad_DCA) continue;

	  //bool b_bad_dEdx  = ( i_nHitsDedx < 20 );
	  //if (b_bad_dEdx) continue;
	  //=========================================================
	  //          End Track QA Cuts
	  //=========================================================

	  // Sort this one track into the proper histogram based on d_rapidity_assumingProton
	  lowerRapidityBound = -0.1;
	  upperRapidityBound = 0.0;

	  for (int jthRapidityBin = 0; jthRapidityBin < 40; jthRapidityBin++)
	    {
	      if (d_rapidity_assumingProton > lowerRapidityBound && d_rapidity_assumingProton < upperRapidityBound)
		{
		  vectorOfHistos.at(jthRapidityBin)->Fill(s_charge * d_mom, TMath::Log(d_dEdx));
		  break;
		}
	      else
		{
		  lowerRapidityBound -= 0.1;
		  upperRapidityBound -= 0.1;
		}
	    }
	}//End track loop
    }// End event loop


  for (int kthRapidityBin = 0; kthRapidityBin < 40; kthRapidityBin++)
    {
      vectorOfHistos.at(kthRapidityBin)->Write();
      delete vectorOfHistos.at(kthRapidityBin);
    }
  
  gROOT->GetListOfFiles()->Remove(outputFile);
  outputFile->Close();

  std::cout << "Done!" << std::endl;
}//End nSigmaInfo()
