// This is an version of Yang Wu's (UCR) tree maker adapted for my own (Cameron Racz) analysis.


// Load header files
#include "TreeMaker.h"

// Config header
#include "../ConfigReader/ConfigReader.h"

// Class implementation in CINT
ClassImp(TreeMaker)

// Constructor
TreeMaker::TreeMaker(StPicoDstMaker* Maker, std::string configFileName, TString JobId, Int_t EventsNumber, Double_t inputParameter) : StMaker()
{
  // Initialize parameters
  mPicoDstMaker = Maker;
  configs.read(configFileName);
  JobIdName = JobId;
  cutTest = inputParameter;
  Y_MID   = configs.y_mid;
  JobIdName.Append(".root"); // Name output file by assigned Job ID

  pileupBoundHigh_3p2GeV = new TF1("pileupBoundHigh_3p2GeV", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x", 0, 100);
  pileupBoundHigh_3p2GeV->SetParameter(0, 19.48);
  pileupBoundHigh_3p2GeV->SetParameter(1, 5.428);
  pileupBoundHigh_3p2GeV->SetParameter(2, -0.007);
  pileupBoundHigh_3p2GeV->SetParameter(3, -2.428e-4);
  pileupBoundHigh_3p2GeV->SetParameter(4, 1.197e-7);

  pileupBoundLow_3p2GeV = new TF1("pileupBoundLow_3p2GeV", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x", 0, 100);
  pileupBoundLow_3p2GeV->SetParameter(0, -13.59);
  pileupBoundLow_3p2GeV->SetParameter(1, 1.515);
  pileupBoundLow_3p2GeV->SetParameter(2, 0.02816);
  pileupBoundLow_3p2GeV->SetParameter(3, -1.195e-4);
  pileupBoundLow_3p2GeV->SetParameter(4, -9.639e-7);
}

// Destructor
TreeMaker::~TreeMaker() {}

Int_t TreeMaker::Init()
{
  if (configs.errorFound()) return kStErr;

  outputFile = new TFile(JobIdName,"recreate");
  
  // New Tree
  FxtTree = new TTree("Autree","TTree to hold FXT events and tracks");
  FxtTree->Branch("runId",&tree_runId,"runId/I");
  FxtTree->Branch("eventId",&tree_eventId,"eventId/I");
  FxtTree->Branch("bField",&tree_bField,"bField/F");
  FxtTree->Branch("Vx",&tree_Vx,"Vx/F");
  FxtTree->Branch("Vy",&tree_Vy,"Vy/F");
  FxtTree->Branch("Vz",&tree_Vz,"Vz/F");
  FxtTree->Branch("centrality",&tree_centrality,"centrality/s");
  //FxtTree->Branch("refmult",&tree_refmult,"refmult/s");
  FxtTree->Branch("tracknumber",&tree_tracknumber,"tracknumber/s");
  FxtTree->Branch("nSigmaPi",tree_nSigmaPi,"nSigmaPi[tracknumber]/F");
  FxtTree->Branch("nSigmaKa",tree_nSigmaKa,"nSigmaKa[tracknumber]/F");
  FxtTree->Branch("nSigmaPr",tree_nSigmaPr,"nSigmaPr[tracknumber]/F");
  FxtTree->Branch("Charge",tree_Charge,"Charge[tracknumber]/S");
  FxtTree->Branch("Px",tree_Px,"Px[tracknumber]/F");
  FxtTree->Branch("Py",tree_Py,"Py[tracknumber]/F");
  FxtTree->Branch("Pz",tree_Pz,"Pz[tracknumber]/F");
  FxtTree->Branch("DCA",tree_DCA,"DCA[tracknumber]/F");
  FxtTree->Branch("nEPDhits",&tree_nEPDhits,"nEPDhits/s");
  FxtTree->Branch("EPDid",tree_EPDid,"EPDid[nEPDhits]/S");
  FxtTree->Branch("EPDnMip",tree_EPDnMip,"EPDnMip[nEPDhits]/F");
  FxtTree->Branch("tofBeta",tree_Beta,"tofBeta[tracknumber]/F");
  FxtTree->Branch("EtofBeta",tree_EtofBeta,"EtofBeta[tracknumber]/F");
  FxtTree->Branch("dEdx",tree_dEdx,"dEdx[tracknumber]/F");
  FxtTree->Branch("dEdxError",tree_dEdxError,"dEdxError[tracknumber]/F");
  FxtTree->Branch("nHits",tree_nHits,"nHits[tracknumber]/I");
  FxtTree->Branch("nHitsFit",tree_nHitsFit,"nHitsFit[tracknumber]/I");
  FxtTree->Branch("nHitsPoss",tree_nHitsPoss,"nHitsPoss[tracknumber]/I");
  FxtTree->Branch("nHitsDedx",tree_nHitsDedx,"nHitsDedx[tracknumber]/I");


  // QA plots

  // temp variables when histogram bins/bounds depend on the energy
  int tempBins1 = 0;
  double tempLowBound1 = 0;
  double tempHighBound1 = 0;

  h_eventCheck = new TH1D("h_eventCheck","Event number after each cut;;Events", 6, 0, 6);
  h_eventCheck->GetXaxis()->SetBinLabel(1,"No Cuts");
  h_eventCheck->GetXaxis()->SetBinLabel(2,"Bad Run Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(3,"Minbias Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(4,"V_{Z} Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(5,"V_{r} Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(6,"Centrality Cut");

  tempBins1      = (configs.fixed_target) ? 500 : 800;
  tempLowBound1  = (configs.fixed_target) ? 190.0 : -210.0;
  tempHighBound1 = 210.0;
  h_zvtx = new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", tempBins1, tempLowBound1, tempHighBound1);

  h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
    
  h_refmult = new TH1D("h_refmult","Reference multiplicity",1001,-0.5,1000.5);
  h_refmult->GetXaxis()->SetTitle("RefMult");
  h_refmult->GetXaxis()->SetTitle("# of events");
    
  h_tofmult = new TH1D("h_tofmult","TOF multiplicity",1001,-0.5,1000.5);
  h_tofmult->GetXaxis()->SetTitle("TofMult");
  h_tofmult->GetXaxis()->SetTitle("# of events");
    
  h_trackmult = new TH1D("h_trackmult","Actual track multiplicity",1501,-0.5,1500.5);
  h_trackmult->GetXaxis()->SetTitle("TrackMult");
  h_trackmult->GetXaxis()->SetTitle("# of events");
    
  h2_refmult_vs_trackmult = new TH2D("h2_refmult_vs_trackmult","RefMult vs. Actual track multiplicity;TrackMult;RefMult",1501,-0.5,1500.5,1001,-0.5,1000.5);    
  h2_tofmult_vs_trackmult = new TH2D("h2_tofmult_vs_trackmult","TofMult vs. Actual track multiplicity;TrackMult;TofMult",1501,-0.5,1500.5,1001,-0.5,1000.5);
  h2_tofmult_vs_refmult = new TH2D("h2_tofmult_vs_refmult","TofMult vs. RefMult;RefMult;TofMult",1001,-0.5,1000.5,1001,-0.5,1000.5);

  return kStOK;
}

Int_t TreeMaker::Finish()
{
  outputFile->Write();
  return kStOK;
}

Int_t TreeMaker::Make()
{
  h_eventCheck->Fill(0); // Count # of events before any cuts

  StPicoEvent* event = mPicoDstMaker->picoDst()->event(); // Get Event pointer

  Int_t nBTOFMatch = (Int_t)event->nBTOFMatch();

  TVector3 pVtx = event->primaryVertex();
  Double_t d_xvtx = pVtx.x();
  Double_t d_yvtx = pVtx.y();
  Double_t d_zvtx = pVtx.z();
  Double_t d_rvtx = TMath::Sqrt(d_xvtx * d_xvtx + (d_yvtx + 2) * (d_yvtx + 2));

  h_zvtx->Fill(d_zvtx);
  h2_trans_vtx->Fill(d_xvtx,d_yvtx);
    
  if(event) // Ensure event pointer is not NULL
    { 
      if( IsGoodRun(event->runId(), configs.sqrt_s_NN) )  // Select only good runs
	{ 
	  h_eventCheck->Fill(1); // Count # of good run events

	  Bool_t b_good_trig = false;
	  std::vector<UInt_t> triggerIDs = event->triggerIds();
	  
	  for (UInt_t i = 0; i < triggerIDs.size(); i++) 
	    { if ( configs.triggersMatch(triggerIDs[i]) ) {b_good_trig = true;} }
	  
	  if(b_good_trig) // Minbias cut
	    { 
	      h_eventCheck->Fill(2); // Count # of events after minbias trigger id cut
	      	      
	      if(d_zvtx > configs.m_z_vtx_low && d_zvtx < configs.m_z_vtx_high) // Vz cut
		{ 
		  if(d_rvtx < configs.m_r_vtx) // Vr cut
		    { 
		      Int_t nTracks = mPicoDstMaker->picoDst()->numberOfTracks(); // Get the number of Primary Tracks
		      Int_t primTracks = 0;

		      // 1st Primary Tracks loop to get centrality starts
		      for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
			{
			  StPicoTrack* track = mPicoDstMaker->picoDst()->track(iTrk); // Get Track pointer
			  if(!track->isPrimary()) continue; // Only Primary Tracks
			  primTracks++;
			} // 1st Primary tracks loop ends
		      
		      h_refmult->Fill(event->refMultHalfEast());
		      h_trackmult->Fill(primTracks);
		      h_tofmult->Fill(event->nBTOFMatch());
		      h2_refmult_vs_trackmult->Fill(primTracks,event->refMultHalfEast());
		      h2_tofmult_vs_trackmult->Fill(primTracks,event->nBTOFMatch());
		      h2_tofmult_vs_refmult->Fill(event->refMultHalfEast(),event->nBTOFMatch());


		      // GET CENTRALITY
		      Int_t centrality = -99;

		      // 3.0 GeV FXT  --  From Zachary Sweger Nov 11, 2020
		      if (configs.sqrt_s_NN == 3.0) // DON'T FORGET TO UPDATE N_track IN TreeAnalyzer.cxx!!
			{
			  if     ( primTracks >=   5 && primTracks <=   6 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   7 && primTracks <=   8 ) centrality =  1;
			  else if( primTracks >=   9 && primTracks <=  11 ) centrality =  2;
			  else if( primTracks >=  12 && primTracks <=  15 ) centrality =  3;
			  else if( primTracks >=  16 && primTracks <=  20 ) centrality =  4;
			  else if( primTracks >=  21 && primTracks <=  25 ) centrality =  5;
			  else if( primTracks >=  26 && primTracks <=  32 ) centrality =  6;
			  else if( primTracks >=  33 && primTracks <=  40 ) centrality =  7;
			  else if( primTracks >=  41 && primTracks <=  49 ) centrality =  8;
			  else if( primTracks >=  50 && primTracks <=  59 ) centrality =  9;
			  else if( primTracks >=  60 && primTracks <=  71 ) centrality = 10;
			  else if( primTracks >=  72 && primTracks <=  85 ) centrality = 11;
			  else if( primTracks >=  86 && primTracks <= 100 ) centrality = 12;
			  else if( primTracks >= 101 && primTracks <= 118 ) centrality = 13;
			  else if( primTracks >= 119 && primTracks <= 141 ) centrality = 14;
			  else if( primTracks >= 142 && primTracks <= 195 ) centrality = 15;  // 0% - 5% (Central)
			}
		      // 3.2 GeV FXT  --  From Eric Loyd 2022
		      // https://drupal.star.bnl.gov/STAR/pwg/common/BES-II-Centrality-Calibration
		      else if (configs.sqrt_s_NN == 3.2) // DON'T FORGET TO UPDATE N_track IN TreeAnalyzer.cxx!!
			{
			  if     ( primTracks >=   5 && primTracks <=   6 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   7 && primTracks <=  10 ) centrality =  1;
			  else if( primTracks >=  11 && primTracks <=  14 ) centrality =  2;
			  else if( primTracks >=  15 && primTracks <=  19 ) centrality =  3;
			  else if( primTracks >=  20 && primTracks <=  25 ) centrality =  4;
			  else if( primTracks >=  26 && primTracks <=  32 ) centrality =  5;
			  else if( primTracks >=  33 && primTracks <=  42 ) centrality =  6;
			  else if( primTracks >=  43 && primTracks <=  52 ) centrality =  7;
			  else if( primTracks >=  53 && primTracks <=  65 ) centrality =  8;
			  else if( primTracks >=  66 && primTracks <=  80 ) centrality =  9;
			  else if( primTracks >=  81 && primTracks <=  97 ) centrality = 10;
			  else if( primTracks >=  98 && primTracks <= 117 ) centrality = 11;
			  else if( primTracks >= 118 && primTracks <= 140 ) centrality = 12;
			  else if( primTracks >= 141 && primTracks <= 165 ) centrality = 13;
			  else if( primTracks >= 166 && primTracks <= 196 ) centrality = 14;
			  else if( primTracks >= 197 && primTracks <= 287 ) centrality = 15;  // 0% - 5% (Central)
			}
		      // 3.5 GeV FXT  --  From Eric Loyd 2023
		      // https://drupal.star.bnl.gov/STAR/pwg/common/BES-II-Centrality-Calibration
		      else if (configs.sqrt_s_NN == 3.5) // DON'T FORGET TO UPDATE N_track IN TreeAnalyzer.cxx!!
			{
			  if     ( primTracks >=   6 && primTracks <=   7 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   8 && primTracks <=  11 ) centrality =  1;
			  else if( primTracks >=  12 && primTracks <=  15 ) centrality =  2;
			  else if( primTracks >=  16 && primTracks <=  20 ) centrality =  3;
			  else if( primTracks >=  21 && primTracks <=  27 ) centrality =  4;
			  else if( primTracks >=  28 && primTracks <=  36 ) centrality =  5;
			  else if( primTracks >=  37 && primTracks <=  46 ) centrality =  6;
			  else if( primTracks >=  47 && primTracks <=  58 ) centrality =  7;
			  else if( primTracks >=  59 && primTracks <=  72 ) centrality =  8;
			  else if( primTracks >=  73 && primTracks <=  88 ) centrality =  9;
			  else if( primTracks >=  89 && primTracks <= 107 ) centrality = 10;
			  else if( primTracks >= 108 && primTracks <= 127 ) centrality = 11;
			  else if( primTracks >= 128 && primTracks <= 152 ) centrality = 12;
			  else if( primTracks >= 153 && primTracks <= 180 ) centrality = 13;
			  else if( primTracks >= 181 && primTracks <= 215 ) centrality = 14;
			  else if( primTracks >= 216 && primTracks <= 325 ) centrality = 15;  // 0% - 5% (Central)
			}
		      // 3.9 GeV FXT ***2019 RUN ONLY!***  --  From Eric Loyd 2022
		      // https://drupal.star.bnl.gov/STAR/pwg/common/BES-II-Centrality-Calibration
		      else if (configs.sqrt_s_NN == 3.9 && ((int)event->runId() / (int)1000000 == 21)) // DON'T FORGET TO UPDATE N_track IN TreeAnalyzer.cxx!!
			{
			  if     ( primTracks >=   7 && primTracks <=   8 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   9 && primTracks <=  12 ) centrality =  1;
			  else if( primTracks >=  13 && primTracks <=  17 ) centrality =  2;
			  else if( primTracks >=  18 && primTracks <=  23 ) centrality =  3;
			  else if( primTracks >=  24 && primTracks <=  31 ) centrality =  4;
			  else if( primTracks >=  32 && primTracks <=  40 ) centrality =  5;
			  else if( primTracks >=  41 && primTracks <=  51 ) centrality =  6;
			  else if( primTracks >=  52 && primTracks <=  64 ) centrality =  7;
			  else if( primTracks >=  65 && primTracks <=  80 ) centrality =  8;
			  else if( primTracks >=  81 && primTracks <=  97 ) centrality =  9;
			  else if( primTracks >=  98 && primTracks <= 118 ) centrality = 10;
			  else if( primTracks >= 119 && primTracks <= 140 ) centrality = 11;
			  else if( primTracks >= 141 && primTracks <= 167 ) centrality = 12;
			  else if( primTracks >= 168 && primTracks <= 197 ) centrality = 13;
			  else if( primTracks >= 198 && primTracks <= 234 ) centrality = 14;
			  else if( primTracks >= 235 && primTracks <= 344 ) centrality = 15;  // 0% - 5% (Central)
			}
		      // 3.9 GeV FXT ***2020 RUN ONLY!***  --  From Eric Loyd 2023
		      // https://drupal.star.bnl.gov/STAR/pwg/common/BES-II-Centrality-Calibration
		      else if (configs.sqrt_s_NN == 3.9 && ((int)event->runId() / (int)1000000 == 20)) // DON'T FORGET TO UPDATE N_track IN TreeAnalyzer.cxx!!
			{
			  if     ( primTracks >=   6 && primTracks <=   8 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   9 && primTracks <=  12 ) centrality =  1;
			  else if( primTracks >=  13 && primTracks <=  16 ) centrality =  2;
			  else if( primTracks >=  17 && primTracks <=  22 ) centrality =  3;
			  else if( primTracks >=  23 && primTracks <=  30 ) centrality =  4;
			  else if( primTracks >=  31 && primTracks <=  39 ) centrality =  5;
			  else if( primTracks >=  40 && primTracks <=  50 ) centrality =  6;
			  else if( primTracks >=  51 && primTracks <=  63 ) centrality =  7;
			  else if( primTracks >=  64 && primTracks <=  79 ) centrality =  8;
			  else if( primTracks >=  80 && primTracks <=  96 ) centrality =  9;
			  else if( primTracks >=  97 && primTracks <= 117 ) centrality = 10;
			  else if( primTracks >= 118 && primTracks <= 140 ) centrality = 11;
			  else if( primTracks >= 141 && primTracks <= 166 ) centrality = 12;
			  else if( primTracks >= 167 && primTracks <= 197 ) centrality = 13;
			  else if( primTracks >= 198 && primTracks <= 235 ) centrality = 14;
			  else if( primTracks >= 236 && primTracks <= 344 ) centrality = 15;  // 0% - 5% (Central)
			}
		      // 4.5 GeV FXT 
		      // https://drupal.star.bnl.gov/STAR/blog/eloyd/Run-20-45GeVFXT-Centrality-Definition
		      else if (configs.sqrt_s_NN == 4.5) // DON'T FORGET TO UPDATE N_track IN TreeAnalyzer.cxx!!
			{
			  if     ( primTracks >=   7 && primTracks <=   9 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=  10 && primTracks <=  13 ) centrality =  1;
			  else if( primTracks >=  14 && primTracks <=  19 ) centrality =  2;
			  else if( primTracks >=  20 && primTracks <=  26 ) centrality =  3;
			  else if( primTracks >=  27 && primTracks <=  34 ) centrality =  4;
			  else if( primTracks >=  35 && primTracks <=  44 ) centrality =  5;
			  else if( primTracks >=  45 && primTracks <=  56 ) centrality =  6;
			  else if( primTracks >=  57 && primTracks <=  71 ) centrality =  7;
			  else if( primTracks >=  72 && primTracks <=  87 ) centrality =  8;
			  else if( primTracks >=  88 && primTracks <= 107 ) centrality =  9;
			  else if( primTracks >= 108 && primTracks <= 129 ) centrality = 10;
			  else if( primTracks >= 130 && primTracks <= 153 ) centrality = 11;
			  else if( primTracks >= 154 && primTracks <= 181 ) centrality = 12;
			  else if( primTracks >= 182 && primTracks <= 215 ) centrality = 13;
			  else if( primTracks >= 216 && primTracks <= 256 ) centrality = 14;
			  else if( primTracks >= 257 && primTracks <= 367 ) centrality = 15;  // 0% - 5% (Central)
			}

		      // 7.2 GeV FXT
		      else if (configs.sqrt_s_NN == 7.2) // DON'T FORGET TO UPDATE N_track IN TreeAnalyzer.cxx!!
			{
			  if     ( primTracks >=   2 && primTracks <=   3 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   4 && primTracks <=   5 ) centrality =  1;
			  else if( primTracks >=   6 && primTracks <=   8 ) centrality =  2;
			  else if( primTracks >=   9 && primTracks <=  11 ) centrality =  3;
			  else if( primTracks >=  12 && primTracks <=  15 ) centrality =  4;
			  else if( primTracks >=  16 && primTracks <=  21 ) centrality =  5;
			  else if( primTracks >=  22 && primTracks <=  29 ) centrality =  6;
			  else if( primTracks >=  30 && primTracks <=  38 ) centrality =  7;
			  else if( primTracks >=  39 && primTracks <=  49 ) centrality =  8;
			  else if( primTracks >=  50 && primTracks <=  63 ) centrality =  9;
			  else if( primTracks >=  64 && primTracks <=  79 ) centrality = 10;
			  else if( primTracks >=  80 && primTracks <=  99 ) centrality = 11;
			  else if( primTracks >= 100 && primTracks <= 123 ) centrality = 12;
			  else if( primTracks >= 124 && primTracks <= 153 ) centrality = 13;
			  else if( primTracks >= 154 && primTracks <= 190 ) centrality = 14;
			  else if( primTracks >= 191 && primTracks <= 240 ) centrality = 15;  // 0% - 5% (Central)
			}

		      // Select only good centrality valued events
		      if(centrality >= FIRST_CENT)
			{
			  // Prepare New Tree parameters
			  tree_Vx = d_xvtx;
			  tree_Vy = d_yvtx;
			  tree_Vz = d_zvtx;
			  tree_runId   = event->runId();
			  tree_eventId = event->eventId();
			  tree_bField  = event->bField();
			  tree_centrality = centrality;
			  //tree_refmult = event->refMultHalfEast();
			  Int_t realTrackIndex = 0;


			  tree_nEPDhits = mPicoDstMaker->picoDst()->numberOfEpdHits();
			  if(tree_nEPDhits)
			    {
			      for(Int_t hit=0; hit<tree_nEPDhits; hit++)
				{
				  StPicoEpdHit *pEpdHit = mPicoDstMaker->picoDst()->epdHit(hit);
				  if(!pEpdHit) continue;
				  if(!pEpdHit->isGood()) continue;
				  tree_EPDid[hit]   = pEpdHit->id();
				  tree_EPDnMip[hit] = pEpdHit->nMIP();
				}
			    }
			    
			  // 2nd Loop through tracks to fill histograms
			  for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
			    {
			      StPicoTrack* track = mPicoDstMaker->picoDst()->track(iTrk); // Get Track pointer
			      if(!track->isPrimary()) continue;

			      Int_t i_nHits = track->nHits();
			      Int_t i_nHitsDedx = track->nHitsDedx();
			      Int_t i_nHitsFit = track->nHitsFit();
			      Int_t i_nHitsPoss = track->nHitsPoss();
			      Double_t d_DCA = track->gDCA(d_xvtx,d_yvtx,d_zvtx);

			      //=========================================================
			      //          Track QA Cuts
			      //=========================================================
			      bool b_bad_hits  = ( i_nHits < configs.m_nHits );
			      bool b_bad_dEdx  = ( i_nHitsDedx <= configs.m_nHits_dEdx );
			      bool b_bad_ratio = ( ((double)i_nHitsFit / (double)i_nHitsPoss) <= configs.m_nHits_ratio );
			      bool b_bad_DCA   = ( d_DCA >= configs.m_dca );
			      if (b_bad_hits || b_bad_dEdx || b_bad_ratio || b_bad_DCA) continue;
			      //=========================================================
			      //          End Track QA Cuts
			      //=========================================================
			      
			      Float_t d_px = track->pMom().X();
			      Float_t d_py = track->pMom().Y();
			      Float_t d_pz = track->pMom().Z();
			      Float_t d_dEdx = track->dEdx();
			      Float_t d_dEdxError = track->dEdxError();
			      Short_t s_charge = track->charge();
			      Float_t d_nSigmaPi = track->nSigmaPion();
			      Float_t d_nSigmaKa = track->nSigmaKaon();
			      Float_t d_nSigmaPr = track->nSigmaProton();
			      Float_t d_tofBeta = -999.0;
			      Float_t d_EtofBeta = -999.0;

			      if (track->isTofTrack())
				{ d_tofBeta = mPicoDstMaker->picoDst()->btofPidTraits(track->bTofPidTraitsIndex())->btofBeta(); }
			      if (track->isETofTrack())
				{ d_EtofBeta = mPicoDstMaker->picoDst()->etofPidTraits(track->eTofPidTraitsIndex())->beta(); }

			      tree_Charge[realTrackIndex] = s_charge;
			      tree_Px[realTrackIndex] = d_px;
			      tree_Py[realTrackIndex] = d_py;
			      tree_Pz[realTrackIndex] = d_pz;
			      tree_DCA[realTrackIndex] = d_DCA;
			      tree_nSigmaPi[realTrackIndex] = d_nSigmaPi;
			      tree_nSigmaKa[realTrackIndex] = d_nSigmaKa;
			      tree_nSigmaPr[realTrackIndex] = d_nSigmaPr;
                              tree_Beta[realTrackIndex] = d_tofBeta;
			      tree_EtofBeta[realTrackIndex] = d_EtofBeta;
			      tree_dEdx[realTrackIndex] = d_dEdx;
			      tree_dEdxError[realTrackIndex] = d_dEdxError;
			      tree_nHits[realTrackIndex]  = i_nHits;
			      tree_nHitsFit[realTrackIndex]  = i_nHitsFit;
		              tree_nHitsPoss[realTrackIndex] = i_nHitsPoss;
		              tree_nHitsDedx[realTrackIndex] = i_nHitsDedx;
			      realTrackIndex += 1;
			    } // 2nd track loop ends

			  tree_tracknumber = realTrackIndex;
			  FxtTree->Fill();			    
			} // Good centrality event selection ends
                    } // Vr cut ends
                } // Vz cut ends
            } // Minbias cut ends
        } // Good run cut ends
    } // Event pointer not NULL cut ends
  return kStOK;
}

Bool_t TreeMaker::IsGoodRun(Int_t runNumber, Double_t sqrt_s_NN)
{
  // From Ben Kimelman Nov 6, 2020
  const Int_t badRuns_3p0GeV = 24;
  Int_t badRunList_3p0GeV[badRuns_3p0GeV] = {19151029, 19151045, 19152001, 19152078, 19153023, 19153032, 19153065, 19154012, 19154013, 19154014, 19154015, 19154016, 
					     19154017, 19154018, 19154019, 19154020, 19154021, 19154022, 19154023, 19154024, 19154026, 19154046, 19154051, 19154056};

  // From Erik Loyd 2022
  const Int_t badRuns_3p22GeV = 10;
  Int_t badRunList_3p22GeV[badRuns_3p22GeV] = {20180005, 20180006, 20180019, 20180025, 20181016, 20182034, 20183001, 20183013, 20183014, 20183019};

  // FXT above 3.2 GeV:
  // https://drupal.star.bnl.gov/STAR/pwg/common/bes-ii-run-qa/FXT-datasets
  const Int_t badRuns_3p5GeV = 8;
  Int_t badRunList_3p5GeV[badRuns_3p5GeV] = {20355020, 20355021, 21044023, 21045024, 21045025, 21044027, 21044035, 21045004};

  //const Int_t badRuns_3p9GeV = 7;
  //Int_t badRunList_3p9GeV[badRuns_3p9GeV] = {20107029, 20113042, 20113043, 20169033, 20169043, 21035011};
  const Int_t badRuns_3p9GeV = 7;
  Int_t badRunList_3p9GeV[badRuns_3p9GeV] = {20107029, 20113042, 20113043, 20169033, 20169043, 21035011, 21036012};

  const Int_t badRuns_4p5GeV = 1;
  Int_t badRunList_4p5GeV[badRuns_4p5GeV] = {21032001};

  Bool_t b_good_run = true;

  // Run is good if NOT found in the bad run list.
  if (sqrt_s_NN == 3.0)
    b_good_run = !foundInArray(runNumber, badRunList_3p0GeV, badRuns_3p0GeV);
  else if (sqrt_s_NN == 3.2)
    b_good_run = !foundInArray(runNumber, badRunList_3p22GeV, badRuns_3p22GeV);
  else if (sqrt_s_NN == 3.5)
    b_good_run = !foundInArray(runNumber, badRunList_3p5GeV, badRuns_3p5GeV);
  else if (sqrt_s_NN == 3.9)
    b_good_run = !foundInArray(runNumber, badRunList_3p9GeV, badRuns_3p9GeV);
  else if (sqrt_s_NN == 4.5)
    b_good_run = !foundInArray(runNumber, badRunList_4p5GeV, badRuns_4p5GeV);

  return b_good_run;
}

bool TreeMaker::foundInArray(int number, int array[], int size)
{
  bool isFound = false;

  for (int i = 0; i < size; i++)
    {
      if (number == array[i])
	{
	  isFound = true;
	  break;
	}
    }

  return isFound;
}
