// This is an version of Yang Wu's (UCR) tree maker adapted for my own (Cameron Racz) analysis.


// Load header files
#include "TreeMaker.h"
#include "TF1.h"

// Bichsel header
#include "StRoot/StBichsel/Bichsel.h"

// Config header
#include "../ConfigReader/ConfigReader.h"

// Bichsel Function Functions
Double_t bichselZ(Double_t *x,Double_t *par) 
{
  Double_t pove   = x[0];
  Double_t poverm = pove/par[0];
  return TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(poverm),par[1]));
}

Double_t bichsel70(Double_t *x,Double_t *par) 
{
  Double_t pove   = x[0];
  Double_t poverm = pove/par[0];
  return TMath::Exp(Bichsel::Instance()->GetI70(TMath::Log10(poverm),par[1]));
}

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
}

// Destructor
TreeMaker::~TreeMaker() {}

Int_t TreeMaker::Init()
{
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
  FxtTree->Branch("tracknumber",&tree_tracknumber,"tracknumber/s");
  FxtTree->Branch("PID",tree_PID,"PID[tracknumber]/S");
  FxtTree->Branch("Charge",tree_Charge,"Charge[tracknumber]/S");
  FxtTree->Branch("Px",tree_Px,"Px[tracknumber]/F");
  FxtTree->Branch("Py",tree_Py,"Py[tracknumber]/F");
  FxtTree->Branch("Pz",tree_Pz,"Pz[tracknumber]/F");
  FxtTree->Branch("nEPDhits",&tree_nEPDhits,"nEPDhits/s");
  FxtTree->Branch("EPDid",tree_EPDid,"EPDid[nEPDhits]/S");
  FxtTree->Branch("EPDnMip",tree_EPDnMip,"EPDnMip[nEPDhits]/F");
  FxtTree->Branch("tofBeta",tree_Beta,"tofBeta[tracknumber]/D");
  FxtTree->Branch("dEdx",tree_dEdx,"dEdx[tracknumber]/F");
  FxtTree->Branch("nHitsFit",tree_nHitsFit,"nHitsFit[tracknumber]/I");
  FxtTree->Branch("nHitsPoss",tree_nHitsPoss,"nHitsPoss[tracknumber]/I");


  // QA plots
  h_eventCheck = new TH1D("h_eventCheck","Event number after each cut;;Events", 6, 0, 6);
  h_eventCheck->GetXaxis()->SetBinLabel(1,"No Cuts");
  h_eventCheck->GetXaxis()->SetBinLabel(2,"Bad Run Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(3,"Minbias Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(4,"V_{Z} Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(5,"V_{r} Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(6,"Centrality Cut");
    
  h_trackCheck = new TH1D("h_trackCheck","Track number after each cut;;Tracks", 3, 0, 3);
  h_trackCheck->GetXaxis()->SetBinLabel(1,"Event Cuts Only");
  h_trackCheck->GetXaxis()->SetBinLabel(2,"QA Cuts");
  h_trackCheck->GetXaxis()->SetBinLabel(3,"PID Cuts");
    
  h_centralities = new TH1D("h_centralities", "Centralities;Centrality ID;Events", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  h_centralities->GetXaxis()->SetTitle("Centrality bin");
  h_centralities->GetYaxis()->SetTitle("# of events");
    
  //h_zvtx = new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", 100, 190, 210);
  h_zvtx = new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", 840, -210, 210);

  h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  h2_trans_vtx_cut = new TH2D("h2_trans_vtx_cut","Final Primary Vertices;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
    
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

  h_nhits       = new TH1D("h_nhits", "nHits;Number of hits;Tracks", 50, 0, 50);
  h_nhits_dEdx  = new TH1D("h_nhits_dEdx","nHitsdEdx;Number of hits;Tracks", 50, 0, 50);
  h_nhits_ratio = new TH1D("h_nhits_ratio","nhitsFit/nhitsPoss;nhitsFit/nhitsPoss;Tracks",200,0.0,2.0);    
  h_DCA         = new TH1D("h_DCA","DCA (cm);DCA (cm);Tracks",100,0.0,10.0);
    
  h_pT = new TH1D("h_pT","p_{T};p_{T};Tracks",1000,0.0,10.0);
  h_eta = new TH1D("h_eta","#eta;#eta;Tracks",600,-6.0,6.0);
  h_phi = new TH1D("h_phi","#phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h2_dEdx_vs_qp = new TH2D("h2_dEdx_vs_qp", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_half = new TH2D("h2_dEdx_vs_qp_half", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, 0, 4, 500, 0, 12);
  h2_beta_vs_qp = new TH2D("h2_beta_vs_qp","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_m2_vs_qp = new TH2D("h2_m2_vs_qp", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);

 
  h_pT_pp = new TH1D("h_pT_pp","#pi^{+} p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_pm = new TH1D("h_pT_pm","#pi^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_kp = new TH1D("h_pT_kp","K^{+} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_km = new TH1D("h_pT_km","K^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_pr = new TH1D("h_pT_pr","Proton p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_de = new TH1D("h_pT_de","Deuteron p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_tr = new TH1D("h_pT_tr","Triton p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  
  h_eta_pp = new TH1D("h_eta_pp","#pi^{+} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_pm = new TH1D("h_eta_pm","#pi^{-} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_kp = new TH1D("h_eta_kp","K^{+} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_km = new TH1D("h_eta_km","K^{-} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_pr = new TH1D("h_eta_pr","Proton #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_de = new TH1D("h_eta_de","Deuteron #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_tr = new TH1D("h_eta_tr","Triton #eta;#eta;Tracks",500,-5.0,5.0);

  h_dndy_pp = new TH1D("h_dndy_pp", "#pi^{+} Raw Rapidity Spectrum;y;dN/dy", 80, -2, 2);
  h_dndy_pm = new TH1D("h_dndy_pm", "#pi^{-} Raw Rapidity Spectrum;y;dN/dy", 80, -2, 2);
  h_dndy_kp = new TH1D("h_dndy_kp", "K^{+} Raw Rapidity Spectrum;y;dN/dy",   80, -2, 2);
  h_dndy_km = new TH1D("h_dndy_km", "K^{-} Raw Rapidity Spectrum;y;dN/dy",   80, -2, 2);
  h_dndy_pr = new TH1D("h_dndy_pr", "Proton Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);
  h_dndy_de = new TH1D("h_dndy_de", "Deuteron Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);
  h_dndy_tr = new TH1D("h_dndy_tr", "Triton Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);
  
  h_phi_pp = new TH1D("h_phi_pp","#pi^{+} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_pm = new TH1D("h_phi_pm","#pi^{-} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_kp = new TH1D("h_phi_kp","K^{+} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_km = new TH1D("h_phi_km","K^{-} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_pr = new TH1D("h_phi_pr","Proton #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_de = new TH1D("h_phi_de","Deuteron #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_tr = new TH1D("h_phi_tr","Triton #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  
  h_mult_pp = new TH1D("h_mult_pp","#pi^{+} track multiplicity;#pi^{+} Mult;Events",1001,-0.5,1000.5);
  h_mult_pm = new TH1D("h_mult_pm","#pi^{-} track multiplicity;#pi^{-} Mult;Events",1001,-0.5,1000.5);
  h_mult_kp = new TH1D("h_mult_kp","K^{#plus} track multiplicity;K^{+} Mult;Events",1001,-0.5,1000.5);
  h_mult_km = new TH1D("h_mult_km","K^{-} track multiplicity;K^{-} Mult;Events",1001,-0.5,1000.5);
  h_mult_pr = new TH1D("h_mult_pr","Proton track multiplicity;Proton Mult;Events",1001,-0.5,1000.5);
  h_mult_de = new TH1D("h_mult_de","Deuteron track multiplicity;Deuteron Mult;Events",1001,-0.5,1000.5);
  h_mult_tr = new TH1D("h_mult_tr","Triton track multiplicity;Triton Mult;Events",1001,-0.5,1000.5);
  
  h2_dEdx_vs_qp_pp = new TH2D("h2_dEdx_vs_qp_pp", "#pi^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_pm = new TH2D("h2_dEdx_vs_qp_pm", "#pi^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_kp = new TH2D("h2_dEdx_vs_qp_kp", "K^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_km = new TH2D("h2_dEdx_vs_qp_km", "K^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_pr = new TH2D("h2_dEdx_vs_qp_pr", "Proton dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_de = new TH2D("h2_dEdx_vs_qp_de", "Deuteron dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_tr = new TH2D("h2_dEdx_vs_qp_tr", "Triton dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  
  h2_beta_vs_qp_pp = new TH2D("h2_beta_vs_qp_pp","#pi^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_pm = new TH2D("h2_beta_vs_qp_pm","#pi^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_kp = new TH2D("h2_beta_vs_qp_kp","K^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_km = new TH2D("h2_beta_vs_qp_km","K^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_pr = new TH2D("h2_beta_vs_qp_pr","Proton 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_de = new TH2D("h2_beta_vs_qp_de","Deuteron 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_tr = new TH2D("h2_beta_vs_qp_tr","Triton 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
 
  h2_m2_vs_qp_pp = new TH2D("h2_m2_vs_qp_pp", "#pi^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_pm = new TH2D("h2_m2_vs_qp_pm", "#pi^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_kp = new TH2D("h2_m2_vs_qp_kp", "K^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_km = new TH2D("h2_m2_vs_qp_km", "K^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_pr = new TH2D("h2_m2_vs_qp_pr", "Proton m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_de = new TH2D("h2_m2_vs_qp_de", "Deuteron m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_tr = new TH2D("h2_m2_vs_qp_tr", "Triton m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);

  h2_pT_vs_yCM_pp = new TH2D("h2_pT_vs_yCM_pp", "#pi^{+};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_pm = new TH2D("h2_pT_vs_yCM_pm", "#pi^{-};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_kp = new TH2D("h2_pT_vs_yCM_kp", "K^{+};y-y_{mid};p_{T} (GeV/c)",   300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_km = new TH2D("h2_pT_vs_yCM_km", "K^{-};y-y_{mid};p_{T} (GeV/c)",   300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_pr = new TH2D("h2_pT_vs_yCM_pr", "Proton;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_de = new TH2D("h2_pT_vs_yCM_de", "Deuteron;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_tr = new TH2D("h2_pT_vs_yCM_tr", "Triton;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);

  return kStOK;
}

Int_t TreeMaker::Finish()
{
  outputFile->Write();
  return kStOK;
}

Int_t TreeMaker::Make()
{
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
  //=========================================================
  //          END Bichsel Function Setup
  //=========================================================


  h_eventCheck->Fill(0); // Count # of events before any cuts

  StPicoEvent* event = mPicoDstMaker->picoDst()->event(); // Get Event pointer
    
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
	      
	      TVector3 pVtx = event->primaryVertex();
	      Double_t d_xvtx = pVtx.x();
	      Double_t d_yvtx = pVtx.y();
	      Double_t d_zvtx = pVtx.z();
	      Double_t d_rvtx = TMath::Sqrt(d_xvtx * d_xvtx + (d_yvtx + 2) * (d_yvtx + 2));
	      
	      h_zvtx->Fill(pVtx.Z());
	      h2_trans_vtx->Fill(pVtx.X(),pVtx.Y());
	      
	      if(d_zvtx > configs.z_vtx_low && d_zvtx < configs.z_vtx_high) // Vz cut
		{ 
		  h_eventCheck->Fill(3);    // Count # of events after Vz cut
		  
		  if(d_rvtx < configs.r_vtx) // Vr cut
		    { 
		      h_eventCheck->Fill(4); // Count # of events after Vr cut
		      
		      Int_t nTracks = mPicoDstMaker->picoDst()->numberOfTracks(); // Get the number of Primary Tracks
		      Int_t primTracks = 0;
		      Int_t N_pp = 0;
		      Int_t N_pm = 0;
		      Int_t N_kp = 0;
		      Int_t N_km = 0;
		      Int_t N_pr = 0;
		      Int_t N_de = 0;
		      Int_t N_tr = 0;

		      // 1st Primary Tracks loop to get centrality starts
		      for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
			{
			  StPicoTrack* track = mPicoDstMaker->picoDst()->track(iTrk); // Get Track pointer
			  if(!track->isPrimary()) continue; // Only Primary Tracks
			  primTracks++;
			  h_nhits->Fill(track->nHits());
			  h_nhits_dEdx->Fill(track->nHitsDedx());
			  h_nhits_ratio->Fill((Double_t)track->nHitsFit()/track->nHitsPoss());
			  h_DCA->Fill(track->gDCA(pVtx).Mag());
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
		      if (configs.sqrt_s_NN == 3.0)
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

		      // 7.2 GeV FXT
		      else if (configs.sqrt_s_NN == 7.2)
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

		      // 19.6 GeV COL
		      else if (configs.sqrt_s_NN == 19.6)
			{
			  if     ( primTracks >=   4 && primTracks <=   5 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   6 && primTracks <=   8 ) centrality =  1;
			  else if( primTracks >=   9 && primTracks <=  11 ) centrality =  2;
			  else if( primTracks >=  12 && primTracks <=  16 ) centrality =  3;
			  else if( primTracks >=  17 && primTracks <=  23 ) centrality =  4;
			  else if( primTracks >=  24 && primTracks <=  31 ) centrality =  5;
			  else if( primTracks >=  32 && primTracks <=  42 ) centrality =  6;
			  else if( primTracks >=  43 && primTracks <=  56 ) centrality =  7;
			  else if( primTracks >=  57 && primTracks <=  72 ) centrality =  8;
			  else if( primTracks >=  73 && primTracks <=  93 ) centrality =  9;
			  else if( primTracks >=  94 && primTracks <= 118 ) centrality = 10;
			  else if( primTracks >= 119 && primTracks <= 149 ) centrality = 11;
			  else if( primTracks >= 150 && primTracks <= 186 ) centrality = 12;
			  else if( primTracks >= 187 && primTracks <= 232 ) centrality = 13;
			  else if( primTracks >= 233 && primTracks <= 290 ) centrality = 14;
			  else if( primTracks >= 291 && primTracks <= 2048) centrality = 15;  // 0% - 5% (Central)
			}


		      // Select only good centrality valued events
		      if(centrality >= FIRST_CENT)
			{
			  h_eventCheck->Fill(5); // Count # of min-bias events after all cuts
			  h_centralities->Fill(centrality);

			  // Prepare New Tree parameters
			  tree_Vx = d_xvtx;
			  tree_Vy = d_yvtx;
			  tree_Vz = d_zvtx;
			  tree_runId   = event->runId();
			  tree_eventId = event->eventId();
			  tree_bField  = event->bField();
			  tree_centrality = centrality;
			  /*tree_tracknumber = primTracks;*/
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

			      h_trackCheck->Fill(0);  // Only event cuts
			      
			      //=========================================================
			      //          Track QA Cuts
			      //=========================================================
			      bool b_bad_hits  = ( track->nHits() < configs.nHits );
			      bool b_bad_dEdx  = ( track->nHitsDedx() <= configs.nHits_dEdx );
			      bool b_bad_ratio = ( ((double)track->nHitsFit() / (double)track->nHitsPoss()) <= configs.nHits_ratio );
			      bool b_bad_DCA   = ( track->gDCA(d_xvtx,d_yvtx,d_zvtx) >= configs.dca );

			      if (b_bad_hits || b_bad_dEdx || b_bad_ratio || b_bad_DCA) continue;
			      //=========================================================
			      //          End Track QA Cuts
			      //=========================================================
			      
			      h_trackCheck->Fill(1); // After QA Cuts
			      
			      Double_t d_pT = track->pPt();
			      Double_t d_px = track->pMom().X();
			      Double_t d_py = track->pMom().Y();
			      Double_t d_pz = track->pMom().Z();
			      Double_t d_eta = track->pMom().Eta();
			      Double_t d_phi = track->pMom().Phi();
			      Double_t d_mom = track->pMom().Mag();
			      Double_t d_dEdx = track->dEdx();
			      Short_t  s_charge = track->charge();
			      Double_t d_TPCnSigmaPion   = track->nSigmaPion();
			      Double_t d_TPCnSigmaKaon   = track->nSigmaKaon();
			      Double_t d_TPCnSigmaProton = track->nSigmaProton();
			      Double_t d_zDeuteron = (s_charge > 0) ? TMath::Log(d_dEdx / bichselZ_de->Eval(d_mom)) : -999.0;
			      Double_t d_zTriton   = (s_charge > 0) ? TMath::Log(d_dEdx / bichselZ_tr->Eval(d_mom)) : -999.0;

			      h_phi->Fill(d_phi);
			      h_eta->Fill(d_eta);
			      h_pT->Fill(d_pT);
			      h2_dEdx_vs_qp->Fill(s_charge * d_mom, track->dEdx());
			      
			      Bool_t tofTrack = false;
			      Double_t d_tofBeta = -999.0;
			      Double_t d_m2 = 0.0;
			      Int_t trackTofIndex = track->bTofPidTraitsIndex();
			      if(trackTofIndex >= 0)
				d_tofBeta = mPicoDstMaker->picoDst()->btofPidTraits(trackTofIndex)->btofBeta();
			      
			      if(d_tofBeta != -999.0)
				{
				  tofTrack = true;
				  d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
				  h2_beta_vs_qp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				  h2_m2_vs_qp->Fill(s_charge*d_mom, d_m2);
				}

			      //=========================================================
			      //          PID Cuts
			      //=========================================================
			      Bool_t pion   = false;
			      Bool_t kaon   = false;
			      Bool_t proton = (d_TPCnSigmaProton > configs.nSig_pr_low) && (d_TPCnSigmaProton < configs.nSig_pr_high) && (s_charge > 0);
			      Bool_t deuteron = false;
			      Bool_t triton   = false;
			      //Bool_t deuteron = (d_zDeuteron > configs.z_de_low) && (d_zDeuteron < configs.z_de_high);
			      //Bool_t triton   = (d_zTriton > configs.z_tr_low) && (d_zTriton < configs.z_tr_high);
			      // d_zDeuteron and d_zDeuteron already ensure that s_charge > 0

			      if (tofTrack)
				{
				  pion = (d_TPCnSigmaPion > configs.nSig_pi_low) &&
				         (d_TPCnSigmaPion < configs.nSig_pi_high) &&
				         (d_m2 > configs.m2_pi_low) &&
				         (d_m2 < configs.m2_pi_high);

				  kaon = (d_TPCnSigmaKaon > configs.nSig_ka_low) &&
				         (d_TPCnSigmaKaon < configs.nSig_ka_high) &&
				         (d_m2 > configs.m2_ka_low) &&
				         (d_m2 < configs.m2_ka_high);

				  deuteron = (d_zDeuteron > configs.z_de_low) &&
				             (d_zDeuteron < configs.z_de_high) &&
				             (d_m2 > configs.m2_de_low) &&
				             (d_m2 < configs.m2_de_high);

				  triton   = (d_zTriton > configs.z_tr_low) &&
				             (d_zTriton < configs.z_tr_high) &&
				             (d_m2 > configs.m2_tr_low) &&
				             (d_m2 < configs.m2_tr_high);
				}

			      //if (deuteron) h2_dEdx_vs_qp_half_postZdCut->Fill(s_charge * d_mom, d_dEdx);
			      //if (triton)   h2_dEdx_vs_qp_half_postZtCut->Fill(s_charge * d_mom, d_dEdx);
	    
			      //if (!pion && !kaon && !proton) continue;

			      if (pion     && proton) { proton = false; }
			      if (kaon     && proton) { proton = false; }
			      if (deuteron && proton) { proton = false; }
			      if (triton   && proton) { proton = false; }

			      //if (pion && deuteron) { deuteron = false; }
			      //if (pion && triton)   { triton   = false; }
			      //if (kaon && deuteron) { deuteron = false; }
			      //if (kaon && triton)   { triton   = false; }

			      /*
				if (deuteron && proton) 
				{ 
				if (TMath::Abs(d_zDeuteron) < TMath::Abs(d_TPCnSigmaProton)) { proton = false; }
				else if (TMath::Abs(d_zDeuteron) == TMath::Abs(d_TPCnSigmaProton)) { proton = false; deuteron = false; }
				else { deuteron = false; }
				}
				if (triton && proton)
				{
				if (TMath::Abs(d_zTriton) < TMath::Abs(d_TPCnSigmaProton)) { proton = false; }
				else if (TMath::Abs(d_zTriton) == TMath::Abs(d_TPCnSigmaProton)) { proton = false; triton = false; }
				else { triton = false; }
				}
				if (deuteron && triton)
				{
				if (TMath::Abs(d_zDeuteron) < TMath::Abs(d_zTriton)) { triton = false; }
				else if (TMath::Abs(d_zDeuteron) == TMath::Abs(d_zTriton)) { triton = false; deuteron = false; }
				else { deuteron = false; }
				}
			      */
			      //=========================================================
			      //          END PID Cuts
			      //=========================================================

			      if (pion || kaon || proton || deuteron || triton) h_trackCheck->Fill(2);

			      Double_t mRapidity = 999.0;

			      if(pion) // PID Pions
				{ 
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_PI);
				  
				  if(s_charge > 0)
				    {
				      N_pp++;
				      h_eta_pp->Fill(d_eta);
				      h_phi_pp->Fill(d_phi);
				      h_pT_pp->Fill(d_pT);
				      h_dndy_pp->Fill(mRapidity);
				      h2_pT_vs_yCM_pp->Fill(mRapidity - Y_MID, d_pT);
				      h2_dEdx_vs_qp_pp->Fill(s_charge*d_mom, track->dEdx());
				      h2_beta_vs_qp_pp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				      h2_m2_vs_qp_pp->Fill(s_charge*d_mom, d_m2);
				    }
				  else if(s_charge < 0)
				    {
				      N_pm++;
				      h_eta_pm->Fill(d_eta);
				      h_phi_pm->Fill(d_phi);
				      h_pT_pm->Fill(d_pT);
				      h_dndy_pm->Fill(mRapidity);
				      h2_pT_vs_yCM_pm->Fill(mRapidity - Y_MID, d_pT);
				      h2_dEdx_vs_qp_pm->Fill(s_charge*d_mom, track->dEdx());
				      h2_beta_vs_qp_pm->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				      h2_m2_vs_qp_pm->Fill(s_charge*d_mom, d_m2);
				    }
				}
			      else if(kaon) // PID Kaons
				{
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_KA);
				  
				  if(s_charge > 0)
				    {
				      N_kp++;
				      h_eta_kp->Fill(d_eta);
				      h_phi_kp->Fill(d_phi);
				      h_pT_kp->Fill(d_pT);
				      h_dndy_kp->Fill(mRapidity);
				      h2_pT_vs_yCM_kp->Fill(mRapidity - Y_MID, d_pT);
				      h2_dEdx_vs_qp_kp->Fill(s_charge*d_mom, track->dEdx());
				      h2_beta_vs_qp_kp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				      h2_m2_vs_qp_kp->Fill(s_charge*d_mom, d_m2);
				    }
				  if(s_charge < 0)
				    {
				      N_km++;
				      h_eta_km->Fill(d_eta);
				      h_phi_km->Fill(d_phi);
				      h_pT_km->Fill(d_pT);
				      h_dndy_km->Fill(mRapidity);
				      h2_pT_vs_yCM_km->Fill(mRapidity - Y_MID, d_pT);
				      h2_dEdx_vs_qp_km->Fill(s_charge*d_mom, track->dEdx());
				      h2_beta_vs_qp_km->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				      h2_m2_vs_qp_km->Fill(s_charge*d_mom, d_m2);
				    }
				}
			      else if(proton) // PID Proton
				{
				  N_pr++;
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_PR);
				  
				  h_eta_pr->Fill(d_eta);
				  h_phi_pr->Fill(d_phi);
				  h_pT_pr->Fill(d_pT);
				  h_dndy_pr->Fill(mRapidity);
				  h2_pT_vs_yCM_pr->Fill(mRapidity - Y_MID, d_pT);
				  h2_dEdx_vs_qp_pr->Fill(s_charge*d_mom, track->dEdx());
				  h2_beta_vs_qp_pr->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				  h2_m2_vs_qp_pr->Fill(s_charge*d_mom, d_m2);
				}		    
			      else if(deuteron) // PID Deuteron
				{
				  N_de++;
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_DE);
				  
				  h_eta_de->Fill(d_eta);
				  h_phi_de->Fill(d_phi);
				  h_pT_de->Fill(d_pT);
				  h_dndy_de->Fill(mRapidity);
				  h2_pT_vs_yCM_de->Fill(mRapidity - Y_MID, d_pT);
				  h2_dEdx_vs_qp_de->Fill(s_charge*d_mom, track->dEdx());
				  h2_beta_vs_qp_de->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				  h2_m2_vs_qp_de->Fill(s_charge*d_mom, d_m2);
				}		    
			      else if(triton) // PID Triton
				{
				  N_tr++;
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_TR);
				  
				  h_eta_tr->Fill(d_eta);
				  h_phi_tr->Fill(d_phi);
				  h_pT_tr->Fill(d_pT);
				  h_dndy_tr->Fill(mRapidity);
				  h2_pT_vs_yCM_tr->Fill(mRapidity - Y_MID, d_pT);
				  h2_dEdx_vs_qp_tr->Fill(s_charge*d_mom, track->dEdx());
				  h2_beta_vs_qp_tr->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				  h2_m2_vs_qp_tr->Fill(s_charge*d_mom, d_m2);
				}		    

			      tree_Charge[realTrackIndex] = s_charge;
			      if(pion)          tree_PID[realTrackIndex] = 0;
			      else if(kaon)     tree_PID[realTrackIndex] = 1;
			      else if(proton)   tree_PID[realTrackIndex] = 2;
			      else if(deuteron) tree_PID[realTrackIndex] = 3;
			      else if(triton)   tree_PID[realTrackIndex] = 4;
			      else              tree_PID[realTrackIndex] = -1;

			      tree_Px[realTrackIndex] = d_px;
			      tree_Py[realTrackIndex] = d_py;
			      tree_Pz[realTrackIndex] = d_pz;
                              tree_Beta[realTrackIndex] = d_tofBeta;
			      tree_dEdx[realTrackIndex] = track->dEdx();
		              tree_nHitsFit[realTrackIndex] = track->nHitsFit();
		              tree_nHitsPoss[realTrackIndex] = track->nHitsPoss();
			      realTrackIndex += 1;
			    } // 2nd track loop ends
			    
			  h_mult_pp->Fill(N_pp);
			  h_mult_pm->Fill(N_pm);
			  h_mult_kp->Fill(N_kp);
			  h_mult_km->Fill(N_km);
			  h_mult_pr->Fill(N_pr);
			  h_mult_de->Fill(N_de);
			  h_mult_tr->Fill(N_tr);
			  //tree_tracknumber = N_pr + N_pp + N_pm + N_kp + N_km + N_de + N_tr;
			  tree_tracknumber = realTrackIndex + 1;
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
  Int_t badRunList_3p0GeV[24] = {19151029, 19151045, 19152001, 19152078, 19153023, 19153032, 19153065, 19154012, 19154013, 19154014, 19154015, 19154016, 
    19154017, 19154018, 19154019, 19154020, 19154021, 19154022, 19154023, 19154024, 19154026, 19154046, 19154051, 19154056};

  Int_t badRunList_19p6GeV[374] = {20057007, 20057025, 20057026, 20057050, 20058001, 20058002, 20058003, 20058004, 20058005, 20060012, 
				   20060022, 20060025, 20060060, 20060061, 20060062, 20062010, 20062011, 20062012, 20062036, 20063011, 
				   20063034, 20063035, 20063036, 20063039, 20064008, 20064009, 20064011, 20064012, 20064040, 20065018, 
				   20067014, 20067023, 20067024, 20067029, 20067030, 20067045, 20067046, 20069030, 20069032, 20069054, 
				   20070042, 20070043, 20070044, 20070047, 20071001, 20071004, 20071005, 20071006, 20071027, 20071037, 
				   20072034, 20072035, 20072036, 20072039, 20072041, 20072045, 20072047, 20073071, 20073072, 20073076, 
				   20074001, 20074003, 20074004, 20074005, 20074007, 20074008, 20074009, 20074012, 20074014, 20074017, 
				   20074018, 20074020, 20074021, 20074026, 20074027, 20074029, 20074032, 20074033, 20074034, 20074044, 
				   20074045, 20075001, 20075002, 20075006, 20075007, 20075009, 20075011, 20075013, 20081002, 20081014, 
				   20082060, 20082065, 20083024, 20086012, 20087007, 20089008, 20090024, 20091011, 20092054,
				   20062007, 20062009, 20065017, 20065056, 20065060, 20066001, 20066008, 20066015, 20066019, 20066023, 
				   20066026, 20066066, 20066067, 20066068, 20066073, 20066078, 20067001, 20067004, 20067009, 20067012, 
				   20067015, 20067016, 20067019, 20067028, 20067038, 20067041, 20067047, 20068001, 20068004, 20068008,
				   20068012, 20068019, 20068026, 20068034, 20068051, 20068055, 20068058, 20068060, 20068064, 20069001, 
				   20069004, 20069007, 20069010, 20069020, 20069023, 20069026, 20069031, 20069033, 20069042, 20069050, 
				   20069053, 20069057, 20069060, 20070002, 20070005, 20070010, 20070013, 20070016, 20070019, 20070041,
				   20070045, 20071003, 20071007, 20071010, 20071013, 20071016, 20071019, 20071024, 20071029, 20071036,
				   20071041, 20071044, 20071047, 20071050, 20071053, 20071056, 20071059, 20071063, 20072002, 20072005,
				   20072009, 20072012, 20072016, 20072037, 20072038, 20072046, 20072050, 20072055, 20073002, 20073006,
				   20073013, 20073017, 20073022, 20073025, 20073074, 20074002, 20074006, 20074010, 20074011, 20074016, 
				   20074019, 20074023, 20074030, 20074043, 20074046, 20075004, 20075008, 20075014, 20075010, 20075015,
				   20075020, 20075025, 20075031, 20075035, 20075039, 20075043, 20075048, 20075054, 20075057, 20075060, 
				   20075066, 20076001, 20076004, 20076007, 20076010, 20076013, 20076017, 20076021, 20076025, 20076028, 
				   20076031, 20076034, 20076037, 20076040, 20076045, 20076048, 20076051, 20076054, 20076059, 20077002, 
				   20077005, 20077008, 20077011, 20077014, 20077017, 20077018, 20078001, 20078007, 20078013, 20078016, 
				   20078019, 20078022, 20078028, 20078032, 20078035, 20078040, 20078043, 20078046, 20078051, 20078054, 
				   20078057, 20078060, 20078063, 20078067, 20079006, 20079009, 20079013, 20079017, 20079020, 20079023, 
				   20079044, 20080006, 20080009, 20080012, 20080016, 20080020, 20081001, 20081004, 20081007, 20081012, 
				   20081015, 20081018, 20081025, 20082002, 20082005, 20082010, 20082013, 20082016, 20082019, 20082024, 
				   20082029, 20082034, 20082038, 20082047, 20082050, 20082053, 20082056, 20082059, 20082063, 20082066, 
				   20083001, 20083004, 20083019, 20083022, 20083025, 20083029, 20083032, 20083074, 20083077, 20083079,
				   20084001, 20084002, 20084005, 20084009, 20084013, 20084016, 20084022, 20085006, 20085009, 20085017, 
				   20086002, 20086005, 20086056, 20086011, 20086015, 20087008, 20087012, 20087021, 20088005, 20088009, 
				   20088012, 20088030, 20088033, 20088037, 20089003, 20089006, 20089009, 20089012, 20089015, 20089018, 
				   20089028, 20090002, 20090005, 20090008, 20090011, 20090014, 20090017, 20090021, 20090031, 20090048, 
				   20091003, 20091006, 20091009, 20091012, 20091016, 20091019, 20091020, 20092005, 20092012, 20092015, 
				   20092018, 20092021, 20092024, 20092027, 20092030, 20092033, 20092038, 20092053, 20092057, 20093001, 
				   20093005, 20093010, 20093016, 20093025, 20093035};

  Bool_t b_good_run = true;
  
  if (sqrt_s_NN == 3.0)
    { 
      for (Int_t i = 0; i < 24; i++) 
	{ if (runNumber == badRunList_3p0GeV[i]) {b_good_run = false; break;} } 
    }
  if (sqrt_s_NN == 19.6)
    { 
      for (Int_t i = 0; i < 374; i++) 
	{ if (runNumber == badRunList_19p6GeV[i]) {b_good_run = false; break;} } 
    }

  //for (Int_t i = 0; i < 24; i++) { if (runNumber == badRunList_3p0GeV[i]) {b_good_run = false; break;} } 
  
  return b_good_run;
}


Double_t TreeMaker::getRapidity(Double_t px, Double_t py, Double_t pz, Double_t mass)
{
  Double_t rapidity, energy, momentum = 0;
  momentum = TMath::Sqrt(px*px + py*py + pz*pz);
  energy   = TMath::Sqrt(momentum*momentum + mass*mass);
  rapidity = TMath::ATanH(pz/energy);
  return rapidity;
}
