#include "PlotUtils.h"

void prelimRapidityPlots(TString jobID, Bool_t useSystematics = false, Double_t sqrt_s_NN = 3.0)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1000, 1000);
  //canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetLeftMargin(0.15);
  canvas->SetTopMargin(0.04);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.1);
  canvas->cd();

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetLineWidth(3);
  //gStyle->SetOptDate();

  TProfile2D *p2_vn_yCM_cent_pr = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr");
  TProfile2D *p2_vn_yCM_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr_symmetry");

  TProfile *p_vn_yCM_00to10_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_00to10_pr", 15, 16);
  TProfile *p_vn_yCM_10to40_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_10to40_pr", 9, 14);
  TProfile *p_vn_yCM_40to60_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_40to60_pr", 5, 8);
  TProfile *p_vn_yCM_00to10_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_00to10_pr_symm", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_10to40_pr_symm", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_40to60_pr_symm", 5, 8);

  TH1D *h_vn_yCM_00to10_pr = new TH1D("h_vn_yCM_00to10_pr", ";y-y_{mid};v_{3}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr = new TH1D("h_vn_yCM_10to40_pr", ";y-y_{mid};v_{3}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr = new TH1D("h_vn_yCM_40to60_pr", ";y-y_{mid};v_{3}", 20, -1, 1);
  TH1D *h_vn_yCM_00to10_pr_symm = new TH1D("h_vn_yCM_00to10_pr_symm", ";y-y_{mid};v_{3}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr_symm = new TH1D("h_vn_yCM_10to40_pr_symm", ";y-y_{mid};v_{3}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr_symm = new TH1D("h_vn_yCM_40to60_pr_symm", ";y-y_{mid};v_{3}", 20, -1, 1);

  // Mirrored Plots
  TH1D *h_vn_yCM_00to10_pr_mirror = new TH1D("h_vn_yCM_00to10_pr_mirror", ";y-y_{mid};v_{3}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr_mirror = new TH1D("h_vn_yCM_10to40_pr_mirror", ";y-y_{mid};v_{3}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr_mirror = new TH1D("h_vn_yCM_40to60_pr_mirror", ";y-y_{mid};v_{3}", 20, -1, 1);

  // Convert profiles to histograms
  h_vn_yCM_00to10_pr = p_vn_yCM_00to10_pr->ProjectionX();
  h_vn_yCM_10to40_pr = p_vn_yCM_10to40_pr->ProjectionX();
  h_vn_yCM_40to60_pr = p_vn_yCM_40to60_pr->ProjectionX();
  h_vn_yCM_00to10_pr_symm = p_vn_yCM_00to10_pr_symm->ProjectionX();
  h_vn_yCM_10to40_pr_symm = p_vn_yCM_10to40_pr_symm->ProjectionX();
  h_vn_yCM_40to60_pr_symm = p_vn_yCM_40to60_pr_symm->ProjectionX();

  // Retrieve systematic uncertainties
  TFile* systematicFile;
  TGraphErrors* sys_yCM_00to10_pr;
  TGraphErrors* sys_yCM_10to40_pr;
  TGraphErrors* sys_yCM_40to60_pr;
  TGraphErrors* sys_yCM_00to10_pr_symm;
  TGraphErrors* sys_yCM_10to40_pr_symm;
  TGraphErrors* sys_yCM_40to60_pr_symm;

  if (useSystematics)
    {
      systematicFile = TFile::Open("systematicErrors_thirdDraft.root", "READ");
      sys_yCM_00to10_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_00to10_pr_px");
      sys_yCM_10to40_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_10to40_pr_px");
      sys_yCM_40to60_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_40to60_pr_px");
      sys_yCM_00to10_pr_symm = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_00to10_pr_symm_px");
      sys_yCM_10to40_pr_symm = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_10to40_pr_symm_px");
      sys_yCM_40to60_pr_symm = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_40to60_pr_symm_px");
    }
  ////


  std::cout << std::endl << "Fitting 10-40% non-symmetric: " << std::endl;

  TF1* function = new TF1("function", "[0]*x + [1]*x*x*x", 0.0, 1.0);
  TH1D* h_10to40 = (TH1D*)h_vn_yCM_10to40_pr->Clone();
  h_10to40->Fit(function, "", "", 0.0, 1.0);
  Double_t slopeMeasured = function->GetParameter(0);

  std::cout << "Measured slope = " << slopeMeasured << std::endl
	    << "Chi^2/NDF = " << function->GetChisquare() / function->GetNDF() << std::endl
	    << std::endl
	    << "done"
	    << std::endl;

  std::cout << std::endl << std::endl << "Now estimating uncertainty from symmetric plot..."
	    << std::endl << "Fitting 10-40% symmetric: " << std::endl;
  
  TF1* function1 = new TF1("function1", "[0]*x + [1]*x*x*x", 0.0, 0.5);
  function1->SetLineColor(1);

  TF1* function2 = new TF1("function2", "[0]*x + [1]*x*x*x", -0.5, 0.0);
  function2->SetParLimits(1, 0.0, 2.0);
  function2->SetLineColor(1);

  TF1* function3 = new TF1("function3", "[0]*x + [1]*x*x*x", -0.5, 0.5);
  function3->SetLineColor(1);

  Double_t slopeLeft = 0.0;
  Double_t slopeRight = 0.0;
  Double_t slopeDiff = 0.0;
  Double_t slopeAvg = 0.0;
  Double_t slopeUncertPercent = 0.0;


  TH1D* h_10to40_symm = (TH1D*)h_vn_yCM_10to40_pr_symm->Clone();

  // Add systematics with statistics
  if (useSystematics)
    {
      for (int i = 6; i <= 15; i++)
	{
	  Double_t currentError1 = h_10to40_symm->GetBinError(i);
	  h_10to40_symm->SetBinError(i, std::sqrt(std::pow(currentError1, 2) +
						  std::pow(sys_yCM_10to40_pr_symm->GetErrorY(i-1), 2)));
	}
    }

  h_10to40_symm->Fit(function1, "", "", 0.0, 0.5);
  slopeRight = function1->GetParameter(0);
  std::cout << "Fitting y in [0.0, 0.5]:" << std::endl
	    << "   Slope = " << function1->GetParameter(0) << std::endl
	    << "   Error = " << function1->GetParError(0) << std::endl
	    << "   Chi^2/NDF = " << function1->GetChisquare()/function1->GetNDF() << std::endl
	    << std::endl;

  h_10to40_symm->Fit(function2, "", "", -0.5, 0.0);
  slopeLeft = function2->GetParameter(0);
  std::cout << "Fitting y in [-0.5, 0.0]:" << std::endl
	    << "   Slope = " << function2->GetParameter(0) << std::endl
	    << "   Error = " << function2->GetParError(0) << std::endl
	    << "   Chi^2/NDF = " << function2->GetChisquare()/function2->GetNDF() << std::endl
	    << std::endl;

  slopeDiff = TMath::Abs(slopeRight - slopeLeft);
  slopeAvg  = TMath::Abs((slopeRight + slopeLeft)/2.0);
  slopeUncertPercent = slopeDiff * 100.0 / (slopeAvg * TMath::Sqrt(12));

  std::cout << "Estimated Percent Uncertainty = " << slopeUncertPercent << "%"
	    << std::endl
	    << "Estimated Uncertainty Value = " << TMath::Abs(slopeMeasured) * (slopeUncertPercent/100.0)
	    << std::endl
	    << std::endl;


  
  std::cout << std::endl << std::endl << "Now fitting all symmetric data with an offset linear fit..." << std::endl;

  TF1* function4 = new TF1("function4", "[0]*x + [1]", -0.5, 0.5);

  std::cout << std::endl << "Fitting 0-10% symmetric: " << std::endl;
  TH1D* h_00to10_symm = (TH1D*)h_vn_yCM_00to10_pr_symm->Clone();
  h_00to10_symm->Fit(function4, "", "", -0.5, 0.5);
  std::cout << "Chi^2/NDF = " << function4->GetChisquare() << "/" << function4->GetNDF() << " = " << function4->GetChisquare() / function4->GetNDF() << std::endl;

  std::cout << std::endl << "Fitting 10-40% symmetric: " << std::endl;
  h_10to40_symm = (TH1D*)h_vn_yCM_10to40_pr_symm->Clone();
  h_10to40_symm->Fit(function4, "", "", -0.5, 0.5);
  std::cout << "Chi^2/NDF = " << function4->GetChisquare() << "/" << function4->GetNDF() << " = " << function4->GetChisquare() / function4->GetNDF() << std::endl;

  std::cout << std::endl << "Fitting 40-60% symmetric: " << std::endl;  
  TH1D* h_40to60_symm = (TH1D*)h_vn_yCM_40to60_pr_symm->Clone();
  h_40to60_symm->Fit(function4, "", "", -0.5, 0.5);
  std::cout << "Chi^2/NDF = " << function4->GetChisquare() << "/" << function4->GetNDF() << " = " << function4->GetChisquare() / function4->GetNDF() << std::endl;
  
  
  // Make mirrored plots
  for (int i = 1; i <= 20; i++)
    {
      int j = 0;  // mirrored bin
      
      switch(i)
	{
	case 11: j = 10; break;
	case 12: j = 9; break;
	case 13: j = 8; break;
	case 14: j = 7; break;
	case 15: j = 6; break;
	case 16: j = 5; break;
	case 17: j = 4; break;
	case 18: j = 3; break;
	case 19: j = 2; break;
	case 20: j = 1; break;
	}
      
      if(j != 0)
	{
	  h_vn_yCM_00to10_pr_mirror->SetBinContent(j, -h_vn_yCM_00to10_pr->GetBinContent(i));
	  h_vn_yCM_00to10_pr_mirror->SetBinError(j, h_vn_yCM_00to10_pr->GetBinError(i));
	  h_vn_yCM_10to40_pr_mirror->SetBinContent(j, -h_vn_yCM_10to40_pr->GetBinContent(i));
	  h_vn_yCM_10to40_pr_mirror->SetBinError(j, h_vn_yCM_10to40_pr->GetBinError(i));
	  h_vn_yCM_40to60_pr_mirror->SetBinContent(j, -h_vn_yCM_40to60_pr->GetBinContent(i));
	  h_vn_yCM_40to60_pr_mirror->SetBinError(j, h_vn_yCM_40to60_pr->GetBinError(i));
	}
    }

  if (useSystematics)
    {
      int j = 19;
      for (int i = 0; i < 10; i++)
	{
	  sys_yCM_00to10_pr->SetPointY(i, -sys_yCM_00to10_pr->GetPointY(j));
	  sys_yCM_00to10_pr->SetPointError(i,0.0, sys_yCM_00to10_pr->GetErrorY(j));
	  sys_yCM_10to40_pr->SetPointY(i, -sys_yCM_10to40_pr->GetPointY(j));
	  sys_yCM_10to40_pr->SetPointError(i,0.0, sys_yCM_10to40_pr->GetErrorY(j));
	  sys_yCM_40to60_pr->SetPointY(i, -sys_yCM_40to60_pr->GetPointY(j));
	  sys_yCM_40to60_pr->SetPointError(i,0.0, sys_yCM_40to60_pr->GetErrorY(j));
	  j--;
	}
    }

  
  // Set various aesthetics
  h_vn_yCM_00to10_pr->SetMarkerStyle(21);
  h_vn_yCM_10to40_pr->SetMarkerStyle(20);
  h_vn_yCM_40to60_pr->SetMarkerStyle(33);
  h_vn_yCM_00to10_pr->SetMarkerColor(kRed-4);
  //h_vn_yCM_10to40_pr->SetMarkerColor(kBlue-4);
  h_vn_yCM_40to60_pr->SetMarkerColor(kBlue-4);
  h_vn_yCM_00to10_pr->SetMarkerSize(2);
  h_vn_yCM_10to40_pr->SetMarkerSize(2);
  h_vn_yCM_40to60_pr->SetMarkerSize(3);
  h_vn_yCM_00to10_pr->SetLineColor(kRed-4);
  //h_vn_yCM_10to40_pr->SetLineColor(kBlue-4);
  h_vn_yCM_40to60_pr->SetLineColor(kBlue-4);
  h_vn_yCM_00to10_pr->SetLineWidth(3);
  h_vn_yCM_10to40_pr->SetLineWidth(3);
  h_vn_yCM_40to60_pr->SetLineWidth(3);

  h_vn_yCM_00to10_pr_symm->SetMarkerStyle(21);
  h_vn_yCM_10to40_pr_symm->SetMarkerStyle(20);
  h_vn_yCM_40to60_pr_symm->SetMarkerStyle(33);
  h_vn_yCM_00to10_pr_symm->SetMarkerColor(kRed-4);
  //h_vn_yCM_10to40_pr_symm->SetMarkerColor(kBlue-4);
  h_vn_yCM_40to60_pr_symm->SetMarkerColor(kBlue-4);
  h_vn_yCM_00to10_pr_symm->SetMarkerSize(2);
  h_vn_yCM_10to40_pr_symm->SetMarkerSize(2);
  h_vn_yCM_40to60_pr_symm->SetMarkerSize(3);
  h_vn_yCM_00to10_pr_symm->SetLineColor(kRed-4);
  //h_vn_yCM_10to40_pr_symm->SetLineColor(kBlue-4);
  h_vn_yCM_40to60_pr_symm->SetLineColor(kBlue-4);
  h_vn_yCM_00to10_pr_symm->SetLineWidth(3);
  h_vn_yCM_10to40_pr_symm->SetLineWidth(3);
  h_vn_yCM_40to60_pr_symm->SetLineWidth(3);

  if (useSystematics)
    {
      sys_yCM_00to10_pr->SetMarkerColor(kRed-4);
      //sys_yCM_10to40_pr->SetMarkerColor(kBlue-4);
      sys_yCM_40to60_pr->SetMarkerColor(kBlue-4);
      sys_yCM_00to10_pr->SetLineColor(kRed-4);
      //sys_yCM_10to40_pr->SetLineColor(kBlue-4);
      sys_yCM_40to60_pr->SetLineColor(kBlue-4);

      sys_yCM_00to10_pr_symm->SetMarkerColor(kRed-4);
      //sys_yCM_10to40_pr_symm->SetMarkerColor(kBlue-4);
      sys_yCM_40to60_pr_symm->SetMarkerColor(kBlue-4);
      sys_yCM_00to10_pr_symm->SetLineColor(kRed-4);
      //sys_yCM_10to40_pr_symm->SetLineColor(kBlue-4);
      sys_yCM_40to60_pr_symm->SetLineColor(kBlue-4);
    }
  
  h_vn_yCM_00to10_pr_mirror->SetMarkerStyle(25);
  h_vn_yCM_10to40_pr_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_pr_mirror->SetMarkerStyle(27);
  h_vn_yCM_00to10_pr_mirror->SetMarkerColor(kRed-4);
  //h_vn_yCM_10to40_pr_mirror->SetMarkerColor(kBlue-4);
  h_vn_yCM_40to60_pr_mirror->SetMarkerColor(kBlue-4);
  h_vn_yCM_00to10_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_pr_mirror->SetMarkerSize(3);
  h_vn_yCM_00to10_pr_mirror->SetLineColor(kRed-4);
  //h_vn_yCM_10to40_pr_mirror->SetLineColor(kBlue-4);
  h_vn_yCM_40to60_pr_mirror->SetLineColor(kBlue-4);
  h_vn_yCM_00to10_pr_mirror->SetLineWidth(3);
  h_vn_yCM_10to40_pr_mirror->SetLineWidth(3);
  h_vn_yCM_40to60_pr_mirror->SetLineWidth(3);
  ////

  // Even and odd components of the symmetric plots
  /*
  TH1D* h_even_00to10 = PlotUtils::getEvenComponent(h_vn_yCM_00to10_pr_symm);
  TH1D* h_even_10to40 = PlotUtils::getEvenComponent(h_vn_yCM_10to40_pr_symm);
  TH1D* h_even_40to60 = PlotUtils::getEvenComponent(h_vn_yCM_40to60_pr_symm);
  TH1D* h_odd_00to10 = PlotUtils::getOddComponent(h_vn_yCM_00to10_pr_symm);
  TH1D* h_odd_10to40 = PlotUtils::getOddComponent(h_vn_yCM_10to40_pr_symm);
  TH1D* h_odd_40to60 = PlotUtils::getOddComponent(h_vn_yCM_40to60_pr_symm);
  */
  TProfileHelper ph_even_00to10;
  TProfileHelper ph_even_10to40;
  TProfileHelper ph_even_40to60;
  TProfileHelper ph_odd_00to10;
  TProfileHelper ph_odd_10to40;
  TProfileHelper ph_odd_40to60;

  ph_even_00to10.initialize(p_vn_yCM_00to10_pr_symm);
  ph_even_10to40.initialize(p_vn_yCM_10to40_pr_symm);
  ph_even_40to60.initialize(p_vn_yCM_40to60_pr_symm);
  ph_odd_00to10.initialize(p_vn_yCM_00to10_pr_symm);
  ph_odd_10to40.initialize(p_vn_yCM_10to40_pr_symm);
  ph_odd_40to60.initialize(p_vn_yCM_40to60_pr_symm);

  ph_even_00to10.removeOddComponent();
  ph_even_10to40.removeOddComponent();
  ph_even_40to60.removeOddComponent();
  ph_odd_00to10.removeEvenComponent();
  ph_odd_10to40.removeEvenComponent();
  ph_odd_40to60.removeEvenComponent();

  TH1D* h_even_00to10 = ph_even_00to10.getHistogram("h_even_00to10");
  TH1D* h_even_10to40 = ph_even_10to40.getHistogram("h_even_10to40");
  TH1D* h_even_40to60 = ph_even_40to60.getHistogram("h_even_40to60");
  TH1D* h_odd_00to10 = ph_odd_00to10.getHistogram("h_odd_00to10");
  TH1D* h_odd_10to40 = ph_odd_10to40.getHistogram("h_odd_10to40");
  TH1D* h_odd_40to60 = ph_odd_40to60.getHistogram("h_odd_40to60");

  h_even_00to10->SetMarkerStyle(21);
  h_even_10to40->SetMarkerStyle(20);
  h_even_40to60->SetMarkerStyle(33);
  h_even_00to10->SetMarkerColor(kRed-4);
  //h_even_10to40->SetMarkerColor(kBlue-4);
  h_even_40to60->SetMarkerColor(kBlue-4);
  h_even_00to10->SetMarkerSize(2);
  h_even_10to40->SetMarkerSize(2);
  h_even_40to60->SetMarkerSize(3);
  h_even_00to10->SetLineColor(kRed-4);
  //h_even_10to40->SetLineColor(kBlue-4);
  h_even_40to60->SetLineColor(kBlue-4);
  h_even_00to10->SetLineWidth(3);
  h_even_10to40->SetLineWidth(3);
  h_even_40to60->SetLineWidth(3);

  h_odd_00to10->SetMarkerStyle(21);
  h_odd_10to40->SetMarkerStyle(20);
  h_odd_40to60->SetMarkerStyle(33);
  h_odd_00to10->SetMarkerColor(kRed-4);
  //h_odd_10to40->SetMarkerColor(kBlue-4);
  h_odd_40to60->SetMarkerColor(kBlue-4);
  h_odd_00to10->SetMarkerSize(2);
  h_odd_10to40->SetMarkerSize(2);
  h_odd_40to60->SetMarkerSize(3);
  h_odd_00to10->SetLineColor(kRed-4);
  //h_odd_10to40->SetLineColor(kBlue-4);
  h_odd_40to60->SetLineColor(kBlue-4);
  h_odd_00to10->SetLineWidth(3);
  h_odd_10to40->SetLineWidth(3);
  h_odd_40to60->SetLineWidth(3);
  ////

  THStack *prRapidityStack   = new THStack("prRapidityStack", ";y-y_{mid};v_{1} {#Psi_{1}}");
  prRapidityStack->Add(h_vn_yCM_00to10_pr);
  if (sqrt_s_NN < 3.5)
    prRapidityStack->Add(h_vn_yCM_40to60_pr);
  prRapidityStack->Add(h_vn_yCM_10to40_pr);
  prRapidityStack->Add(h_vn_yCM_00to10_pr_mirror);
  if (sqrt_s_NN < 3.5)
    prRapidityStack->Add(h_vn_yCM_40to60_pr_mirror);
  prRapidityStack->Add(h_vn_yCM_10to40_pr_mirror);

  THStack *prRapidityStack_symm = new THStack("prRapidityStack_symm", ";y-y_{mid};v_{1} {#Psi_{1}}");
  prRapidityStack_symm->Add(h_vn_yCM_00to10_pr_symm);
  if (sqrt_s_NN < 3.5)
    prRapidityStack_symm->Add(h_vn_yCM_40to60_pr_symm);
  prRapidityStack_symm->Add(h_vn_yCM_10to40_pr_symm);

  THStack *prEvenStack = new THStack("prEvenStack", ";y-y_{mid};v_{1}^{even} {#Psi_{1}}");
  if (sqrt_s_NN < 3.5)
    prEvenStack->Add(h_even_40to60);
  prEvenStack->Add(h_even_00to10);
  prEvenStack->Add(h_even_10to40);
  
  THStack *prOddStack = new THStack("prOddStack", ";y-y_{mid};v_{1}^{odd} {#Psi_{1}}");
  if (sqrt_s_NN < 3.5)
    prOddStack->Add(h_odd_40to60);
  prOddStack->Add(h_odd_00to10);
  prOddStack->Add(h_odd_10to40);
  

  // Make text boxes, legends, and line at zero
  TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
  prLegend->AddEntry(h_vn_yCM_00to10_pr, "0 - 10%");
  prLegend->AddEntry(h_vn_yCM_10to40_pr, "10 - 40%");
  if (sqrt_s_NN < 3.5)
    prLegend->AddEntry(h_vn_yCM_40to60_pr, "40 - 60%");
  prLegend->SetBorderSize(0);
  prLegend->SetFillColorAlpha(0,0);
  prLegend->SetTextSize(0.04);

  TLegend *prLegend_symm = new TLegend(0.19, 0.15, 0.39, 0.3);
  prLegend_symm->AddEntry(h_vn_yCM_00to10_pr_symm, "0 - 10%");
  prLegend_symm->AddEntry(h_vn_yCM_10to40_pr_symm, "10 - 40%");
  if (sqrt_s_NN < 3.5)
    prLegend_symm->AddEntry(h_vn_yCM_40to60_pr_symm, "40 - 60%");
  prLegend_symm->SetBorderSize(0);
  prLegend_symm->SetFillColorAlpha(0,0);
  prLegend_symm->SetTextSize(0.04);

  TPaveText *prText_y = new TPaveText(-0.35, 0.03, 0.75, 0.055, "NB");
  prText_y->AddText("Proton");
  if (sqrt_s_NN == 3.0)
    prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT (year 2018)");
  else if (sqrt_s_NN == 3.2 || sqrt_s_NN == 3.22)
    prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.2 GeV FXT (year 2019)");
  else if (sqrt_s_NN == 3.5)
    prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.5 GeV FXT (year 2020)");
  else if (sqrt_s_NN == 3.9)
    prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.9 GeV FXT (years 2019, 2020)");
  else if (sqrt_s_NN == 4.5)
    prText_y->AddText("Au+Au #sqrt{s_{NN}} = 4.5 GeV FXT (year 2020)");
  prText_y->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
  prText_y->SetFillColorAlpha(0,0);
  prText_y->SetLineColorAlpha(0,0);
  prText_y->SetTextSize(.035);
 
  //TPaveText *prText_y_symm = new TPaveText(-0.6, 0.037, 0.6, 0.075, "NB");
  TPaveText *prText_y_symm = new TPaveText(-0.8, 0.03, 1.0, 0.065, "NB");
  prText_y_symm->SetTextAlign(32);
  if (sqrt_s_NN == 3.0)
    prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT (year 2018)");
  else if (sqrt_s_NN == 3.2 || sqrt_s_NN == 3.22)
    prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 3.2 GeV FXT (year 2019)");
  else if (sqrt_s_NN == 3.5)
    prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 3.5 GeV FXT (year 2020)");
  else if (sqrt_s_NN == 3.9)
    prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 3.9 GeV FXT (years 2019, 2020)");
  else if (sqrt_s_NN == 4.5)
    prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 4.5 GeV FXT (year 2020)");
  prText_y_symm->AddText("Proton");
  prText_y_symm->AddText("1.0 #leq p_{T} #leq 2.5 GeV");
  prText_y_symm->SetFillColorAlpha(0,0);
  prText_y_symm->SetLineColorAlpha(0,0);
  prText_y_symm->SetTextSize(.045);

  TPaveText *prText_odd = new TPaveText(-0.8, 0.04, 1.0, 0.075, "NB");
  prText_odd->SetTextAlign(32);
  if (sqrt_s_NN == 3.0)
    prText_odd->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT (year 2018)");
  else if (sqrt_s_NN == 3.2 || sqrt_s_NN == 3.22)
    prText_odd->AddText("Au+Au #sqrt{s_{NN}} = 3.2 GeV FXT (year 2019)");
  else if (sqrt_s_NN == 3.5)
    prText_odd->AddText("Au+Au #sqrt{s_{NN}} = 3.5 GeV FXT (year 2020)");
  else if (sqrt_s_NN == 3.9)
    prText_odd->AddText("Au+Au #sqrt{s_{NN}} = 3.9 GeV FXT (years 2019, 2020)");
  else if (sqrt_s_NN == 4.5)
    prText_odd->AddText("Au+Au #sqrt{s_{NN}} = 4.5 GeV FXT (year 2020)");
  prText_odd->AddText("Proton");
  prText_odd->AddText("1.0 #leq p_{T} #leq 2.5 GeV");
  prText_odd->SetFillColorAlpha(0,0);
  prText_odd->SetLineColorAlpha(0,0);
  prText_odd->SetTextSize(.04);

  TPaveText* prelimText_symm = new TPaveText(-0.9, -0.05, -0.1, -0.04, "NB");
  prelimText_symm->AddText("STAR Preliminary");
  prelimText_symm->SetTextColor(kRed);
  prelimText_symm->SetFillColorAlpha(0,0);
  prelimText_symm->SetTextSize(0.04);

  TPaveText* prelimText = new TPaveText(-0.2, 0.04, 0.6, 0.045, "NB");
  prelimText->AddText("STAR Preliminary");
  prelimText->SetTextColor(kRed);
  prelimText->SetFillColorAlpha(0,0);
  prelimText->SetTextSize(0.04);

  TLine *zeroLine_y_pr = new TLine(-1, 0, 1, 0);
  zeroLine_y_pr->SetLineStyle(9);
  zeroLine_y_pr->SetLineWidth(3);
  ////

  prRapidityStack->Draw();
  prRapidityStack->GetYaxis()->SetLabelSize(0.043);
  prRapidityStack->GetXaxis()->SetLabelSize(0.043);
  prRapidityStack->GetYaxis()->SetTitleOffset(1.4);
  prRapidityStack->GetXaxis()->SetTitleOffset(1.0);
  prRapidityStack->GetXaxis()->SetNdivisions(210);
  prRapidityStack->GetXaxis()->SetTitleSize(0.045);
  prRapidityStack->GetYaxis()->SetTitleSize(0.05);
  prRapidityStack->Draw();
  prRapidityStack->SetMaximum(0.06);
  prRapidityStack->SetMinimum(-0.06);
  prRapidityStack->Draw("NOSTACK EP");
  zeroLine_y_pr->Draw("SAME");
  if (useSystematics)
    {
      sys_yCM_40to60_pr->Draw("[]");
      sys_yCM_00to10_pr->Draw("[]");
      sys_yCM_10to40_pr->Draw("[]");
    }
  prRapidityStack->Draw("NOSTACK EP SAME");
  prLegend->Draw();
  //prText_y->Draw();

  //prelimText->Draw();
  //function->Draw("C SAME");
  canvas->SaveAs("v3_prRapidityStack.pdf");
  //canvas->SaveAs("v3_prRapidityStack.png");
  canvas->Clear();


  prRapidityStack_symm->Draw();
  prRapidityStack_symm->GetYaxis()->SetLabelSize(0.043);
  prRapidityStack_symm->GetXaxis()->SetLabelSize(0.043);
  prRapidityStack_symm->GetYaxis()->SetTitleOffset(1.4);
  prRapidityStack_symm->GetXaxis()->SetTitleOffset(1.0);
  prRapidityStack_symm->GetXaxis()->SetNdivisions(210);
  prRapidityStack_symm->GetXaxis()->SetTitleSize(0.045);
  prRapidityStack_symm->GetYaxis()->SetTitleSize(0.05);
  prRapidityStack_symm->Draw();
  prRapidityStack_symm->SetMaximum(0.07);
  prRapidityStack_symm->SetMinimum(-0.09);
  prRapidityStack_symm->Draw("NOSTACK EP");
  zeroLine_y_pr->Draw("SAME");
  //h_00to10_symm->Draw("SAME");
  //h_10to40_symm->Draw("SAME");
  //h_40to60_symm->Draw("SAME");
  if (useSystematics)
    {
      sys_yCM_00to10_pr_symm->Draw("[]");
      sys_yCM_40to60_pr_symm->Draw("[]");
      sys_yCM_10to40_pr_symm->Draw("[]");
    }
  prRapidityStack_symm->Draw("NOSTACK EP SAME");
  prLegend_symm->Draw();
  //prText_y_symm->Draw();

  //prelimText_symm->Draw();
  //function1->Draw("C SAME");
  //function2->Draw("C SAME");
  //function3->Draw("C SAME");
  canvas->SaveAs("v3_prRapidityStack_symm.pdf");
  //canvas->SaveAs("v3_prRapidityStack_symm.png");
  canvas->Clear();

  prEvenStack->Draw();
  prEvenStack->GetYaxis()->SetLabelSize(0.043);
  prEvenStack->GetXaxis()->SetLabelSize(0.043);
  prEvenStack->GetYaxis()->SetTitleOffset(1.4);
  prEvenStack->GetXaxis()->SetTitleOffset(1.0);
  prEvenStack->GetXaxis()->SetNdivisions(210);
  prEvenStack->GetXaxis()->SetTitleSize(0.045);
  prEvenStack->GetYaxis()->SetTitleSize(0.05);
  prEvenStack->Draw();
  prEvenStack->SetMaximum(0.02);
  prEvenStack->SetMinimum(-0.03);
  prEvenStack->Draw("NOSTACK EP");
  zeroLine_y_pr->Draw("SAME");
  prEvenStack->Draw("NOSTACK EP SAME");
  prLegend_symm->Draw();
  //prText_y_symm->Draw();
  //prelimText_symm->Draw();
  canvas->SaveAs("v3_prEvenStack.pdf");
  //canvas->SaveAs("v3_prEvenStack.png");
  canvas->Clear();

  prOddStack->Draw();
  prOddStack->GetYaxis()->SetLabelSize(0.043);
  prOddStack->GetXaxis()->SetLabelSize(0.043);
  prOddStack->GetYaxis()->SetTitleOffset(1.4);
  prOddStack->GetXaxis()->SetTitleOffset(1.0);
  prOddStack->GetXaxis()->SetNdivisions(210);
  prOddStack->GetXaxis()->SetTitleSize(0.045);
  prOddStack->GetYaxis()->SetTitleSize(0.05);
  prOddStack->Draw();
  prOddStack->SetMaximum(0.08);
  prOddStack->SetMinimum(-0.08);
  prOddStack->Draw("NOSTACK EP");
  zeroLine_y_pr->Draw("SAME");
  prOddStack->Draw("NOSTACK EP SAME");
  prLegend_symm->Draw();
  //prText_odd->Draw();
  //prelimText_symm->Draw();
  canvas->SaveAs("v3_prOddStack.pdf");
  //canvas->SaveAs("v3_prOddStack.png");
  canvas->Clear();
  
  if (useSystematics)
    systematicFile->Close();
  
  file->Close();
}
