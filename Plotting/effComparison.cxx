#include <iostream>

void effComparison()
{
  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 800);
  canvas->SetTicks();
  
  TFile *normalFile = TFile::Open("Normal.picoDst.result.combined.root", "READ");
  TFile *efficiencyFile = TFile::Open("Efficiency.picoDst.result.combined.root", "READ");

  // Normal profiles
  TProfile *p_vn_pp_norm = (TProfile*)normalFile->Get("p_vn_pp");
  TProfile *p_vn_pm_norm = (TProfile*)normalFile->Get("p_vn_pm");
  TProfile *p_vn_kp_norm = (TProfile*)normalFile->Get("p_vn_kp");
  TProfile *p_vn_km_norm = (TProfile*)normalFile->Get("p_vn_km");
  TProfile *p_vn_pr_norm = (TProfile*)normalFile->Get("p_vn_pr");
  p_vn_kp_norm->Rebin();
  p_vn_km_norm->Rebin();

  TProfile *p_vn_pp_ext_norm = (TProfile*)normalFile->Get("p_vn_pp_ext");
  TProfile *p_vn_pm_ext_norm = (TProfile*)normalFile->Get("p_vn_pm_ext");
  TProfile *p_vn_kp_ext_norm = (TProfile*)normalFile->Get("p_vn_kp_ext");
  TProfile *p_vn_km_ext_norm = (TProfile*)normalFile->Get("p_vn_km_ext");
  TProfile *p_vn_pr_ext_norm = (TProfile*)normalFile->Get("p_vn_pr_ext");

  TProfile *p_vn_pr_for_norm = (TProfile*)normalFile->Get("p_vn_pr_for");

  p_vn_kp_ext_norm->Rebin();
  p_vn_km_ext_norm->Rebin();
  
  TProfile2D *p2_vn_yCM_cent_pp_norm = (TProfile2D*)normalFile->Get("p2_vn_yCM_cent_pp");
  TProfile *p_vn_yCM_00to10_pp_norm = p2_vn_yCM_cent_pp_norm->ProfileY("p_vn_yCM_00to10_pp_norm", 15, 16);
  TProfile *p_vn_yCM_10to40_pp_norm = p2_vn_yCM_cent_pp_norm->ProfileY("p_vn_yCM_10to40_pp_norm", 9, 14);
  TProfile *p_vn_yCM_40to60_pp_norm = p2_vn_yCM_cent_pp_norm->ProfileY("p_vn_yCM_40to60_pp_norm", 5, 8);

  TProfile2D *p2_vn_yCM_cent_pm_norm = (TProfile2D*)normalFile->Get("p2_vn_yCM_cent_pm");
  TProfile *p_vn_yCM_00to10_pm_norm = p2_vn_yCM_cent_pm_norm->ProfileY("p_vn_yCM_00to10_pm_norm", 15, 16);
  TProfile *p_vn_yCM_10to40_pm_norm = p2_vn_yCM_cent_pm_norm->ProfileY("p_vn_yCM_10to40_pm_norm", 9, 14);
  TProfile *p_vn_yCM_40to60_pm_norm = p2_vn_yCM_cent_pm_norm->ProfileY("p_vn_yCM_40to60_pm_norm", 5, 8);

  TProfile2D *p2_vn_yCM_cent_kp_norm = (TProfile2D*)normalFile->Get("p2_vn_yCM_cent_kp");
  p2_vn_yCM_cent_kp_norm->RebinY();
  TProfile *p_vn_yCM_00to10_kp_norm = p2_vn_yCM_cent_kp_norm->ProfileY("p_vn_yCM_00to10_kp_norm", 15, 16);
  TProfile *p_vn_yCM_10to40_kp_norm = p2_vn_yCM_cent_kp_norm->ProfileY("p_vn_yCM_10to40_kp_norm", 9, 14);
  TProfile *p_vn_yCM_40to60_kp_norm = p2_vn_yCM_cent_kp_norm->ProfileY("p_vn_yCM_40to60_kp_norm", 5, 8);

  TProfile2D *p2_vn_yCM_cent_km_norm = (TProfile2D*)normalFile->Get("p2_vn_yCM_cent_km");
  p2_vn_yCM_cent_km_norm->RebinY();
  TProfile *p_vn_yCM_00to10_km_norm = p2_vn_yCM_cent_km_norm->ProfileY("p_vn_yCM_00to10_km_norm", 15, 16);
  TProfile *p_vn_yCM_10to40_km_norm = p2_vn_yCM_cent_km_norm->ProfileY("p_vn_yCM_10to40_km_norm", 9, 14);
  TProfile *p_vn_yCM_40to60_km_norm = p2_vn_yCM_cent_km_norm->ProfileY("p_vn_yCM_40to60_km_norm", 5, 8);

  TProfile2D *p2_vn_yCM_cent_pr_norm = (TProfile2D*)normalFile->Get("p2_vn_yCM_cent_pr");
  TProfile *p_vn_yCM_00to10_pr_norm = p2_vn_yCM_cent_pr_norm->ProfileY("p_vn_yCM_00to10_pr_norm", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_norm = p2_vn_yCM_cent_pr_norm->ProfileY("p_vn_yCM_10to40_pr_norm", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_norm = p2_vn_yCM_cent_pr_norm->ProfileY("p_vn_yCM_40to60_pr_norm", 5, 8);

  TProfile2D *p2_vn_yCM_cent_pr_norm_symmetry = (TProfile2D*)normalFile->Get("p2_vn_yCM_cent_pr_symmetry");
  TProfile *p_vn_yCM_00to10_pr_norm_symm = p2_vn_yCM_cent_pr_norm_symmetry->ProfileY("p_vn_yCM_00to10_pr_norm_symm", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_norm_symm = p2_vn_yCM_cent_pr_norm_symmetry->ProfileY("p_vn_yCM_10to40_pr_norm_symm", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_norm_symm = p2_vn_yCM_cent_pr_norm_symmetry->ProfileY("p_vn_yCM_40to60_pr_norm_symm", 5, 8);


  TProfile2D *p2_vn_pT_cent_pp_norm = (TProfile2D*)normalFile->Get("p2_vn_pT_cent_pp");
  TProfile *p_vn_pT_00to10_pp_norm = p2_vn_pT_cent_pp_norm->ProfileY("p_vn_pT_00to10_pp_norm", 15, 16);
  TProfile *p_vn_pT_10to40_pp_norm = p2_vn_pT_cent_pp_norm->ProfileY("p_vn_pT_10to40_pp_norm", 9, 14);
  TProfile *p_vn_pT_40to60_pp_norm = p2_vn_pT_cent_pp_norm->ProfileY("p_vn_pT_40to60_pp_norm", 5, 8);
  
  TProfile2D *p2_vn_pT_cent_pm_norm = (TProfile2D*)normalFile->Get("p2_vn_pT_cent_pm");
  TProfile *p_vn_pT_00to10_pm_norm = p2_vn_pT_cent_pm_norm->ProfileY("p_vn_pT_00to10_pm_norm", 15, 16);
  TProfile *p_vn_pT_10to40_pm_norm = p2_vn_pT_cent_pm_norm->ProfileY("p_vn_pT_10to40_pm_norm", 9, 14);
  TProfile *p_vn_pT_40to60_pm_norm = p2_vn_pT_cent_pm_norm->ProfileY("p_vn_pT_40to60_pm_norm", 5, 8);
  
  TProfile2D *p2_vn_pT_cent_kp_norm = (TProfile2D*)normalFile->Get("p2_vn_pT_cent_kp");
  p2_vn_pT_cent_kp_norm->RebinY();
  TProfile *p_vn_pT_00to10_kp_norm = p2_vn_pT_cent_kp_norm->ProfileY("p_vn_pT_00to10_kp_norm", 15, 16);
  TProfile *p_vn_pT_10to40_kp_norm = p2_vn_pT_cent_kp_norm->ProfileY("p_vn_pT_10to40_kp_norm", 9, 14);
  TProfile *p_vn_pT_40to60_kp_norm = p2_vn_pT_cent_kp_norm->ProfileY("p_vn_pT_40to60_kp_norm", 5, 8);

  TProfile2D *p2_vn_pT_cent_km_norm = (TProfile2D*)normalFile->Get("p2_vn_pT_cent_km");
  p2_vn_pT_cent_km_norm->RebinY();
  TProfile *p_vn_pT_00to10_km_norm = p2_vn_pT_cent_km_norm->ProfileY("p_vn_pT_00to10_km_norm", 15, 16);
  TProfile *p_vn_pT_10to40_km_norm = p2_vn_pT_cent_km_norm->ProfileY("p_vn_pT_10to40_km_norm", 9, 14);
  TProfile *p_vn_pT_40to60_km_norm = p2_vn_pT_cent_km_norm->ProfileY("p_vn_pT_40to60_km_norm", 5, 8);

  TProfile2D *p2_vn_pT_cent_pr_norm = (TProfile2D*)normalFile->Get("p2_vn_pT_cent_pr");
  TProfile *p_vn_pT_00to10_pr_norm = p2_vn_pT_cent_pr_norm->ProfileY("p_vn_pT_00to10_pr_norm", 15, 16);
  TProfile *p_vn_pT_10to40_pr_norm = p2_vn_pT_cent_pr_norm->ProfileY("p_vn_pT_10to40_pr_norm", 9, 14);
  TProfile *p_vn_pT_40to60_pr_norm = p2_vn_pT_cent_pr_norm->ProfileY("p_vn_pT_40to60_pr_norm", 5, 8);

  
  
  // Efficiency profiles
  TProfile *p_vn_pp_eff = (TProfile*)efficiencyFile->Get("p_vn_pp");
  TProfile *p_vn_pm_eff = (TProfile*)efficiencyFile->Get("p_vn_pm");
  TProfile *p_vn_kp_eff = (TProfile*)efficiencyFile->Get("p_vn_kp");
  TProfile *p_vn_km_eff = (TProfile*)efficiencyFile->Get("p_vn_km");
  TProfile *p_vn_pr_eff = (TProfile*)efficiencyFile->Get("p_vn_pr");
  p_vn_kp_eff->Rebin();
  p_vn_km_eff->Rebin();

  TProfile *p_vn_pp_ext_eff = (TProfile*)efficiencyFile->Get("p_vn_pp_ext");
  TProfile *p_vn_pm_ext_eff = (TProfile*)efficiencyFile->Get("p_vn_pm_ext");
  TProfile *p_vn_kp_ext_eff = (TProfile*)efficiencyFile->Get("p_vn_kp_ext");
  TProfile *p_vn_km_ext_eff = (TProfile*)efficiencyFile->Get("p_vn_km_ext");
  TProfile *p_vn_pr_ext_eff = (TProfile*)efficiencyFile->Get("p_vn_pr_ext");

  TProfile *p_vn_pr_for_eff = (TProfile*)efficiencyFile->Get("p_vn_pr_for");

  p_vn_kp_ext_eff->Rebin();
  p_vn_km_ext_eff->Rebin();

  TProfile2D *p2_vn_yCM_cent_pp_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_yCM_cent_pp");
  TProfile *p_vn_yCM_00to10_pp_eff = p2_vn_yCM_cent_pp_eff->ProfileY("p_vn_yCM_00to10_pp_eff", 15, 16);
  TProfile *p_vn_yCM_10to40_pp_eff = p2_vn_yCM_cent_pp_eff->ProfileY("p_vn_yCM_10to40_pp_eff", 9, 14);
  TProfile *p_vn_yCM_40to60_pp_eff = p2_vn_yCM_cent_pp_eff->ProfileY("p_vn_yCM_40to60_pp_eff", 5, 8);

  TProfile2D *p2_vn_yCM_cent_pm_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_yCM_cent_pm");
  TProfile *p_vn_yCM_00to10_pm_eff = p2_vn_yCM_cent_pm_eff->ProfileY("p_vn_yCM_00to10_pm_eff", 15, 16);
  TProfile *p_vn_yCM_10to40_pm_eff = p2_vn_yCM_cent_pm_eff->ProfileY("p_vn_yCM_10to40_pm_eff", 9, 14);
  TProfile *p_vn_yCM_40to60_pm_eff = p2_vn_yCM_cent_pm_eff->ProfileY("p_vn_yCM_40to60_pm_eff", 5, 8);

  TProfile2D *p2_vn_yCM_cent_kp_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_yCM_cent_kp");
  p2_vn_yCM_cent_kp_eff->RebinY();
  TProfile *p_vn_yCM_00to10_kp_eff = p2_vn_yCM_cent_kp_eff->ProfileY("p_vn_yCM_00to10_kp_eff", 15, 16);
  TProfile *p_vn_yCM_10to40_kp_eff = p2_vn_yCM_cent_kp_eff->ProfileY("p_vn_yCM_10to40_kp_eff", 9, 14);
  TProfile *p_vn_yCM_40to60_kp_eff = p2_vn_yCM_cent_kp_eff->ProfileY("p_vn_yCM_40to60_kp_eff", 5, 8);

  TProfile2D *p2_vn_yCM_cent_km_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_yCM_cent_km");
  p2_vn_yCM_cent_km_eff->RebinY();
  TProfile *p_vn_yCM_00to10_km_eff = p2_vn_yCM_cent_km_eff->ProfileY("p_vn_yCM_00to10_km_eff", 15, 16);
  TProfile *p_vn_yCM_10to40_km_eff = p2_vn_yCM_cent_km_eff->ProfileY("p_vn_yCM_10to40_km_eff", 9, 14);
  TProfile *p_vn_yCM_40to60_km_eff = p2_vn_yCM_cent_km_eff->ProfileY("p_vn_yCM_40to60_km_eff", 5, 8);
  
  TProfile2D *p2_vn_yCM_cent_pr_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_yCM_cent_pr");
  TProfile *p_vn_yCM_00to10_pr_eff = p2_vn_yCM_cent_pr_eff->ProfileY("p_vn_yCM_00to10_pr_eff", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_eff = p2_vn_yCM_cent_pr_eff->ProfileY("p_vn_yCM_10to40_pr_eff", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_eff = p2_vn_yCM_cent_pr_eff->ProfileY("p_vn_yCM_40to60_pr_eff", 5, 8);

  TProfile2D *p2_vn_yCM_cent_pr_eff_symmetry = (TProfile2D*)efficiencyFile->Get("p2_vn_yCM_cent_pr_symmetry");
  TProfile *p_vn_yCM_00to10_pr_eff_symm = p2_vn_yCM_cent_pr_eff_symmetry->ProfileY("p_vn_yCM_00to10_pr_eff_symm", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_eff_symm = p2_vn_yCM_cent_pr_eff_symmetry->ProfileY("p_vn_yCM_10to40_pr_eff_symm", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_eff_symm = p2_vn_yCM_cent_pr_eff_symmetry->ProfileY("p_vn_yCM_40to60_pr_eff_symm", 5, 8);


  TProfile2D *p2_vn_pT_cent_pp_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_pT_cent_pp");
  TProfile *p_vn_pT_00to10_pp_eff = p2_vn_pT_cent_pp_eff->ProfileY("p_vn_pT_00to10_pp_eff", 15, 16);
  TProfile *p_vn_pT_10to40_pp_eff = p2_vn_pT_cent_pp_eff->ProfileY("p_vn_pT_10to40_pp_eff", 9, 14);
  TProfile *p_vn_pT_40to60_pp_eff = p2_vn_pT_cent_pp_eff->ProfileY("p_vn_pT_40to60_pp_eff", 5, 8);
  
  TProfile2D *p2_vn_pT_cent_pm_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_pT_cent_pm");
  TProfile *p_vn_pT_00to10_pm_eff = p2_vn_pT_cent_pm_eff->ProfileY("p_vn_pT_00to10_pm_eff", 15, 16);
  TProfile *p_vn_pT_10to40_pm_eff = p2_vn_pT_cent_pm_eff->ProfileY("p_vn_pT_10to40_pm_eff", 9, 14);
  TProfile *p_vn_pT_40to60_pm_eff = p2_vn_pT_cent_pm_eff->ProfileY("p_vn_pT_40to60_pm_eff", 5, 8);
  
  TProfile2D *p2_vn_pT_cent_kp_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_pT_cent_kp");
  p2_vn_pT_cent_kp_eff->RebinY();
  TProfile *p_vn_pT_00to10_kp_eff = p2_vn_pT_cent_kp_eff->ProfileY("p_vn_pT_00to10_kp_eff", 15, 16);
  TProfile *p_vn_pT_10to40_kp_eff = p2_vn_pT_cent_kp_eff->ProfileY("p_vn_pT_10to40_kp_eff", 9, 14);
  TProfile *p_vn_pT_40to60_kp_eff = p2_vn_pT_cent_kp_eff->ProfileY("p_vn_pT_40to60_kp_eff", 5, 8);

  TProfile2D *p2_vn_pT_cent_km_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_pT_cent_km");
  p2_vn_pT_cent_km_eff->RebinY();
  TProfile *p_vn_pT_00to10_km_eff = p2_vn_pT_cent_km_eff->ProfileY("p_vn_pT_00to10_km_eff", 15, 16);
  TProfile *p_vn_pT_10to40_km_eff = p2_vn_pT_cent_km_eff->ProfileY("p_vn_pT_10to40_km_eff", 9, 14);
  TProfile *p_vn_pT_40to60_km_eff = p2_vn_pT_cent_km_eff->ProfileY("p_vn_pT_40to60_km_eff", 5, 8);

  TProfile2D *p2_vn_pT_cent_pr_eff = (TProfile2D*)efficiencyFile->Get("p2_vn_pT_cent_pr");
  TProfile *p_vn_pT_00to10_pr_eff = p2_vn_pT_cent_pr_eff->ProfileY("p_vn_pT_00to10_pr_eff", 15, 16);
  TProfile *p_vn_pT_10to40_pr_eff = p2_vn_pT_cent_pr_eff->ProfileY("p_vn_pT_10to40_pr_eff", 9, 14);
  TProfile *p_vn_pT_40to60_pr_eff = p2_vn_pT_cent_pr_eff->ProfileY("p_vn_pT_40to60_pr_eff", 5, 8);

  
  // Normal histograms
  TH1D *h_vn_pp_norm = p_vn_pp_norm->ProjectionX("h_vn_pp_norm");
  TH1D *h_vn_pm_norm = p_vn_pm_norm->ProjectionX("h_vn_pm_norm");
  TH1D *h_vn_kp_norm = p_vn_kp_norm->ProjectionX("h_vn_kp_norm");
  TH1D *h_vn_km_norm = p_vn_km_norm->ProjectionX("h_vn_km_norm");
  TH1D *h_vn_pr_norm = p_vn_pr_norm->ProjectionX("h_vn_pr_norm");

  TH1D *h_vn_pp_ext_norm = p_vn_pp_ext_norm->ProjectionX("h_vn_pp_ext_norm");
  TH1D *h_vn_pm_ext_norm = p_vn_pm_ext_norm->ProjectionX("h_vn_pm_ext_norm");
  TH1D *h_vn_kp_ext_norm = p_vn_kp_ext_norm->ProjectionX("h_vn_kp_ext_norm");
  TH1D *h_vn_km_ext_norm = p_vn_km_ext_norm->ProjectionX("h_vn_km_ext_norm");
  TH1D *h_vn_pr_ext_norm = p_vn_pr_ext_norm->ProjectionX("h_vn_pr_ext_norm");

  TH1D *h_vn_pr_for_norm = p_vn_pr_for_norm->ProjectionX("h_vn_pr_for_norm");

  TH1D *h_vn_yCM_00to10_pp_norm = p_vn_yCM_00to10_pp_norm->ProjectionX("h_vn_yCM_00to10_pp_norm");
  TH1D *h_vn_yCM_10to40_pp_norm = p_vn_yCM_10to40_pp_norm->ProjectionX("h_vn_yCM_10to40_pp_norm");
  TH1D *h_vn_yCM_40to60_pp_norm = p_vn_yCM_40to60_pp_norm->ProjectionX("h_vn_yCM_40to60_pp_norm");
  h_vn_yCM_00to10_pp_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_pp_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_pp_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_pp_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_pp_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_pp_norm->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_pm_norm = p_vn_yCM_00to10_pm_norm->ProjectionX("h_vn_yCM_00to10_pm_norm");
  TH1D *h_vn_yCM_10to40_pm_norm = p_vn_yCM_10to40_pm_norm->ProjectionX("h_vn_yCM_10to40_pm_norm");
  TH1D *h_vn_yCM_40to60_pm_norm = p_vn_yCM_40to60_pm_norm->ProjectionX("h_vn_yCM_40to60_pm_norm");
  h_vn_yCM_00to10_pm_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_pm_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_pm_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_pm_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_pm_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_pm_norm->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_kp_norm = p_vn_yCM_00to10_kp_norm->ProjectionX("h_vn_yCM_00to10_kp_norm");
  TH1D *h_vn_yCM_10to40_kp_norm = p_vn_yCM_10to40_kp_norm->ProjectionX("h_vn_yCM_10to40_kp_norm");
  TH1D *h_vn_yCM_40to60_kp_norm = p_vn_yCM_40to60_kp_norm->ProjectionX("h_vn_yCM_40to60_kp_norm");
  h_vn_yCM_00to10_kp_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_kp_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_kp_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_kp_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_kp_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_kp_norm->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_km_norm = p_vn_yCM_00to10_km_norm->ProjectionX("h_vn_yCM_00to10_km_norm");
  TH1D *h_vn_yCM_10to40_km_norm = p_vn_yCM_10to40_km_norm->ProjectionX("h_vn_yCM_10to40_km_norm");
  TH1D *h_vn_yCM_40to60_km_norm = p_vn_yCM_40to60_km_norm->ProjectionX("h_vn_yCM_40to60_km_norm");
  h_vn_yCM_00to10_km_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_km_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_km_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_km_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_km_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_km_norm->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_pr_norm = p_vn_yCM_00to10_pr_norm->ProjectionX("h_vn_yCM_00to10_pr_norm");
  TH1D *h_vn_yCM_10to40_pr_norm = p_vn_yCM_10to40_pr_norm->ProjectionX("h_vn_yCM_10to40_pr_norm");
  TH1D *h_vn_yCM_40to60_pr_norm = p_vn_yCM_40to60_pr_norm->ProjectionX("h_vn_yCM_40to60_pr_norm");
  h_vn_yCM_00to10_pr_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_pr_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_pr_norm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_pr_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_pr_norm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_pr_norm->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_pr_norm_symm = p_vn_yCM_00to10_pr_norm_symm->ProjectionX("h_vn_yCM_00to10_pr_norm_symm");
  TH1D *h_vn_yCM_10to40_pr_norm_symm = p_vn_yCM_10to40_pr_norm_symm->ProjectionX("h_vn_yCM_10to40_pr_norm_symm");
  TH1D *h_vn_yCM_40to60_pr_norm_symm = p_vn_yCM_40to60_pr_norm_symm->ProjectionX("h_vn_yCM_40to60_pr_norm_symm");
  h_vn_yCM_00to10_pr_norm_symm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_pr_norm_symm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_pr_norm_symm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_pr_norm_symm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_pr_norm_symm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_pr_norm_symm->GetYaxis()->SetTitle("v_{3}");


  TH1D *h_vn_pT_00to10_pp_norm = new TH1D("h_vn_pT_00to10_pp_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pp_norm = new TH1D("h_vn_pT_10to40_pp_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pp_norm = new TH1D("h_vn_pT_40to60_pp_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  h_vn_pT_00to10_pp_norm = p_vn_pT_00to10_pp_norm->ProjectionX();
  h_vn_pT_10to40_pp_norm = p_vn_pT_10to40_pp_norm->ProjectionX();
  h_vn_pT_40to60_pp_norm = p_vn_pT_40to60_pp_norm->ProjectionX();

  TH1D *h_vn_pT_00to10_pm_norm = new TH1D("h_vn_pT_00to10_pm_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pm_norm = new TH1D("h_vn_pT_10to40_pm_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pm_norm = new TH1D("h_vn_pT_40to60_pm_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  h_vn_pT_00to10_pm_norm = p_vn_pT_00to10_pm_norm->ProjectionX();
  h_vn_pT_10to40_pm_norm = p_vn_pT_10to40_pm_norm->ProjectionX();
  h_vn_pT_40to60_pm_norm = p_vn_pT_40to60_pm_norm->ProjectionX();

  TH1D *h_vn_pT_00to10_kp_norm = new TH1D("h_vn_pT_00to10_kp_norm", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  TH1D *h_vn_pT_10to40_kp_norm = new TH1D("h_vn_pT_10to40_kp_norm", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  TH1D *h_vn_pT_40to60_kp_norm = new TH1D("h_vn_pT_40to60_kp_norm", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  h_vn_pT_00to10_kp_norm = p_vn_pT_00to10_kp_norm->ProjectionX();
  h_vn_pT_10to40_kp_norm = p_vn_pT_10to40_kp_norm->ProjectionX();
  h_vn_pT_40to60_kp_norm = p_vn_pT_40to60_kp_norm->ProjectionX();

  TH1D *h_vn_pT_00to10_km_norm = new TH1D("h_vn_pT_00to10_km_norm", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  TH1D *h_vn_pT_10to40_km_norm = new TH1D("h_vn_pT_10to40_km_norm", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  TH1D *h_vn_pT_40to60_km_norm = new TH1D("h_vn_pT_40to60_km_norm", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  h_vn_pT_00to10_km_norm = p_vn_pT_00to10_km_norm->ProjectionX();
  h_vn_pT_10to40_km_norm = p_vn_pT_10to40_km_norm->ProjectionX();
  h_vn_pT_40to60_km_norm = p_vn_pT_40to60_km_norm->ProjectionX();

  TH1D *h_vn_pT_00to10_pr_norm = new TH1D("h_vn_pT_00to10_pr_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pr_norm = new TH1D("h_vn_pT_10to40_pr_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pr_norm = new TH1D("h_vn_pT_40to60_pr_norm", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  h_vn_pT_00to10_pr_norm = p_vn_pT_00to10_pr_norm->ProjectionX();
  h_vn_pT_10to40_pr_norm = p_vn_pT_10to40_pr_norm->ProjectionX();
  h_vn_pT_40to60_pr_norm = p_vn_pT_40to60_pr_norm->ProjectionX();


  
  // Efficiency histograms
  TH1D *h_vn_pp_eff = p_vn_pp_eff->ProjectionX("h_vn_pp_eff");
  TH1D *h_vn_pm_eff = p_vn_pm_eff->ProjectionX("h_vn_pm_eff");
  TH1D *h_vn_kp_eff = p_vn_kp_eff->ProjectionX("h_vn_kp_eff");
  TH1D *h_vn_km_eff = p_vn_km_eff->ProjectionX("h_vn_km_eff");
  TH1D *h_vn_pr_eff = p_vn_pr_eff->ProjectionX("h_vn_pr_eff");

  TH1D *h_vn_pp_ext_eff = p_vn_pp_ext_eff->ProjectionX("h_vn_pp_ext_eff");
  TH1D *h_vn_pm_ext_eff = p_vn_pm_ext_eff->ProjectionX("h_vn_pm_ext_eff");
  TH1D *h_vn_kp_ext_eff = p_vn_kp_ext_eff->ProjectionX("h_vn_kp_ext_eff");
  TH1D *h_vn_km_ext_eff = p_vn_km_ext_eff->ProjectionX("h_vn_km_ext_eff");
  TH1D *h_vn_pr_ext_eff = p_vn_pr_ext_eff->ProjectionX("h_vn_pr_ext_eff");

  TH1D *h_vn_pr_for_eff = p_vn_pr_for_eff->ProjectionX("h_vn_pr_for_eff");

  TH1D *h_vn_yCM_00to10_pp_eff = p_vn_yCM_00to10_pp_eff->ProjectionX("h_vn_yCM_00to10_pp_eff");
  TH1D *h_vn_yCM_10to40_pp_eff = p_vn_yCM_10to40_pp_eff->ProjectionX("h_vn_yCM_10to40_pp_eff");
  TH1D *h_vn_yCM_40to60_pp_eff = p_vn_yCM_40to60_pp_eff->ProjectionX("h_vn_yCM_40to60_pp_eff");
  h_vn_yCM_00to10_pp_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_pp_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_pp_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_pp_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_pp_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_pp_eff->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_pm_eff = p_vn_yCM_00to10_pm_eff->ProjectionX("h_vn_yCM_00to10_pm_eff");
  TH1D *h_vn_yCM_10to40_pm_eff = p_vn_yCM_10to40_pm_eff->ProjectionX("h_vn_yCM_10to40_pm_eff");
  TH1D *h_vn_yCM_40to60_pm_eff = p_vn_yCM_40to60_pm_eff->ProjectionX("h_vn_yCM_40to60_pm_eff");
  h_vn_yCM_00to10_pm_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_pm_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_pm_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_pm_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_pm_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_pm_eff->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_kp_eff = p_vn_yCM_00to10_kp_eff->ProjectionX("h_vn_yCM_00to10_kp_eff");
  TH1D *h_vn_yCM_10to40_kp_eff = p_vn_yCM_10to40_kp_eff->ProjectionX("h_vn_yCM_10to40_kp_eff");
  TH1D *h_vn_yCM_40to60_kp_eff = p_vn_yCM_40to60_kp_eff->ProjectionX("h_vn_yCM_40to60_kp_eff");
  h_vn_yCM_00to10_kp_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_kp_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_kp_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_kp_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_kp_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_kp_eff->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_km_eff = p_vn_yCM_00to10_km_eff->ProjectionX("h_vn_yCM_00to10_km_eff");
  TH1D *h_vn_yCM_10to40_km_eff = p_vn_yCM_10to40_km_eff->ProjectionX("h_vn_yCM_10to40_km_eff");
  TH1D *h_vn_yCM_40to60_km_eff = p_vn_yCM_40to60_km_eff->ProjectionX("h_vn_yCM_40to60_km_eff");
  h_vn_yCM_00to10_km_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_km_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_km_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_km_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_km_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_km_eff->GetYaxis()->SetTitle("v_{3}");
  
  TH1D *h_vn_yCM_00to10_pr_eff = p_vn_yCM_00to10_pr_eff->ProjectionX("h_vn_yCM_00to10_pr_eff");
  TH1D *h_vn_yCM_10to40_pr_eff = p_vn_yCM_10to40_pr_eff->ProjectionX("h_vn_yCM_10to40_pr_eff");
  TH1D *h_vn_yCM_40to60_pr_eff = p_vn_yCM_40to60_pr_eff->ProjectionX("h_vn_yCM_40to60_pr_eff");
  h_vn_yCM_00to10_pr_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_pr_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_pr_eff->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_pr_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_pr_eff->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_pr_eff->GetYaxis()->SetTitle("v_{3}");

  TH1D *h_vn_yCM_00to10_pr_eff_symm = p_vn_yCM_00to10_pr_eff_symm->ProjectionX("h_vn_yCM_00to10_pr_eff_symm");
  TH1D *h_vn_yCM_10to40_pr_eff_symm = p_vn_yCM_10to40_pr_eff_symm->ProjectionX("h_vn_yCM_10to40_pr_eff_symm");
  TH1D *h_vn_yCM_40to60_pr_eff_symm = p_vn_yCM_40to60_pr_eff_symm->ProjectionX("h_vn_yCM_40to60_pr_eff_symm");
  h_vn_yCM_00to10_pr_eff_symm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_10to40_pr_eff_symm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_40to60_pr_eff_symm->GetXaxis()->SetTitle("y-y_{mid}");
  h_vn_yCM_00to10_pr_eff_symm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_10to40_pr_eff_symm->GetYaxis()->SetTitle("v_{3}");
  h_vn_yCM_40to60_pr_eff_symm->GetYaxis()->SetTitle("v_{3}");


  TH1D *h_vn_pT_00to10_pp_eff = new TH1D("h_vn_pT_00to10_pp_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pp_eff = new TH1D("h_vn_pT_10to40_pp_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pp_eff = new TH1D("h_vn_pT_40to60_pp_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  h_vn_pT_00to10_pp_eff = p_vn_pT_00to10_pp_eff->ProjectionX();
  h_vn_pT_10to40_pp_eff = p_vn_pT_10to40_pp_eff->ProjectionX();
  h_vn_pT_40to60_pp_eff = p_vn_pT_40to60_pp_eff->ProjectionX();

  TH1D *h_vn_pT_00to10_pm_eff = new TH1D("h_vn_pT_00to10_pm_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pm_eff = new TH1D("h_vn_pT_10to40_pm_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pm_eff = new TH1D("h_vn_pT_40to60_pm_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  h_vn_pT_00to10_pm_eff = p_vn_pT_00to10_pm_eff->ProjectionX();
  h_vn_pT_10to40_pm_eff = p_vn_pT_10to40_pm_eff->ProjectionX();
  h_vn_pT_40to60_pm_eff = p_vn_pT_40to60_pm_eff->ProjectionX();

  TH1D *h_vn_pT_00to10_kp_eff = new TH1D("h_vn_pT_00to10_kp_eff", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  TH1D *h_vn_pT_10to40_kp_eff = new TH1D("h_vn_pT_10to40_kp_eff", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  TH1D *h_vn_pT_40to60_kp_eff = new TH1D("h_vn_pT_40to60_kp_eff", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  h_vn_pT_00to10_kp_eff = p_vn_pT_00to10_kp_eff->ProjectionX();
  h_vn_pT_10to40_kp_eff = p_vn_pT_10to40_kp_eff->ProjectionX();
  h_vn_pT_40to60_kp_eff = p_vn_pT_40to60_kp_eff->ProjectionX();

  TH1D *h_vn_pT_00to10_km_eff = new TH1D("h_vn_pT_00to10_km_eff", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  TH1D *h_vn_pT_10to40_km_eff = new TH1D("h_vn_pT_10to40_km_eff", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  TH1D *h_vn_pT_40to60_km_eff = new TH1D("h_vn_pT_40to60_km_eff", ";p_{T} (GeV);v_{3}", 5, 0, 2);
  h_vn_pT_00to10_km_eff = p_vn_pT_00to10_km_eff->ProjectionX();
  h_vn_pT_10to40_km_eff = p_vn_pT_10to40_km_eff->ProjectionX();
  h_vn_pT_40to60_km_eff = p_vn_pT_40to60_km_eff->ProjectionX();

  TH1D *h_vn_pT_00to10_pr_eff = new TH1D("h_vn_pT_00to10_pr_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pr_eff = new TH1D("h_vn_pT_10to40_pr_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pr_eff = new TH1D("h_vn_pT_40to60_pr_eff", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  h_vn_pT_00to10_pr_eff = p_vn_pT_00to10_pr_eff->ProjectionX();
  h_vn_pT_10to40_pr_eff = p_vn_pT_10to40_pr_eff->ProjectionX();
  h_vn_pT_40to60_pr_eff = p_vn_pT_40to60_pr_eff->ProjectionX();

  
  // Flip and truncate empty bins
  TH1D *h_vn_pp_norm_flip = new TH1D("h_vn_pp_norm_flip",h_vn_pp_norm->GetTitle(),16,0,16);
  h_vn_pp_norm_flip->GetXaxis()->SetTitle((TString)h_vn_pp_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pp_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pp_norm->GetYaxis()->GetTitle());

  TH1D *h_vn_pm_norm_flip = new TH1D("h_vn_pm_norm_flip",h_vn_pm_norm->GetTitle(),16,0,16);
  h_vn_pm_norm_flip->GetXaxis()->SetTitle((TString)h_vn_pm_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pm_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pm_norm->GetYaxis()->GetTitle());

  TH1D *h_vn_kp_norm_flip = new TH1D("h_vn_kp_norm_flip",h_vn_kp_norm->GetTitle(),8,0,16);
  h_vn_kp_norm_flip->GetXaxis()->SetTitle((TString)h_vn_kp_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_kp_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_kp_norm->GetYaxis()->GetTitle());

  TH1D *h_vn_km_norm_flip = new TH1D("h_vn_km_norm_flip",h_vn_km_norm->GetTitle(),8,0,16);
  h_vn_km_norm_flip->GetXaxis()->SetTitle((TString)h_vn_km_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_km_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_km_norm->GetYaxis()->GetTitle());

  TH1D *h_vn_pr_norm_flip = new TH1D("h_vn_pr_norm_flip",h_vn_pr_norm->GetTitle(),16,0,16);
  h_vn_pr_norm_flip->GetXaxis()->SetTitle((TString)h_vn_pr_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pr_norm->GetYaxis()->GetTitle());

  
  TH1D *h_vn_pp_eff_flip = new TH1D("h_vn_pp_eff_flip",h_vn_pp_eff->GetTitle(),16,0,16);
  h_vn_pp_eff_flip->GetXaxis()->SetTitle((TString)h_vn_pp_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_pp_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pp_eff->GetYaxis()->GetTitle());

  TH1D *h_vn_pm_eff_flip = new TH1D("h_vn_pm_eff_flip",h_vn_pm_eff->GetTitle(),16,0,16);
  h_vn_pm_eff_flip->GetXaxis()->SetTitle((TString)h_vn_pm_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_pm_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pm_eff->GetYaxis()->GetTitle());

  TH1D *h_vn_kp_eff_flip = new TH1D("h_vn_kp_eff_flip",h_vn_kp_eff->GetTitle(),8,0,16);
  h_vn_kp_eff_flip->GetXaxis()->SetTitle((TString)h_vn_kp_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_kp_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_kp_eff->GetYaxis()->GetTitle());

  TH1D *h_vn_km_eff_flip = new TH1D("h_vn_km_eff_flip",h_vn_km_eff->GetTitle(),8,0,16);
  h_vn_km_eff_flip->GetXaxis()->SetTitle((TString)h_vn_km_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_km_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_km_eff->GetYaxis()->GetTitle());

  TH1D *h_vn_pr_eff_flip = new TH1D("h_vn_pr_eff_flip",h_vn_pr_eff->GetTitle(),16,0,16);
  h_vn_pr_eff_flip->GetXaxis()->SetTitle((TString)h_vn_pr_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pr_eff->GetYaxis()->GetTitle());


  TH1D *h_vn_pp_ext_norm_flip = new TH1D("h_vn_pp_ext_norm_flip",h_vn_pp_ext_norm->GetTitle(),16,0,16);
  h_vn_pp_ext_norm_flip->GetXaxis()->SetTitle((TString)h_vn_pp_ext_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pp_ext_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pp_ext_norm->GetYaxis()->GetTitle());

  TH1D *h_vn_pm_ext_norm_flip = new TH1D("h_vn_pm_ext_norm_flip",h_vn_pm_ext_norm->GetTitle(),16,0,16);
  h_vn_pm_ext_norm_flip->GetXaxis()->SetTitle((TString)h_vn_pm_ext_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pm_ext_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pm_ext_norm->GetYaxis()->GetTitle());

  TH1D *h_vn_kp_ext_norm_flip = new TH1D("h_vn_kp_ext_norm_flip",h_vn_kp_ext_norm->GetTitle(),8,0,16);
  h_vn_kp_ext_norm_flip->GetXaxis()->SetTitle((TString)h_vn_kp_ext_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_kp_ext_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_kp_ext_norm->GetYaxis()->GetTitle());

  TH1D *h_vn_km_ext_norm_flip = new TH1D("h_vn_km_ext_norm_flip",h_vn_km_ext_norm->GetTitle(),8,0,16);
  h_vn_km_ext_norm_flip->GetXaxis()->SetTitle((TString)h_vn_km_ext_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_km_ext_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_km_ext_norm->GetYaxis()->GetTitle());

  TH1D *h_vn_pr_ext_norm_flip = new TH1D("h_vn_pr_ext_norm_flip",h_vn_pr_ext_norm->GetTitle(),16,0,16);
  h_vn_pr_ext_norm_flip->GetXaxis()->SetTitle((TString)h_vn_pr_ext_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_ext_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pr_ext_norm->GetYaxis()->GetTitle());


  TH1D *h_vn_pr_for_norm_flip = new TH1D("h_vn_pr_for_norm_flip",h_vn_pr_for_norm->GetTitle(),16,0,16);
  h_vn_pr_for_norm_flip->GetXaxis()->SetTitle((TString)h_vn_pr_for_norm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_for_norm_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pr_for_norm->GetYaxis()->GetTitle());


  
  TH1D *h_vn_pp_ext_eff_flip = new TH1D("h_vn_pp_ext_eff_flip",h_vn_pp_ext_eff->GetTitle(),16,0,16);
  h_vn_pp_ext_eff_flip->GetXaxis()->SetTitle((TString)h_vn_pp_ext_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_pp_ext_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pp_ext_eff->GetYaxis()->GetTitle());

  TH1D *h_vn_pm_ext_eff_flip = new TH1D("h_vn_pm_ext_eff_flip",h_vn_pm_ext_eff->GetTitle(),16,0,16);
  h_vn_pm_ext_eff_flip->GetXaxis()->SetTitle((TString)h_vn_pm_ext_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_pm_ext_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pm_ext_eff->GetYaxis()->GetTitle());

  TH1D *h_vn_kp_ext_eff_flip = new TH1D("h_vn_kp_ext_eff_flip",h_vn_kp_ext_eff->GetTitle(),8,0,16);
  h_vn_kp_ext_eff_flip->GetXaxis()->SetTitle((TString)h_vn_kp_ext_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_kp_ext_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_kp_ext_eff->GetYaxis()->GetTitle());

  TH1D *h_vn_km_ext_eff_flip = new TH1D("h_vn_km_ext_eff_flip",h_vn_km_ext_eff->GetTitle(),8,0,16);
  h_vn_km_ext_eff_flip->GetXaxis()->SetTitle((TString)h_vn_km_ext_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_km_ext_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_km_ext_eff->GetYaxis()->GetTitle());

  TH1D *h_vn_pr_ext_eff_flip = new TH1D("h_vn_pr_ext_eff_flip",h_vn_pr_ext_eff->GetTitle(),16,0,16);
  h_vn_pr_ext_eff_flip->GetXaxis()->SetTitle((TString)h_vn_pr_ext_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_ext_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pr_ext_eff->GetYaxis()->GetTitle());


  TH1D *h_vn_pr_for_eff_flip = new TH1D("h_vn_pr_for_eff_flip",h_vn_pr_for_eff->GetTitle(),16,0,16);
  h_vn_pr_for_eff_flip->GetXaxis()->SetTitle((TString)h_vn_pr_for_eff->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_for_eff_flip->GetYaxis()->SetTitle("v_{3}");//h_vn_pr_for_eff->GetYaxis()->GetTitle());



  // Flip the bin contents into the new histograms
  int j = 1;
  for (int i = 8; i >= 1; i--)     // KAONS
    {
      h_vn_kp_norm_flip->SetBinContent(j, h_vn_kp_norm->GetBinContent(i));
      h_vn_kp_norm_flip->SetBinError(j, h_vn_kp_norm->GetBinError(i));

      h_vn_kp_eff_flip->SetBinContent(j, h_vn_kp_eff->GetBinContent(i));
      h_vn_kp_eff_flip->SetBinError(j, h_vn_kp_eff->GetBinError(i));

      h_vn_km_norm_flip->SetBinContent(j, h_vn_km_norm->GetBinContent(i));
      h_vn_km_norm_flip->SetBinError(j, h_vn_km_norm->GetBinError(i));

      h_vn_km_eff_flip->SetBinContent(j, h_vn_km_eff->GetBinContent(i));
      h_vn_km_eff_flip->SetBinError(j, h_vn_km_eff->GetBinError(i));

      h_vn_kp_ext_norm_flip->SetBinContent(j, h_vn_kp_ext_norm->GetBinContent(i));
      h_vn_kp_ext_norm_flip->SetBinError(j, h_vn_kp_ext_norm->GetBinError(i));

      h_vn_kp_ext_eff_flip->SetBinContent(j, h_vn_kp_ext_eff->GetBinContent(i));
      h_vn_kp_ext_eff_flip->SetBinError(j, h_vn_kp_ext_eff->GetBinError(i));

      h_vn_km_ext_norm_flip->SetBinContent(j, h_vn_km_ext_norm->GetBinContent(i));
      h_vn_km_ext_norm_flip->SetBinError(j, h_vn_km_ext_norm->GetBinError(i));

      h_vn_km_ext_eff_flip->SetBinContent(j, h_vn_km_ext_eff->GetBinContent(i));
      h_vn_km_ext_eff_flip->SetBinError(j, h_vn_km_ext_eff->GetBinError(i));

      j++;
    }

  j = 1;
  for (int i = 16; i >= 1; i--)
    {
      /*
      h_vn_EpdE_flip->SetBinContent(j, h_vn_EpdE->GetBinContent(i));
      h_vn_EpdE_flip->SetBinError(j, h_vn_EpdE->GetBinError(i));

      h_vn_EpdF_flip->SetBinContent(j, h_vn_EpdF->GetBinContent(i));
      h_vn_EpdF_flip->SetBinError(j, h_vn_EpdF->GetBinError(i));

      h_vn_TpcB_flip->SetBinContent(j, h_vn_TpcB->GetBinContent(i));
      h_vn_TpcB_flip->SetBinError(j, h_vn_TpcB->GetBinError(i));

      h_vn_Tpc_flip->SetBinContent(j, h_vn_Tpc->GetBinContent(i));
      h_vn_Tpc_flip->SetBinError(j, h_vn_Tpc->GetBinError(i));
      */
      h_vn_pp_norm_flip->SetBinContent(j, h_vn_pp_norm->GetBinContent(i));
      h_vn_pp_norm_flip->SetBinError(j, h_vn_pp_norm->GetBinError(i));

      h_vn_pp_eff_flip->SetBinContent(j, h_vn_pp_eff->GetBinContent(i));
      h_vn_pp_eff_flip->SetBinError(j, h_vn_pp_eff->GetBinError(i));

      h_vn_pm_norm_flip->SetBinContent(j, h_vn_pm_norm->GetBinContent(i));
      h_vn_pm_norm_flip->SetBinError(j, h_vn_pm_norm->GetBinError(i));

      h_vn_pm_eff_flip->SetBinContent(j, h_vn_pm_eff->GetBinContent(i));
      h_vn_pm_eff_flip->SetBinError(j, h_vn_pm_eff->GetBinError(i));

      h_vn_pr_norm_flip->SetBinContent(j, h_vn_pr_norm->GetBinContent(i));
      h_vn_pr_norm_flip->SetBinError(j, h_vn_pr_norm->GetBinError(i));

      h_vn_pr_eff_flip->SetBinContent(j, h_vn_pr_eff->GetBinContent(i));
      h_vn_pr_eff_flip->SetBinError(j, h_vn_pr_eff->GetBinError(i));


      h_vn_pp_ext_norm_flip->SetBinContent(j, h_vn_pp_ext_norm->GetBinContent(i));
      h_vn_pp_ext_norm_flip->SetBinError(j, h_vn_pp_ext_norm->GetBinError(i));

      h_vn_pp_ext_eff_flip->SetBinContent(j, h_vn_pp_ext_eff->GetBinContent(i));
      h_vn_pp_ext_eff_flip->SetBinError(j, h_vn_pp_ext_eff->GetBinError(i));

      h_vn_pm_ext_norm_flip->SetBinContent(j, h_vn_pm_ext_norm->GetBinContent(i));
      h_vn_pm_ext_norm_flip->SetBinError(j, h_vn_pm_ext_norm->GetBinError(i));

      h_vn_pm_ext_eff_flip->SetBinContent(j, h_vn_pm_ext_eff->GetBinContent(i));
      h_vn_pm_ext_eff_flip->SetBinError(j, h_vn_pm_ext_eff->GetBinError(i));

      h_vn_pr_ext_norm_flip->SetBinContent(j, h_vn_pr_ext_norm->GetBinContent(i));
      h_vn_pr_ext_norm_flip->SetBinError(j, h_vn_pr_ext_norm->GetBinError(i));

      h_vn_pr_ext_eff_flip->SetBinContent(j, h_vn_pr_ext_eff->GetBinContent(i));
      h_vn_pr_ext_eff_flip->SetBinError(j, h_vn_pr_ext_eff->GetBinError(i));


      h_vn_pr_for_norm_flip->SetBinContent(j, h_vn_pr_for_norm->GetBinContent(i));
      h_vn_pr_for_norm_flip->SetBinError(j, h_vn_pr_for_norm->GetBinError(i));

      h_vn_pr_for_eff_flip->SetBinContent(j, h_vn_pr_for_eff->GetBinContent(i));
      h_vn_pr_for_eff_flip->SetBinError(j, h_vn_pr_for_eff->GetBinError(i));

      
      j++;
    }


  // Use these to change the x axis of centrality plots without fighting with TAxis.
  /*
  TH1D *vn_EpdE = new TH1D("vn_EpdE", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_EpdF = new TH1D("vn_EpdF", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_TpcB = new TH1D("vn_TpcB", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_Tpc = new TH1D("vn_Tpc", ";Centrality (%);v_{3}", 12, 0, 60);
  */
  TH1D *vn_pp_norm = new TH1D("vn_pp_norm", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_pm_norm = new TH1D("vn_pm_norm", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_kp_norm = new TH1D("vn_kp_norm", ";Centrality (%);v_{3}", 6, 0, 60);
  TH1D *vn_km_norm = new TH1D("vn_km_norm", ";Centrality (%);v_{3}", 6, 0, 60);
  TH1D *vn_pr_norm = new TH1D("vn_pr_norm", ";Centrality (%);v_{3}", 12, 0, 60);

  TH1D *vn_pp_eff = new TH1D("vn_pp_eff", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_pm_eff = new TH1D("vn_pm_eff", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_kp_eff = new TH1D("vn_kp_eff", ";Centrality (%);v_{3}", 6, 0, 60);
  TH1D *vn_km_eff = new TH1D("vn_km_eff", ";Centrality (%);v_{3}", 6, 0, 60);
  TH1D *vn_pr_eff = new TH1D("vn_pr_eff", ";Centrality (%);v_{3}", 12, 0, 60);

  TH1D *vn_pp_ext_norm = new TH1D("vn_pp_ext_norm", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_pm_ext_norm = new TH1D("vn_pm_ext_norm", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_kp_ext_norm = new TH1D("vn_kp_ext_norm", ";Centrality (%);v_{3}", 6, 0, 60);
  TH1D *vn_km_ext_norm = new TH1D("vn_km_ext_norm", ";Centrality (%);v_{3}", 6, 0, 60);
  TH1D *vn_pr_ext_norm = new TH1D("vn_pr_ext_norm", ";Centrality (%);v_{3}", 12, 0, 60);

  TH1D *vn_pp_ext_eff = new TH1D("vn_pp_ext_eff", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_pm_ext_eff = new TH1D("vn_pm_ext_eff", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_kp_ext_eff = new TH1D("vn_kp_ext_eff", ";Centrality (%);v_{3}", 6, 0, 60);
  TH1D *vn_km_ext_eff = new TH1D("vn_km_ext_eff", ";Centrality (%);v_{3}", 6, 0, 60);
  TH1D *vn_pr_ext_eff = new TH1D("vn_pr_ext_eff", ";Centrality (%);v_{3}", 12, 0, 60);

  TH1D *vn_pr_for_norm = new TH1D("vn_pr_for_norm", ";Centrality (%);v_{3}", 12, 0, 60);
  TH1D *vn_pr_for_eff = new TH1D("vn_pr_for_eff", ";Centrality (%);v_{3}", 12, 0, 60);

  
  // Move to final histos to truncate empty bins
  for (int i = 1; i <= 6; i++)     // KAONS
    {
      vn_kp_norm->SetBinContent(i, h_vn_kp_norm_flip->GetBinContent(i));
      vn_kp_norm->SetBinError(i, h_vn_kp_norm_flip->GetBinError(i));

      vn_km_norm->SetBinContent(i, h_vn_km_norm_flip->GetBinContent(i));
      vn_km_norm->SetBinError(i, h_vn_km_norm_flip->GetBinError(i));

      vn_kp_eff->SetBinContent(i, h_vn_kp_eff_flip->GetBinContent(i));
      vn_kp_eff->SetBinError(i, h_vn_kp_eff_flip->GetBinError(i));

      vn_km_eff->SetBinContent(i, h_vn_km_eff_flip->GetBinContent(i));
      vn_km_eff->SetBinError(i, h_vn_km_eff_flip->GetBinError(i));


      vn_kp_ext_norm->SetBinContent(i, h_vn_kp_ext_norm_flip->GetBinContent(i));
      vn_kp_ext_norm->SetBinError(i, h_vn_kp_ext_norm_flip->GetBinError(i));

      vn_km_ext_norm->SetBinContent(i, h_vn_km_ext_norm_flip->GetBinContent(i));
      vn_km_ext_norm->SetBinError(i, h_vn_km_ext_norm_flip->GetBinError(i));

      vn_kp_ext_eff->SetBinContent(i, h_vn_kp_ext_eff_flip->GetBinContent(i));
      vn_kp_ext_eff->SetBinError(i, h_vn_kp_ext_eff_flip->GetBinError(i));

      vn_km_ext_eff->SetBinContent(i, h_vn_km_ext_eff_flip->GetBinContent(i));
      vn_km_ext_eff->SetBinError(i, h_vn_km_ext_eff_flip->GetBinError(i));
    }

  for (int i = 1; i <= 12; i++)
    {
      /*
      vn_EpdE->SetBinContent(i, h_vn_EpdE_flip->GetBinContent(i));
      vn_EpdE->SetBinError(i, h_vn_EpdE_flip->GetBinError(i));

      vn_EpdF->SetBinContent(i, h_vn_EpdF_flip->GetBinContent(i));
      vn_EpdF->SetBinError(i, h_vn_EpdF_flip->GetBinError(i));

      vn_TpcB->SetBinContent(i, h_vn_TpcB_flip->GetBinContent(i));
      vn_TpcB->SetBinError(i, h_vn_TpcB_flip->GetBinError(i));

      vn_Tpc->SetBinContent(i, h_vn_Tpc_flip->GetBinContent(i));
      vn_Tpc->SetBinError(i, h_vn_Tpc_flip->GetBinError(i));
      */
      vn_pp_norm->SetBinContent(i, h_vn_pp_norm_flip->GetBinContent(i));
      vn_pp_norm->SetBinError(i, h_vn_pp_norm_flip->GetBinError(i));

      vn_pm_norm->SetBinContent(i, h_vn_pm_norm_flip->GetBinContent(i));
      vn_pm_norm->SetBinError(i, h_vn_pm_norm_flip->GetBinError(i));
      
      vn_pr_norm->SetBinContent(i, h_vn_pr_norm_flip->GetBinContent(i));
      vn_pr_norm->SetBinError(i, h_vn_pr_norm_flip->GetBinError(i));

      vn_pp_eff->SetBinContent(i, h_vn_pp_eff_flip->GetBinContent(i));
      vn_pp_eff->SetBinError(i, h_vn_pp_eff_flip->GetBinError(i));
      
      vn_pm_eff->SetBinContent(i, h_vn_pm_eff_flip->GetBinContent(i));
      vn_pm_eff->SetBinError(i, h_vn_pm_eff_flip->GetBinError(i));

      vn_pr_eff->SetBinContent(i, h_vn_pr_eff_flip->GetBinContent(i));
      vn_pr_eff->SetBinError(i, h_vn_pr_eff_flip->GetBinError(i));


      vn_pp_ext_norm->SetBinContent(i, h_vn_pp_ext_norm_flip->GetBinContent(i));
      vn_pp_ext_norm->SetBinError(i, h_vn_pp_ext_norm_flip->GetBinError(i));

      vn_pm_ext_norm->SetBinContent(i, h_vn_pm_ext_norm_flip->GetBinContent(i));
      vn_pm_ext_norm->SetBinError(i, h_vn_pm_ext_norm_flip->GetBinError(i));

      vn_pr_ext_norm->SetBinContent(i, h_vn_pr_ext_norm_flip->GetBinContent(i));
      vn_pr_ext_norm->SetBinError(i, h_vn_pr_ext_norm_flip->GetBinError(i));

      vn_pr_for_norm->SetBinContent(i, h_vn_pr_for_norm_flip->GetBinContent(i));
      vn_pr_for_norm->SetBinError(i, h_vn_pr_for_norm_flip->GetBinError(i));

      vn_pp_ext_eff->SetBinContent(i, h_vn_pp_ext_eff_flip->GetBinContent(i));
      vn_pp_ext_eff->SetBinError(i, h_vn_pp_ext_eff_flip->GetBinError(i));

      vn_pm_ext_eff->SetBinContent(i, h_vn_pm_ext_eff_flip->GetBinContent(i));
      vn_pm_ext_eff->SetBinError(i, h_vn_pm_ext_eff_flip->GetBinError(i));

      vn_pr_ext_eff->SetBinContent(i, h_vn_pr_ext_eff_flip->GetBinContent(i));
      vn_pr_ext_eff->SetBinError(i, h_vn_pr_ext_eff_flip->GetBinError(i));

      vn_pr_for_eff->SetBinContent(i, h_vn_pr_for_eff_flip->GetBinContent(i));
      vn_pr_for_eff->SetBinError(i, h_vn_pr_for_eff_flip->GetBinError(i));

    }

  vn_pp_eff->SetLineColor(kRed);
  vn_pp_eff->SetMarkerColor(kRed);
  vn_pm_eff->SetLineColor(kRed);
  vn_pm_eff->SetMarkerColor(kRed);
  vn_kp_eff->SetLineColor(kRed);
  vn_kp_eff->SetMarkerColor(kRed);
  vn_km_eff->SetLineColor(kRed);
  vn_km_eff->SetMarkerColor(kRed);
  vn_pr_eff->SetLineColor(kRed);
  vn_pr_eff->SetMarkerColor(kRed);

  vn_pp_ext_eff->SetLineColor(kRed);
  vn_pp_ext_eff->SetMarkerColor(kRed);
  vn_pm_ext_eff->SetLineColor(kRed);
  vn_pm_ext_eff->SetMarkerColor(kRed);
  vn_kp_ext_eff->SetLineColor(kRed);
  vn_kp_ext_eff->SetMarkerColor(kRed);
  vn_km_ext_eff->SetLineColor(kRed);
  vn_km_ext_eff->SetMarkerColor(kRed);
  vn_pr_ext_eff->SetLineColor(kRed);
  vn_pr_ext_eff->SetMarkerColor(kRed);

  vn_pr_for_eff->SetLineColor(kRed);
  vn_pr_for_eff->SetMarkerColor(kRed);

  h_vn_yCM_00to10_pp_eff->SetLineColor(kRed);
  h_vn_yCM_00to10_pp_eff->SetMarkerColor(kRed);
  h_vn_yCM_10to40_pp_eff->SetLineColor(kRed);
  h_vn_yCM_10to40_pp_eff->SetMarkerColor(kRed);
  h_vn_yCM_40to60_pp_eff->SetLineColor(kRed);
  h_vn_yCM_40to60_pp_eff->SetMarkerColor(kRed);

  h_vn_yCM_00to10_pm_eff->SetLineColor(kRed);
  h_vn_yCM_00to10_pm_eff->SetMarkerColor(kRed);
  h_vn_yCM_10to40_pm_eff->SetLineColor(kRed);
  h_vn_yCM_10to40_pm_eff->SetMarkerColor(kRed);
  h_vn_yCM_40to60_pm_eff->SetLineColor(kRed);
  h_vn_yCM_40to60_pm_eff->SetMarkerColor(kRed);

  h_vn_yCM_00to10_kp_eff->SetLineColor(kRed);
  h_vn_yCM_00to10_kp_eff->SetMarkerColor(kRed);
  h_vn_yCM_10to40_kp_eff->SetLineColor(kRed);
  h_vn_yCM_10to40_kp_eff->SetMarkerColor(kRed);
  h_vn_yCM_40to60_kp_eff->SetLineColor(kRed);
  h_vn_yCM_40to60_kp_eff->SetMarkerColor(kRed);

  h_vn_yCM_10to40_km_eff->SetLineColor(kRed);
  h_vn_yCM_10to40_km_eff->SetMarkerColor(kRed);

  h_vn_yCM_00to10_pr_eff->SetLineColor(kRed);
  h_vn_yCM_00to10_pr_eff->SetMarkerColor(kRed);
  h_vn_yCM_10to40_pr_eff->SetLineColor(kRed);
  h_vn_yCM_10to40_pr_eff->SetMarkerColor(kRed);
  h_vn_yCM_40to60_pr_eff->SetLineColor(kRed);
  h_vn_yCM_40to60_pr_eff->SetMarkerColor(kRed);

  h_vn_yCM_00to10_pr_eff_symm->SetLineColor(kRed);
  h_vn_yCM_00to10_pr_eff_symm->SetMarkerColor(kRed);
  h_vn_yCM_10to40_pr_eff_symm->SetLineColor(kRed);
  h_vn_yCM_10to40_pr_eff_symm->SetMarkerColor(kRed);
  h_vn_yCM_40to60_pr_eff_symm->SetLineColor(kRed);
  h_vn_yCM_40to60_pr_eff_symm->SetMarkerColor(kRed);


  h_vn_pT_00to10_pp_eff->SetLineColor(kRed);
  h_vn_pT_00to10_pp_eff->SetMarkerColor(kRed);
  h_vn_pT_10to40_pp_eff->SetLineColor(kRed);
  h_vn_pT_10to40_pp_eff->SetMarkerColor(kRed);
  h_vn_pT_40to60_pp_eff->SetLineColor(kRed);
  h_vn_pT_40to60_pp_eff->SetMarkerColor(kRed);

  h_vn_pT_00to10_pm_eff->SetLineColor(kRed);
  h_vn_pT_00to10_pm_eff->SetMarkerColor(kRed);
  h_vn_pT_10to40_pm_eff->SetLineColor(kRed);
  h_vn_pT_10to40_pm_eff->SetMarkerColor(kRed);
  h_vn_pT_40to60_pm_eff->SetLineColor(kRed);
  h_vn_pT_40to60_pm_eff->SetMarkerColor(kRed);

  h_vn_pT_00to10_kp_eff->SetLineColor(kRed);
  h_vn_pT_00to10_kp_eff->SetMarkerColor(kRed);
  h_vn_pT_10to40_kp_eff->SetLineColor(kRed);
  h_vn_pT_10to40_kp_eff->SetMarkerColor(kRed);
  h_vn_pT_40to60_kp_eff->SetLineColor(kRed);
  h_vn_pT_40to60_kp_eff->SetMarkerColor(kRed);

  h_vn_pT_10to40_km_eff->SetLineColor(kRed);
  h_vn_pT_10to40_km_eff->SetMarkerColor(kRed);

  h_vn_pT_00to10_pr_eff->SetLineColor(kRed);
  h_vn_pT_00to10_pr_eff->SetMarkerColor(kRed);
  h_vn_pT_10to40_pr_eff->SetLineColor(kRed);
  h_vn_pT_10to40_pr_eff->SetMarkerColor(kRed);
  h_vn_pT_40to60_pr_eff->SetLineColor(kRed);
  h_vn_pT_40to60_pr_eff->SetMarkerColor(kRed);

  
  TRatioPlot *rp_pp = new TRatioPlot(vn_pp_norm, vn_pp_eff, "diffsig");
  rp_pp->Draw();
  rp_pp->GetUpperRefYaxis()->SetRangeUser(-0.015,0.02);
  rp_pp->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pp->SetLeftMargin(0.15);
  rp_pp->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pp.png");
  canvas->Clear();
  delete rp_pp;

  TRatioPlot *rp_pm = new TRatioPlot(vn_pm_norm, vn_pm_eff, "diffsig");
  rp_pm->Draw();
  rp_pm->GetUpperRefYaxis()->SetRangeUser(-0.015,0.02);
  rp_pm->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pm->SetLeftMargin(0.13);
  rp_pm->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pm.png");
  canvas->Clear();
  delete rp_pm;

  TRatioPlot *rp_kp = new TRatioPlot(vn_kp_norm, vn_kp_eff, "diffsig");
  rp_kp->Draw();
  rp_kp->GetUpperRefYaxis()->SetRangeUser(-0.15,0.2);
  rp_kp->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_kp->SetLeftMargin(0.13);
  rp_kp->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_kp.png");
  canvas->Clear();
  delete rp_kp;

  TRatioPlot *rp_km = new TRatioPlot(vn_km_norm, vn_km_eff, "diffsig");
  rp_km->Draw();
  rp_km->GetUpperRefYaxis()->SetRangeUser(-0.15,0.2);
  rp_km->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_km->SetLeftMargin(0.13);
  rp_km->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_km.png");
  canvas->Clear();
  delete rp_km;

  TRatioPlot *rp_pr = new TRatioPlot(vn_pr_norm, vn_pr_eff, "diffsig");
  rp_pr->Draw();
  rp_pr->GetUpperRefYaxis()->SetRangeUser(-0.03,0.01);
  rp_pr->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr->SetLeftMargin(0.15);
  rp_pr->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr.png");
  canvas->Clear();
  delete rp_pr;


  TRatioPlot *rp_pp_ext = new TRatioPlot(vn_pp_ext_norm, vn_pp_ext_eff, "diffsig");
  rp_pp_ext->Draw();
  rp_pp_ext->GetUpperRefYaxis()->SetRangeUser(-0.015,0.04);
  rp_pp_ext->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pp_ext->SetLeftMargin(0.13);
  rp_pp_ext->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pp_ext.png");
  canvas->Clear();
  delete rp_pp_ext;

  TRatioPlot *rp_pm_ext = new TRatioPlot(vn_pm_ext_norm, vn_pm_ext_eff, "diffsig");
  rp_pm_ext->Draw();
  rp_pm_ext->GetUpperRefYaxis()->SetRangeUser(-0.015,0.02);
  rp_pm_ext->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pm_ext->SetLeftMargin(0.13);
  rp_pm_ext->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pm_ext.png");
  canvas->Clear();
  delete rp_pm_ext;

  TRatioPlot *rp_kp_ext = new TRatioPlot(vn_kp_ext_norm, vn_kp_ext_eff, "diffsig");
  rp_kp_ext->Draw();
  rp_kp_ext->GetUpperRefYaxis()->SetRangeUser(-0.1,0.1);
  rp_kp_ext->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_kp_ext->SetLeftMargin(0.13);
  rp_kp_ext->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_kp_ext.png");
  canvas->Clear();
  delete rp_kp_ext;

  TRatioPlot *rp_km_ext = new TRatioPlot(vn_km_ext_norm, vn_km_ext_eff, "diffsig");
  rp_km_ext->Draw();
  rp_km_ext->GetUpperRefYaxis()->SetRangeUser(-0.3,0.3);
  rp_km_ext->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_km_ext->SetLeftMargin(0.13);
  rp_km_ext->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_km_ext.png");
  canvas->Clear();
  delete rp_km_ext;

  TRatioPlot *rp_pr_ext = new TRatioPlot(vn_pr_ext_norm, vn_pr_ext_eff, "diffsig");
  rp_pr_ext->Draw();
  rp_pr_ext->GetUpperRefYaxis()->SetRangeUser(-0.08,0.01);
  rp_pr_ext->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_ext->SetLeftMargin(0.13);
  rp_pr_ext->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_ext.png");
  canvas->Clear();
  delete rp_pr_ext;

  TRatioPlot *rp_pr_for = new TRatioPlot(vn_pr_for_norm, vn_pr_for_eff, "diffsig");
  rp_pr_for->Draw();
  rp_pr_for->GetUpperRefYaxis()->SetRangeUser(-0.02,0.05);
  rp_pr_for->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_for->SetLeftMargin(0.13);
  rp_pr_for->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_for.png");
  canvas->Clear();
  delete rp_pr_for;


  TRatioPlot *rp_pp_yCM_00to10 = new TRatioPlot(h_vn_yCM_00to10_pp_norm, h_vn_yCM_00to10_pp_eff, "diffsig");
  rp_pp_yCM_00to10->Draw();
  rp_pp_yCM_00to10->GetUpperRefYaxis()->SetRangeUser(-0.02,0.04);
  rp_pp_yCM_00to10->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pp_yCM_00to10->SetLeftMargin(0.13);
  rp_pp_yCM_00to10->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pp_yCM_00to10.png");
  canvas->Clear();
  delete rp_pp_yCM_00to10;

  TRatioPlot *rp_pp_yCM_10to40 = new TRatioPlot(h_vn_yCM_10to40_pp_norm, h_vn_yCM_10to40_pp_eff, "diffsig");
  rp_pp_yCM_10to40->Draw();
  rp_pp_yCM_10to40->GetUpperRefYaxis()->SetRangeUser(-0.02,0.04);
  rp_pp_yCM_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pp_yCM_10to40->SetLeftMargin(0.13);
  rp_pp_yCM_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pp_yCM_10to40.png");
  canvas->Clear();
  delete rp_pp_yCM_10to40;

  TRatioPlot *rp_pp_yCM_40to60 = new TRatioPlot(h_vn_yCM_40to60_pp_norm, h_vn_yCM_40to60_pp_eff, "diffsig");
  rp_pp_yCM_40to60->Draw();
  rp_pp_yCM_40to60->GetUpperRefYaxis()->SetRangeUser(-0.02,0.04);
  rp_pp_yCM_40to60->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pp_yCM_40to60->SetLeftMargin(0.13);
  rp_pp_yCM_40to60->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pp_yCM_40to60.png");
  canvas->Clear();
  delete rp_pp_yCM_40to60;

  TRatioPlot *rp_pm_yCM_00to10 = new TRatioPlot(h_vn_yCM_00to10_pm_norm, h_vn_yCM_00to10_pm_eff, "diffsig");
  rp_pm_yCM_00to10->Draw();
  rp_pm_yCM_00to10->GetUpperRefYaxis()->SetRangeUser(-0.02,0.04);
  rp_pm_yCM_00to10->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pm_yCM_00to10->SetLeftMargin(0.13);
  rp_pm_yCM_00to10->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pm_yCM_00to10.png");
  canvas->Clear();
  delete rp_pm_yCM_00to10;

  TRatioPlot *rp_pm_yCM_10to40 = new TRatioPlot(h_vn_yCM_10to40_pm_norm, h_vn_yCM_10to40_pm_eff, "diffsig");
  rp_pm_yCM_10to40->Draw();
  rp_pm_yCM_10to40->GetUpperRefYaxis()->SetRangeUser(-0.02,0.04);
  rp_pm_yCM_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pm_yCM_10to40->SetLeftMargin(0.13);
  rp_pm_yCM_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pm_yCM_10to40.png");
  canvas->Clear();
  delete rp_pm_yCM_10to40;

  TRatioPlot *rp_pm_yCM_40to60 = new TRatioPlot(h_vn_yCM_40to60_pm_norm, h_vn_yCM_40to60_pm_eff, "diffsig");
  rp_pm_yCM_40to60->Draw();
  rp_pm_yCM_40to60->GetUpperRefYaxis()->SetRangeUser(-0.02,0.04);
  rp_pm_yCM_40to60->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pm_yCM_40to60->SetLeftMargin(0.13);
  rp_pm_yCM_40to60->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pm_yCM_40to60.png");
  canvas->Clear();
  delete rp_pm_yCM_40to60;

  TRatioPlot *rp_kp_yCM_00to10 = new TRatioPlot(h_vn_yCM_00to10_kp_norm, h_vn_yCM_00to10_kp_eff, "diffsig");
  rp_kp_yCM_00to10->Draw();
  rp_kp_yCM_00to10->GetUpperRefYaxis()->SetRangeUser(-0.2,0.15);
  rp_kp_yCM_00to10->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_kp_yCM_00to10->SetLeftMargin(0.13);
  rp_kp_yCM_00to10->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_kp_yCM_00to10.png");
  canvas->Clear();
  delete rp_kp_yCM_00to10;

  TRatioPlot *rp_kp_yCM_10to40 = new TRatioPlot(h_vn_yCM_10to40_kp_norm, h_vn_yCM_10to40_kp_eff, "diffsig");
  rp_kp_yCM_10to40->Draw();
  rp_kp_yCM_10to40->GetUpperRefYaxis()->SetRangeUser(-0.2,0.15);
  rp_kp_yCM_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_kp_yCM_10to40->SetLeftMargin(0.13);
  rp_kp_yCM_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_kp_yCM_10to40.png");
  canvas->Clear();
  delete rp_kp_yCM_10to40;

  TRatioPlot *rp_kp_yCM_40to60 = new TRatioPlot(h_vn_yCM_40to60_kp_norm, h_vn_yCM_40to60_kp_eff, "diffsig");
  rp_kp_yCM_40to60->Draw();
  rp_kp_yCM_40to60->GetUpperRefYaxis()->SetRangeUser(-0.2,0.15);
  rp_kp_yCM_40to60->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_kp_yCM_40to60->SetLeftMargin(0.13);
  rp_kp_yCM_40to60->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_kp_yCM_40to60.png");
  canvas->Clear();
  delete rp_kp_yCM_40to60;

  TRatioPlot *rp_km_yCM_10to40 = new TRatioPlot(h_vn_yCM_10to40_km_norm, h_vn_yCM_10to40_km_eff, "diffsig");
  rp_km_yCM_10to40->Draw();
  rp_km_yCM_10to40->GetUpperRefYaxis()->SetRangeUser(-0.2,0.15);
  rp_km_yCM_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_km_yCM_10to40->SetLeftMargin(0.13);
  rp_km_yCM_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_km_yCM_10to40.png");
  canvas->Clear();
  delete rp_km_yCM_10to40;

  TRatioPlot *rp_pr_yCM_00to10 = new TRatioPlot(h_vn_yCM_00to10_pr_norm, h_vn_yCM_00to10_pr_eff, "diffsig");
  rp_pr_yCM_00to10->Draw();
  rp_pr_yCM_00to10->GetUpperRefYaxis()->SetRangeUser(-0.07,0.05);
  rp_pr_yCM_00to10->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_yCM_00to10->SetLeftMargin(0.13);
  rp_pr_yCM_00to10->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_yCM_00to10.png");
  canvas->Clear();
  delete rp_pr_yCM_00to10;

  TRatioPlot *rp_pr_yCM_10to40 = new TRatioPlot(h_vn_yCM_10to40_pr_norm, h_vn_yCM_10to40_pr_eff, "diffsig");
  rp_pr_yCM_10to40->Draw();
  rp_pr_yCM_10to40->GetUpperRefYaxis()->SetRangeUser(-0.07,0.05);
  rp_pr_yCM_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_yCM_10to40->SetLeftMargin(0.13);
  rp_pr_yCM_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_yCM_10to40.png");
  canvas->Clear();
  delete rp_pr_yCM_10to40;

  TRatioPlot *rp_pr_yCM_40to60 = new TRatioPlot(h_vn_yCM_40to60_pr_norm, h_vn_yCM_40to60_pr_eff, "diffsig");
  rp_pr_yCM_40to60->Draw();
  rp_pr_yCM_40to60->GetUpperRefYaxis()->SetRangeUser(-0.07,0.05);
  rp_pr_yCM_40to60->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_yCM_40to60->SetLeftMargin(0.13);
  rp_pr_yCM_40to60->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_yCM_40to60.png");
  canvas->Clear();
  delete rp_pr_yCM_40to60;

  TRatioPlot *rp_pr_yCM_00to10_symm = new TRatioPlot(h_vn_yCM_00to10_pr_norm_symm, h_vn_yCM_00to10_pr_eff_symm, "diffsig");
  rp_pr_yCM_00to10_symm->Draw();
  rp_pr_yCM_00to10_symm->GetUpperRefYaxis()->SetRangeUser(-0.09,0.06);
  rp_pr_yCM_00to10_symm->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_yCM_00to10_symm->SetLeftMargin(0.13);
  rp_pr_yCM_00to10_symm->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_yCM_00to10_symm.png");
  canvas->Clear();
  delete rp_pr_yCM_00to10_symm;

  TRatioPlot *rp_pr_yCM_10to40_symm = new TRatioPlot(h_vn_yCM_10to40_pr_norm_symm, h_vn_yCM_10to40_pr_eff_symm, "diffsig");
  rp_pr_yCM_10to40_symm->Draw();
  rp_pr_yCM_10to40_symm->GetUpperRefYaxis()->SetRangeUser(-0.09,0.06);
  rp_pr_yCM_10to40_symm->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_yCM_10to40_symm->SetLeftMargin(0.13);
  rp_pr_yCM_10to40_symm->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_yCM_10to40_symm.png");
  canvas->Clear();
  delete rp_pr_yCM_10to40_symm;

  TRatioPlot *rp_pr_yCM_40to60_symm = new TRatioPlot(h_vn_yCM_40to60_pr_norm_symm, h_vn_yCM_40to60_pr_eff_symm, "diffsig");
  rp_pr_yCM_40to60_symm->Draw();
  rp_pr_yCM_40to60_symm->GetUpperRefYaxis()->SetRangeUser(-0.09,0.06);
  rp_pr_yCM_40to60_symm->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_yCM_40to60_symm->SetLeftMargin(0.13);
  rp_pr_yCM_40to60_symm->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_yCM_40to60_symm.png");
  canvas->Clear();
  delete rp_pr_yCM_40to60_symm;




    TRatioPlot *rp_pp_pT_00to10 = new TRatioPlot(h_vn_pT_00to10_pp_norm, h_vn_pT_00to10_pp_eff, "diffsig");
  rp_pp_pT_00to10->Draw();
  rp_pp_pT_00to10->GetUpperRefYaxis()->SetRangeUser(-0.1,0.2);
  rp_pp_pT_00to10->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pp_pT_00to10->SetLeftMargin(0.13);
  rp_pp_pT_00to10->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pp_pT_00to10.png");
  canvas->Clear();
  delete rp_pp_pT_00to10;

  TRatioPlot *rp_pp_pT_10to40 = new TRatioPlot(h_vn_pT_10to40_pp_norm, h_vn_pT_10to40_pp_eff, "diffsig");
  rp_pp_pT_10to40->Draw();
  rp_pp_pT_10to40->GetUpperRefYaxis()->SetRangeUser(-0.1,0.2);
  rp_pp_pT_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pp_pT_10to40->SetLeftMargin(0.13);
  rp_pp_pT_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pp_pT_10to40.png");
  canvas->Clear();
  delete rp_pp_pT_10to40;

  TRatioPlot *rp_pp_pT_40to60 = new TRatioPlot(h_vn_pT_40to60_pp_norm, h_vn_pT_40to60_pp_eff, "diffsig");
  rp_pp_pT_40to60->Draw();
  rp_pp_pT_40to60->GetUpperRefYaxis()->SetRangeUser(-0.1,0.2);
  rp_pp_pT_40to60->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pp_pT_40to60->SetLeftMargin(0.13);
  rp_pp_pT_40to60->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pp_pT_40to60.png");
  canvas->Clear();
  delete rp_pp_pT_40to60;

  TRatioPlot *rp_pm_pT_00to10 = new TRatioPlot(h_vn_pT_00to10_pm_norm, h_vn_pT_00to10_pm_eff, "diffsig");
  rp_pm_pT_00to10->Draw();
  rp_pm_pT_00to10->GetUpperRefYaxis()->SetRangeUser(-0.12,0.12);
  rp_pm_pT_00to10->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pm_pT_00to10->SetLeftMargin(0.13);
  rp_pm_pT_00to10->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pm_pT_00to10.png");
  canvas->Clear();
  delete rp_pm_pT_00to10;

  TRatioPlot *rp_pm_pT_10to40 = new TRatioPlot(h_vn_pT_10to40_pm_norm, h_vn_pT_10to40_pm_eff, "diffsig");
  rp_pm_pT_10to40->Draw();
  rp_pm_pT_10to40->GetUpperRefYaxis()->SetRangeUser(-0.12,0.12);
  rp_pm_pT_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pm_pT_10to40->SetLeftMargin(0.13);
  rp_pm_pT_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pm_pT_10to40.png");
  canvas->Clear();
  delete rp_pm_pT_10to40;

  TRatioPlot *rp_pm_pT_40to60 = new TRatioPlot(h_vn_pT_40to60_pm_norm, h_vn_pT_40to60_pm_eff, "diffsig");
  rp_pm_pT_40to60->Draw();
  rp_pm_pT_40to60->GetUpperRefYaxis()->SetRangeUser(-0.12,0.12);
  rp_pm_pT_40to60->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pm_pT_40to60->SetLeftMargin(0.13);
  rp_pm_pT_40to60->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pm_pT_40to60.png");
  canvas->Clear();
  delete rp_pm_pT_40to60;

  TRatioPlot *rp_kp_pT_00to10 = new TRatioPlot(h_vn_pT_00to10_kp_norm, h_vn_pT_00to10_kp_eff, "diffsig");
  rp_kp_pT_00to10->Draw();
  rp_kp_pT_00to10->GetUpperRefYaxis()->SetRangeUser(-0.1,0.15);
  rp_kp_pT_00to10->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_kp_pT_00to10->SetLeftMargin(0.13);
  rp_kp_pT_00to10->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_kp_pT_00to10.png");
  canvas->Clear();
  delete rp_kp_pT_00to10;

  TRatioPlot *rp_kp_pT_10to40 = new TRatioPlot(h_vn_pT_10to40_kp_norm, h_vn_pT_10to40_kp_eff, "diffsig");
  rp_kp_pT_10to40->Draw();
  rp_kp_pT_10to40->GetUpperRefYaxis()->SetRangeUser(-0.1,0.15);
  rp_kp_pT_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_kp_pT_10to40->SetLeftMargin(0.13);
  rp_kp_pT_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_kp_pT_10to40.png");
  canvas->Clear();
  delete rp_kp_pT_10to40;

  TRatioPlot *rp_kp_pT_40to60 = new TRatioPlot(h_vn_pT_40to60_kp_norm, h_vn_pT_40to60_kp_eff, "diffsig");
  rp_kp_pT_40to60->Draw();
  rp_kp_pT_40to60->GetUpperRefYaxis()->SetRangeUser(-0.1,0.15);
  rp_kp_pT_40to60->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_kp_pT_40to60->SetLeftMargin(0.13);
  rp_kp_pT_40to60->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_kp_pT_40to60.png");
  canvas->Clear();
  delete rp_kp_pT_40to60;

  TRatioPlot *rp_km_pT_10to40 = new TRatioPlot(h_vn_pT_10to40_km_norm, h_vn_pT_10to40_km_eff, "diffsig");
  rp_km_pT_10to40->Draw();
  rp_km_pT_10to40->GetUpperRefYaxis()->SetRangeUser(-0.1,0.01);
  rp_km_pT_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_km_pT_10to40->SetLeftMargin(0.13);
  rp_km_pT_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_km_pT_10to40.png");
  canvas->Clear();
  delete rp_km_pT_10to40;

  TRatioPlot *rp_pr_pT_00to10 = new TRatioPlot(h_vn_pT_00to10_pr_norm, h_vn_pT_00to10_pr_eff, "diffsig");
  rp_pr_pT_00to10->Draw();
  rp_pr_pT_00to10->GetUpperRefYaxis()->SetRangeUser(-0.12,0.03);
  rp_pr_pT_00to10->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_pT_00to10->SetLeftMargin(0.13);
  rp_pr_pT_00to10->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_pT_00to10.png");
  canvas->Clear();
  delete rp_pr_pT_00to10;

  TRatioPlot *rp_pr_pT_10to40 = new TRatioPlot(h_vn_pT_10to40_pr_norm, h_vn_pT_10to40_pr_eff, "diffsig");
  rp_pr_pT_10to40->Draw();
  rp_pr_pT_10to40->GetUpperRefYaxis()->SetRangeUser(-0.12,0.03);
  rp_pr_pT_10to40->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_pT_10to40->SetLeftMargin(0.13);
  rp_pr_pT_10to40->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_pT_10to40.png");
  canvas->Clear();
  delete rp_pr_pT_10to40;

  TRatioPlot *rp_pr_pT_40to60 = new TRatioPlot(h_vn_pT_40to60_pr_norm, h_vn_pT_40to60_pr_eff, "diffsig");
  rp_pr_pT_40to60->Draw();
  rp_pr_pT_40to60->GetUpperRefYaxis()->SetRangeUser(-0.12,0.03);
  rp_pr_pT_40to60->GetLowerRefYaxis()->SetTitle("h_{1}-h_{2}/#deltah_{1}");
  rp_pr_pT_40to60->SetLeftMargin(0.13);
  rp_pr_pT_40to60->SetSeparationMargin(0.01);
  canvas->Update();
  canvas->SaveAs("ratioPlot_pr_pT_40to60.png");
  canvas->Clear();
  delete rp_pr_pT_40to60;

}
