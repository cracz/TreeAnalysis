#include "PlotUtils.h"

void vnVsYscan(TString jobID, TString order_n_str)
{
  //TH1::SetDefaultSumw2();
  
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

  
  TProfile2D *p2_vn_pT_vs_yCM_pr = (TProfile2D*)file->Get("p2_vn_pT_vs_yCM_pr");

  // Loop over y axis bins starts here??
  
  TProfile *p_vn_yCM_pr_1 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_1", 1, 2);
  TProfile *p_vn_yCM_pr_2 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_2", 3, 4);
  TProfile *p_vn_yCM_pr_3 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_3", 5, 6);
  TProfile *p_vn_yCM_pr_4 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_4", 7, 10);
  //TProfile *p_vn_yCM_pr_4 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_4", 7, 8);
  //TProfile *p_vn_yCM_pr_5 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_5", 9, 10);

  /*
  TProfile *p_vn_yCM_pr_bin6 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_bin6", 6, 6);
  TProfile *p_vn_yCM_pr_bin7 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_bin7", 7, 7);
  TProfile *p_vn_yCM_pr_bin8 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_bin8", 8, 8);
  TProfile *p_vn_yCM_pr_bin9 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_bin9", 9, 9);
  TProfile *p_vn_yCM_pr_bin10 = p2_vn_pT_vs_yCM_pr->ProfileX("p_vn_yCM_pr_bin10", 10, 10);
  */

  TH1D *h_vn_yCM_pr_1 = new TH1D("h_vn_yCM_pr_1", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_pr_2 = new TH1D("h_vn_yCM_pr_2", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_pr_3 = new TH1D("h_vn_yCM_pr_3", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_pr_4 = new TH1D("h_vn_yCM_pr_4", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  //TH1D *h_vn_yCM_pr_5 = new TH1D("h_vn_yCM_pr_5", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  
  h_vn_yCM_pr_1 = p_vn_yCM_pr_1->ProjectionX();
  h_vn_yCM_pr_2 = p_vn_yCM_pr_2->ProjectionX();
  h_vn_yCM_pr_3 = p_vn_yCM_pr_3->ProjectionX();
  h_vn_yCM_pr_4 = p_vn_yCM_pr_4->ProjectionX();
  //h_vn_yCM_pr_5 = p_vn_yCM_pr_5->ProjectionX();  
  
  h_vn_yCM_pr_1->SetTitle(";y-y_{mid};v_{3}");
  h_vn_yCM_pr_2->SetTitle(";y-y_{mid};v_{3}");
  h_vn_yCM_pr_3->SetTitle(";y-y_{mid};v_{3}");
  h_vn_yCM_pr_4->SetTitle(";y-y_{mid};v_{3}");
  //h_vn_yCM_pr_5->SetTitle(";y-y_{mid};v_{3}");

  h_vn_yCM_pr_1->SetMarkerStyle(20);
  h_vn_yCM_pr_2->SetMarkerStyle(33);
  h_vn_yCM_pr_3->SetMarkerStyle(34);
  h_vn_yCM_pr_4->SetMarkerStyle(20);
  //h_vn_yCM_pr_5->SetMarkerStyle(20);

  h_vn_yCM_pr_1->SetMarkerSize(2);
  h_vn_yCM_pr_2->SetMarkerSize(4);
  h_vn_yCM_pr_3->SetMarkerSize(3);
  h_vn_yCM_pr_4->SetMarkerSize(3);

  h_vn_yCM_pr_1->SetLineWidth(3);
  h_vn_yCM_pr_2->SetLineWidth(3);
  h_vn_yCM_pr_3->SetLineWidth(3);
  h_vn_yCM_pr_4->SetLineWidth(3);

  h_vn_yCM_pr_1->SetMarkerColor(1);
  h_vn_yCM_pr_2->SetMarkerColor(kOrange+7);
  h_vn_yCM_pr_3->SetMarkerColor(kGreen+3);
  h_vn_yCM_pr_4->SetMarkerColor(kGreen+1);

  h_vn_yCM_pr_1->SetLineColor(1);
  h_vn_yCM_pr_2->SetLineColor(kOrange+7);
  h_vn_yCM_pr_3->SetLineColor(kGreen+3);
  h_vn_yCM_pr_4->SetLineColor(kGreen+1);

  
  // OMIT BAD SPOTS
  h_vn_yCM_pr_1->SetBinContent(7,0.0);
  h_vn_yCM_pr_1->SetBinError(7,0.0);
  h_vn_yCM_pr_1->SetBinContent(8,0.0);
  h_vn_yCM_pr_1->SetBinError(8,0.0);
  h_vn_yCM_pr_1->SetBinContent(9,0.0);
  h_vn_yCM_pr_1->SetBinError(9,0.0);
  h_vn_yCM_pr_1->SetBinContent(10,0.0);
  h_vn_yCM_pr_1->SetBinError(10,0.0);

  h_vn_yCM_pr_2->SetBinContent(3,0.0);
  h_vn_yCM_pr_2->SetBinError(3,0.0);
  h_vn_yCM_pr_2->SetBinContent(4,0.0);
  h_vn_yCM_pr_2->SetBinError(4,0.0);
  h_vn_yCM_pr_2->SetBinContent(5,0.0);
  h_vn_yCM_pr_2->SetBinError(5,0.0);
  h_vn_yCM_pr_2->SetBinContent(6,0.0);
  h_vn_yCM_pr_2->SetBinError(6,0.0);
  h_vn_yCM_pr_2->SetBinContent(7,0.0);
  h_vn_yCM_pr_2->SetBinError(7,0.0);

  h_vn_yCM_pr_3->SetBinContent(1,0.0);
  h_vn_yCM_pr_3->SetBinError(1,0.0);
  h_vn_yCM_pr_3->SetBinContent(2,0.0);
  h_vn_yCM_pr_3->SetBinError(2,0.0);
  h_vn_yCM_pr_3->SetBinContent(3,0.0);
  h_vn_yCM_pr_3->SetBinError(3,0.0);

  h_vn_yCM_pr_4->SetBinContent(1,0.0);
  h_vn_yCM_pr_4->SetBinError(1,0.0);
  h_vn_yCM_pr_4->SetBinContent(2,0.0);
  h_vn_yCM_pr_4->SetBinError(2,0.0);
  h_vn_yCM_pr_4->SetBinContent(3,0.0);
  h_vn_yCM_pr_4->SetBinError(3,0.0);
  ////

  Double_t max = 0.04;//0.08;
  Double_t min = -0.06;//-0.09;
  
  h_vn_yCM_pr_1->SetMaximum(max);
  h_vn_yCM_pr_1->SetMinimum(min);
  h_vn_yCM_pr_1->Draw("E1P");
  //canvas->SaveAs("h_vn_yCM_pr_1.png");
  canvas->Clear();

  h_vn_yCM_pr_2->SetMaximum(max);
  h_vn_yCM_pr_2->SetMinimum(min);
  h_vn_yCM_pr_2->Draw("E1P");
  //canvas->SaveAs("h_vn_yCM_pr_2.png");
  canvas->Clear();

  h_vn_yCM_pr_3->SetMaximum(max);
  h_vn_yCM_pr_3->SetMinimum(min);
  h_vn_yCM_pr_3->Draw("E1P");
  //canvas->SaveAs("h_vn_yCM_pr_3.png");
  canvas->Clear();

  h_vn_yCM_pr_4->SetMaximum(max);
  h_vn_yCM_pr_4->SetMinimum(min);
  h_vn_yCM_pr_4->Draw("E1P");
  //canvas->SaveAs("h_vn_yCM_pr_4.png");
  canvas->Clear();

  /*
  h_vn_yCM_pr_5->SetMaximum(max);
  h_vn_yCM_pr_5->SetMinimum(min);
  h_vn_yCM_pr_5->Draw("E1P");
  canvas->SaveAs("h_vn_yCM_pr_5.png");
  canvas->Clear();
  */

  TFile* systematicFile = TFile::Open("systematicErrors.root", "READ");
  TGraphErrors* sys_yCM_pr_1 = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_pr_1_px");
  TGraphErrors* sys_yCM_pr_2 = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_pr_2_px");
  TGraphErrors* sys_yCM_pr_3 = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_yCM_pr_3_px");

  //sys_yCM_pr_1->SetMarkerStyle(20);
  //sys_yCM_pr_2->SetMarkerStyle(20);
  //sys_yCM_pr_3->SetMarkerStyle(20);

  sys_yCM_pr_1->SetMarkerColor(1);
  sys_yCM_pr_2->SetMarkerColor(kOrange+7);
  sys_yCM_pr_3->SetMarkerColor(kGreen+3);

  sys_yCM_pr_1->SetLineColor(1);
  sys_yCM_pr_2->SetLineColor(kOrange+7);
  sys_yCM_pr_3->SetLineColor(kGreen+3);
  
  TLegend *prLegend_symm = new TLegend(0.19, 0.15, 0.39, 0.3);
  prLegend_symm->AddEntry(h_vn_yCM_pr_1, "0 #leq p_{T} < 0.5");
  prLegend_symm->AddEntry(h_vn_yCM_pr_2, "0.5 #leq p_{T} < 1.0");
  prLegend_symm->AddEntry(h_vn_yCM_pr_3, "1.0 #leq p_{T} < 1.5");
  prLegend_symm->SetBorderSize(0);
  prLegend_symm->SetFillColorAlpha(0,0);
  prLegend_symm->SetTextSize(0.04);

  //TPaveText *prText_y_symm = new TPaveText(-0.4, 0.05, 0.8, 0.075, "NB");
  TPaveText *prText_y_symm = new TPaveText(-0.4, 0.025, 0.8, 0.04, "NB");
  prText_y_symm->AddText("Proton");
  prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT (year 2018)");
  prText_y_symm->AddText("0 - 60% Centrality");
  prText_y_symm->SetFillColorAlpha(0,0);
  prText_y_symm->SetLineColorAlpha(0,0);
  prText_y_symm->SetTextSize(.035);

  //TPaveText* prelimText_symm = new TPaveText(-0.17, 0.035, 0.63, 0.055, "NB");
  TPaveText* prelimText_symm = new TPaveText(-0.2, 0.02, 0.6, 0.025, "NB");
  prelimText_symm->AddText("STAR Preliminary");
  prelimText_symm->SetTextColor(kRed);
  prelimText_symm->SetFillColorAlpha(0,0);
  prelimText_symm->SetTextSize(0.04);
  
  TLine *zeroLine_y_pr = new TLine(-1, 0, 1, 0);
  zeroLine_y_pr->SetLineStyle(9);
  zeroLine_y_pr->SetLineWidth(3);

  THStack* stack = new THStack("stack", ";y-y_{mid};v_{3} {#psi_{1} EP}");
  stack->Add(h_vn_yCM_pr_3);
  stack->Add(h_vn_yCM_pr_2);
  stack->Add(h_vn_yCM_pr_1);


  stack->Draw();
  stack->GetYaxis()->SetLabelSize(0.043);
  stack->GetXaxis()->SetLabelSize(0.043);
  stack->GetYaxis()->SetTitleOffset(1.4);
  stack->GetXaxis()->SetTitleOffset(1.0);
  stack->GetXaxis()->SetNdivisions(210);
  stack->GetXaxis()->SetTitleSize(0.045);
  stack->GetYaxis()->SetTitleSize(0.05);
  stack->Draw("E1P NOSTACK");
  zeroLine_y_pr->Draw("SAME");
  stack->Draw("E1P NOSTACK SAME");
  sys_yCM_pr_1->Draw("[]");
  sys_yCM_pr_2->Draw("[]");
  sys_yCM_pr_3->Draw("[]");
  prLegend_symm->Draw();
  prText_y_symm->Draw();
  prelimText_symm->Draw();
  /*
  h_vn_yCM_pr_3->Draw("e1p");
  h_vn_yCM_pr_2->Draw("e1p same");
  h_vn_yCM_pr_1->Draw("e1p same");
  */
  canvas->SaveAs("v3_yCM_pr.png");
}
