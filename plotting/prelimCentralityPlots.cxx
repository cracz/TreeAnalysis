#include "PlotUtils.h"

void prelimCentralityPlots()
{
  TString fileName = "Normal.picoDst.result.combined.root";
  TFile *file = TFile::Open(fileName);
  if(!file) {std::cout << "Wrong file!" << std::endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1000, 1000);//1200, 1200);
  //canvas->SetGridx();
  //canvas->SetGridy();
  canvas->SetLogy(0);
  canvas->SetTicks();
  canvas->SetTopMargin(0.04);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.1);
  canvas->SetLeftMargin(0.15);
  canvas->cd();
  
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetLineWidth(3);


  TProfile *p_vn_pp   = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm   = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_pr   = (TProfile*)file->Get("p_vn_pr");

  // Convert profiles to histograms
  TH1D *h_vn_pp = p_vn_pp->ProjectionX();
  TH1D *h_vn_pm = p_vn_pm->ProjectionX();
  TH1D *h_vn_pr = p_vn_pr->ProjectionX();

  // Flip centrality plots
  h_vn_pp = PlotUtils::flipHisto(h_vn_pp);
  h_vn_pm = PlotUtils::flipHisto(h_vn_pm);
  h_vn_pr = PlotUtils::flipHisto(h_vn_pr);

  // Trim and clean up x-axis
  h_vn_pp = PlotUtils::trimCentralityPlot(h_vn_pp);
  h_vn_pm = PlotUtils::trimCentralityPlot(h_vn_pm);
  h_vn_pr = PlotUtils::trimCentralityPlot(h_vn_pr);

  // Retrieve systematic uncertainties
  TFile* systematicFile = TFile::Open("systematicErrors.root", "READ");
  TGraphErrors* sys_pp = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pp_Normal_flip");
  TGraphErrors* sys_pm = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pm_Normal_flip");
  TGraphErrors* sys_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pr_Normal_flip");
  ////
  

  // Set various aesthetics
  sys_pp->SetMarkerColor(1);
  sys_pp->SetLineColor(1);

  sys_pm->SetMarkerColor(4);
  sys_pm->SetLineColor(4);

  //sys_pr should be fine as it is
  
  h_vn_pp->SetMarkerStyle(21);
  h_vn_pp->SetMarkerSize(2.5);
  h_vn_pp->SetMarkerColor(1);
  h_vn_pp->SetLineColor(1);
  h_vn_pp->SetLineWidth(3);
  h_vn_pp->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pm->SetMarkerStyle(22);
  h_vn_pm->SetMarkerSize(3);
  h_vn_pm->SetMarkerColor(4);
  h_vn_pm->SetLineColor(4);
  h_vn_pm->SetLineWidth(3);
  h_vn_pm->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pr->SetMarkerStyle(20);
  h_vn_pr->SetMarkerSize(2.5);
  h_vn_pr->SetMarkerColor(kRed-4);
  h_vn_pr->SetLineColor(kRed-4);
  h_vn_pr->SetLineWidth(3);
  h_vn_pr->GetYaxis()->SetTitleOffset(1.7);
  ////

  THStack *allCentralityStack = new THStack("allCentralityStack", ";Centrality (%);v_{3} {#psi_{1}}");
  allCentralityStack->Add(h_vn_pr);
  allCentralityStack->Add(h_vn_pp);
  allCentralityStack->Add(h_vn_pm);
  
  // Make text boxes, legends, and line at zero
  TPaveText* prelimText = new TPaveText(15, 0.01, 45, 0.014, "NB");
  prelimText->AddText("STAR Preliminary");
  prelimText->SetTextColor(kRed);
  prelimText->SetFillColorAlpha(0,0);
  prelimText->SetTextSize(0.04);

  TPaveText *allText = new TPaveText(15, 0.016, 45, 0.027, "NB");
  allText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT (year 2018)");
  allText->AddText("0 < y_{CM} < 0.5 GeV");
  //allText->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
  allText->SetFillColorAlpha(0,0);
  allText->SetLineColorAlpha(0,0);
  allText->SetTextSize(0.045);

  TLegend *allLegend = new TLegend(0.28, 0.14, 0.6, 0.3);
  allLegend->AddEntry(h_vn_pp,"#pi^{+}, 0.18 #leq p_{T} #leq 1.6 GeV");
  allLegend->AddEntry(h_vn_pm,"#pi^{-}, 0.18 #leq p_{T} #leq 1.6 GeV");
  allLegend->AddEntry(h_vn_pr,"p, 0.4 #leq p_{T} #leq 2.0 GeV");
  allLegend->SetFillColorAlpha(0,0);
  allLegend->SetLineColorAlpha(0,0);
  allLegend->SetTextSize(0.04);

  TLine *zeroLine = new TLine(0, 0, 60, 0);
  zeroLine->SetLineStyle(9);
  zeroLine->SetLineWidth(3);
  ////

  allCentralityStack->Draw();
  allCentralityStack->GetXaxis()->SetLabelSize(0.045);
  allCentralityStack->GetYaxis()->SetLabelSize(0.045);
  allCentralityStack->GetXaxis()->SetTitleOffset(1.0);
  allCentralityStack->GetYaxis()->SetTitleOffset(1.4);
  allCentralityStack->GetXaxis()->SetTitleSize(0.045);
  allCentralityStack->GetYaxis()->SetTitleSize(0.05);
  allCentralityStack->GetXaxis()->SetNdivisions(210);
  allCentralityStack->SetMaximum(0.03);
  allCentralityStack->SetMinimum(-0.045);
  allCentralityStack->Draw("NOSTACK E1P");
  zeroLine->Draw("SAME");
  allCentralityStack->Draw("NOSTACK E1P SAME");
  sys_pp->Draw("[]");
  sys_pm->Draw("[]");
  sys_pr->Draw("[]");
  allLegend->Draw();
  allText->Draw();
  prelimText->Draw();
  canvas->SaveAs("v3_allCentralityStack.png");
  canvas->Clear();

  systematicFile->Close();
  file->Close();
}
