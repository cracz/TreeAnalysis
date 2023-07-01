#include "PlotUtils.h"

void resolutionPlot()
{
  //TFile* file = TFile::Open("eventPlaneSystematics_3p9GeV_1to5_max10.root","READ");
  TFile* file = TFile::Open("eventPlaneSystematics_3p0GeV_v2.root","READ");
  

  TH1D* h_resolutionsWithStats = (TH1D*)file->Get("h_resolutionsWithStats");
  TH1D* h_resolutionswithSysts = (TH1D*)file->Get("h_resolutionsWithSysts");
  //h_resolutionswithSysts = PlotUtils::trimCentralityPlot(h_resolutionswithSysts);
  
  h_resolutionsWithStats->SetTitle("");
  h_resolutionsWithStats->GetYaxis()->SetTitle("R_{31}");
  h_resolutionsWithStats->GetYaxis()->SetTitleSize(0.055);
  h_resolutionsWithStats->GetYaxis()->SetTitleOffset(1.15);
  h_resolutionsWithStats->GetYaxis()->SetLabelSize(0.045);
  h_resolutionsWithStats->GetYaxis()->SetRangeUser(0.0, 1.0);
  h_resolutionsWithStats->GetXaxis()->SetTitle("Centrality (%)");

  //h_resolutionsWithStats->GetXaxis()->SetRangeUser(0.0, 40.0);
  h_resolutionsWithStats->GetXaxis()->SetRangeUser(0.0, 60.0);
  
  h_resolutionsWithStats->GetXaxis()->SetTitleSize(0.045);
  h_resolutionsWithStats->GetXaxis()->SetTitleOffset(1.1);
  h_resolutionsWithStats->GetXaxis()->SetLabelSize(0.05);
  h_resolutionsWithStats->SetLineWidth(2);
  h_resolutionsWithStats->SetLineColor(kBlack);
  h_resolutionsWithStats->SetMarkerStyle(20);
  h_resolutionsWithStats->SetMarkerSize(2);
  h_resolutionsWithStats->SetMarkerColor(kBlue);

  TGraphErrors* finalGraph = new TGraphErrors(h_resolutionswithSysts);
  finalGraph->SetTitle("");
  finalGraph->GetYaxis()->SetTitle("R_{31}");
  finalGraph->GetYaxis()->SetTitleSize(0.055);
  finalGraph->GetYaxis()->SetTitleOffset(1.15);
  finalGraph->GetYaxis()->SetLabelSize(0.045);
  finalGraph->GetYaxis()->SetRangeUser(0.0, 1.0);
  finalGraph->GetXaxis()->SetTitle("Centrality (%)");
  
  //finalGraph->GetXaxis()->SetRangeUser(0.0, 40.0);
  finalGraph->GetXaxis()->SetRangeUser(0.0, 60.0);

  finalGraph->GetXaxis()->SetTitleSize(0.045);
  finalGraph->GetXaxis()->SetTitleOffset(1.1);
  finalGraph->GetXaxis()->SetLabelSize(0.05);
  finalGraph->SetLineWidth(2);
  finalGraph->SetLineColor(kBlack);
  finalGraph->SetMarkerStyle(20);
  finalGraph->SetMarkerSize(2);
  finalGraph->SetMarkerColor(kBlue);

  // Remove x errors
  for (int i = 0; i < finalGraph->GetN(); i++)
    finalGraph->SetPointError(i, 0.0, finalGraph->GetErrorY(i));
  
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 1000);
  canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetTopMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.13);
  canvas->cd();

  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetErrorX(0);

  h_resolutionsWithStats->Draw("EP");
  canvas->Update();
  finalGraph->Draw("[]");
  canvas->SaveAs("resolutionPlot.pdf");

  delete finalGraph;
  delete canvas;
}
