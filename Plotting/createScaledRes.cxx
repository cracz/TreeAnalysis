#include "PlotUtils.h"

Double_t weightedAvgPeakResNoFlip(TH1D* resolutions)
{
  Double_t content8 = resolutions->GetBinContent(8);
  Double_t content9 = resolutions->GetBinContent(9);
  Double_t content10 = resolutions->GetBinContent(10);
  Double_t content11 = resolutions->GetBinContent(11);
  Double_t content12 = resolutions->GetBinContent(12);

  Double_t entries8 = 1.0 / resolutions->GetBinError(8);
  Double_t entries9 = 1.0 / resolutions->GetBinError(9);
  Double_t entries10 = 1.0 / resolutions->GetBinError(10);
  Double_t entries11 = 1.0 / resolutions->GetBinError(11);
  Double_t entries12 = 1.0 / resolutions->GetBinError(12);

  Double_t numerator = entries8 * content8 +
    entries9 * content9 +
    entries10 * content10 +
    entries11 * content11 +
    entries12 * content12;

  Double_t denominator = entries8 +
    entries9 +
    entries10 +
    entries12 +
    entries12;

  /*
  Double_t numerator = content5 +
    content6 +
    content7 +
    content8 +
    content9;
  
  Double_t denominator = 5.0;
  */
  return numerator/denominator;
}



void createScaledRes(TString fileName_9to16  = "resolutionInfo_INPUT_3p0GeV_EPDB9to16.root",
		     TString fileName_10to16 = "resolutionInfo_INPUT_3p0GeV_EPDB10to16.root",
		     TString fileName_13to16 = "resolutionInfo_INPUT_3p0GeV_EPDB13to16.root")
{
  TFile *file_9to16 = TFile::Open(fileName_9to16);
  if(!file_9to16) {std::cout << "Wrong 9to16 file!" << std::endl; return;}

  TFile *file_10to16 = TFile::Open(fileName_10to16);
  if(!file_10to16) {std::cout << "Wrong 10to16 file!" << std::endl; return;}

  TFile *file_13to16 = TFile::Open(fileName_13to16);
  if(!file_13to16) {std::cout << "Wrong 13to16 file!" << std::endl; return;}

  TH1D* h_9to16 = (TH1D*)file_9to16->Get("h_resolutions");
  h_9to16->SetName("h_resolutions_9to16");
  TH1D* h_10to16 = (TH1D*)file_10to16->Get("h_resolutions");
  h_10to16->SetName("h_resolutions_10to16");
  TH1D* h_13to16 = (TH1D*)file_13to16->Get("h_resolutions");
  h_13to16->SetName("h_resolutions_13to16");

  TH1D* h_9to16_scaled = (TH1D*)h_9to16->Clone("h_9to16_scaled");
  h_9to16_scaled->SetName("h_resolutions");
  h_9to16_scaled->Scale(weightedAvgPeakResNoFlip(h_13to16) / weightedAvgPeakResNoFlip(h_9to16));

  
  TFile* newFile = TFile::Open("resolutionInfo_INPUT_3p0GeV_9to16scaled.root", "RECREATE");
  h_9to16_scaled->Write();
  newFile->Close();

  std::cout << "Wrote scaled resolutions to resolutionInfo_INPUT_3p0GeV_9to16scaled.root." << std::endl;


  TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 800);
  canvas->SetTicks();
  canvas->SetLeftMargin(0.17);
  canvas->SetRightMargin(0.05);
  canvas->cd();

  h_9to16  = PlotUtils::flipHisto(h_9to16);
  h_10to16 = PlotUtils::flipHisto(h_10to16);
  h_13to16 = PlotUtils::flipHisto(h_13to16);
  h_9to16_scaled = PlotUtils::flipHisto(h_9to16_scaled);

  h_9to16  = PlotUtils::trimCentralityPlot(h_9to16);
  h_10to16 = PlotUtils::trimCentralityPlot(h_10to16);
  h_13to16 = PlotUtils::trimCentralityPlot(h_13to16);
  h_9to16_scaled = PlotUtils::trimCentralityPlot(h_9to16_scaled);

  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);


  TLine *zeroLine = new TLine(0., 0., 60., 0.);
  zeroLine->SetLineStyle(9);
  zeroLine->SetLineWidth(3);


  h_9to16->SetMarkerStyle(20);
  h_10to16->SetMarkerStyle(20);
  h_13to16->SetMarkerStyle(20);
  h_9to16_scaled->SetMarkerStyle(20);


  h_9to16->SetMarkerSize(1.5);
  h_10to16->SetMarkerSize(1.5);
  h_13to16->SetMarkerSize(1.5);
  h_9to16_scaled->SetMarkerSize(1.5);


  h_10to16->SetMarkerColor(kRed);
  h_13to16->SetMarkerColor(kGreen+3);
  h_9to16_scaled->SetMarkerColor(kBlue);

  h_10to16->SetLineColor(kRed);
  h_13to16->SetLineColor(kGreen+3);
  h_9to16_scaled->SetLineColor(kBlue);


  TLegend* legend = new TLegend(0.4, 0.15, 0.55, 0.4);
  //legend->SetHeader("EPD B Rings","C");
  legend->SetHeader("EPD B Rings");
  legend->AddEntry(h_9to16, "9 - 16");
  legend->AddEntry(h_10to16, "10 - 16");
  legend->AddEntry(h_9to16_scaled, "9 - 16, scaled");
  legend->SetFillColorAlpha(0,0);
  legend->SetLineColorAlpha(0,0);
  legend->SetTextSize(0.04);

  h_9to16->SetMinimum(0.0);
  h_9to16->GetYaxis()->SetTitleOffset(1.7);

  h_9to16->Draw("E1P");
  h_10to16->Draw("E1PSAME");
  h_9to16_scaled->Draw("E1PSAME");
  legend->Draw();

  canvas->SaveAs("resolutionsAt3p0GeV.png");
  delete canvas;
}
