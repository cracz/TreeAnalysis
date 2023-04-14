#include "PlotUtils.h"

void correlations(TString jobID, TString order_n_str)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {std::cout << "Wrong file!" << std::endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 1200);
  //canvas->SetGridx();
  //canvas->SetGridy();
  canvas->SetLogy(0);
  canvas->SetTicks();
  canvas->SetLeftMargin(0.17);
  canvas->SetRightMargin(0.05);
  canvas->cd();
  
  TProfile *p_EpdAEpdB = (TProfile*)file->Get("p_EpdAEpdB");
  TProfile *p_TpcAB    = (TProfile*)file->Get("p_TpcAB");
  TProfile *p_TpcAEpdA = (TProfile*)file->Get("p_TpcAEpdA");
  TProfile *p_TpcAEpdB = (TProfile*)file->Get("p_TpcAEpdB");
  TProfile *p_TpcBEpdA = (TProfile*)file->Get("p_TpcBEpdA");
  TProfile *p_TpcBEpdB = (TProfile*)file->Get("p_TpcBEpdB");

  TH1D *h_EpdAEpdB = p_EpdAEpdB->ProjectionX();
  TH1D *h_TpcAB    = p_TpcAB->ProjectionX();
  TH1D *h_TpcAEpdA = p_TpcAEpdA->ProjectionX();
  TH1D *h_TpcAEpdB = p_TpcAEpdB->ProjectionX();
  TH1D *h_TpcBEpdA = p_TpcBEpdA->ProjectionX();
  TH1D *h_TpcBEpdB = p_TpcBEpdB->ProjectionX();

  h_EpdAEpdB = PlotUtils::flipHisto(h_EpdAEpdB);
  h_TpcAB    = PlotUtils::flipHisto(h_TpcAB);
  h_TpcAEpdA = PlotUtils::flipHisto(h_TpcAEpdA);
  h_TpcAEpdB = PlotUtils::flipHisto(h_TpcAEpdB);
  h_TpcBEpdA = PlotUtils::flipHisto(h_TpcBEpdA);
  h_TpcBEpdB = PlotUtils::flipHisto(h_TpcBEpdB);

  h_EpdAEpdB = PlotUtils::trimCentralityPlot(h_EpdAEpdB);
  h_TpcAB    = PlotUtils::trimCentralityPlot(h_TpcAB);
  h_TpcAEpdA = PlotUtils::trimCentralityPlot(h_TpcAEpdA);
  h_TpcAEpdB = PlotUtils::trimCentralityPlot(h_TpcAEpdB);
  h_TpcBEpdA = PlotUtils::trimCentralityPlot(h_TpcBEpdA);
  h_TpcBEpdB = PlotUtils::trimCentralityPlot(h_TpcBEpdB);

  h_EpdAEpdB->SetTitleOffset(2.2, "y");
  h_TpcAB->SetTitleOffset(2.2, "y");
  h_TpcAEpdA->SetTitleOffset(2.2, "y");
  h_TpcAEpdB->SetTitleOffset(2.2, "y");
  h_TpcBEpdA->SetTitleOffset(2.2, "y");
  h_TpcBEpdB->SetTitleOffset(2.2, "y");

  h_EpdAEpdB->SetMarkerStyle(20);
  h_TpcAB->SetMarkerStyle(20);
  h_TpcAEpdA->SetMarkerStyle(20);
  h_TpcAEpdB->SetMarkerStyle(20);
  h_TpcBEpdA->SetMarkerStyle(20);
  h_TpcBEpdB->SetMarkerStyle(20);

  h_EpdAEpdB->SetMarkerSize(2.5);
  h_TpcAB->SetMarkerSize(2.5);
  h_TpcAEpdA->SetMarkerSize(2.5);
  h_TpcAEpdB->SetMarkerSize(2.5);
  h_TpcBEpdA->SetMarkerSize(2.5);
  h_TpcBEpdB->SetMarkerSize(2.5);

  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);

  TLine *zeroLine = new TLine(0., 0., 60., 0.);
  zeroLine->SetLineStyle(9);
  zeroLine->SetLineWidth(3);
  
  h_EpdAEpdB->Draw("E1P");
  zeroLine->Draw("same");
  canvas->SaveAs(jobID + "_correlEpdAEpdB.png");
  canvas->Clear();
  /*
  h_TpcAB->Draw("E1P");
  canvas->SaveAs(jobID + "_correlTpcAB.png");
  canvas->Clear();

  h_TpcAEpdA->Draw("E1P");
  canvas->SaveAs(jobID + "_correlTpcAEpdA.png");
  canvas->Clear();

  h_TpcAEpdB->Draw("E1P");
  canvas->SaveAs(jobID + "_correlTpcAEpdB.png");
  canvas->Clear();
  */
  h_TpcBEpdA->Draw("E1P");
  zeroLine->Draw("same");
  canvas->SaveAs(jobID + "_correlTpcBEpdA.png");
  canvas->Clear();

  h_TpcBEpdB->Draw("E1P");
  zeroLine->Draw("same");
  canvas->SaveAs(jobID + "_correlTpcBEpdB.png");
  canvas->Clear();

}
