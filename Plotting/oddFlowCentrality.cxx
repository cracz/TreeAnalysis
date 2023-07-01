#include "PlotUtils.h"
#include "TProfileHelper.h"

void oddFlowCentrality(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {std::cout << "Wrong file!" << std::endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1000, 1000);
  //canvas->SetGridx();
  //canvas->SetGridy();
  canvas->SetLogy(0);
  canvas->SetTicks();
  canvas->SetTopMargin(0.04);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.1);
  canvas->SetLeftMargin(0.14);
  canvas->cd();
  
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetLineWidth(3);
  //gStyle->SetOptDate();

  TProfile2D *p2_vn_yCM_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr_symmetry");  
  
  TH1D *h_oddVsCent = new TH1D("h_oddVsCent", ";Centrality (%);v_{3} {#Psi_{1}}", 16, 0, 16);

  TProfile* p1;
  TString name;
  
  for (int centBin = 5; centBin <= 16; centBin++)
    {
      name.Form("p1_%d", centBin);
      p1 = p2_vn_yCM_cent_pr_symmetry->ProfileY(name, centBin, centBin);

      TProfileHelper ph;
      ph.initialize(p1);
      ph.removeEvenComponent();
      
      h_oddVsCent->SetBinContent(centBin, ph.oddAverage);
      h_oddVsCent->SetBinError(centBin, ph.oddAverageUncertainty);
    }


  h_oddVsCent = PlotUtils::flipHisto(h_oddVsCent);
  h_oddVsCent = PlotUtils::trimCentralityPlot(h_oddVsCent);
  h_oddVsCent->SetMarkerStyle(20);
  h_oddVsCent->SetMarkerSize(2.5);
  h_oddVsCent->SetMarkerColor(kRed-4);
  h_oddVsCent->SetLineColor(kRed-4);
  h_oddVsCent->SetLineWidth(3);
  h_oddVsCent->GetYaxis()->SetTitleOffset(1.7);

  TLine *zeroLine = new TLine(0, 0, 60, 0);
  zeroLine->SetLineStyle(9);
  zeroLine->SetLineWidth(3);
  /*
  TProfile* p_vn_pr = (TProfile*)file->Get("p_vn_pr");
  TH1D* h_vn_pr = p_vn_pr->ProjectionX();
  h_vn_pr = PlotUtils::flipHisto(h_vn_pr);
  h_vn_pr = PlotUtils::trimCentralityPlot(h_vn_pr);
  h_vn_pr->SetMarkerStyle(20);
  h_vn_pr->SetMarkerSize(2.5);
  h_vn_pr->SetMarkerColor(kBlue);
  h_vn_pr->SetLineColor(kBlue);
  h_vn_pr->SetLineWidth(3);
  */
  
  h_oddVsCent->GetXaxis()->SetLabelSize(0.045);
  h_oddVsCent->GetYaxis()->SetLabelSize(0.045);
  h_oddVsCent->GetXaxis()->SetTitleOffset(1.0);
  h_oddVsCent->GetYaxis()->SetTitleOffset(1.4);
  h_oddVsCent->GetXaxis()->SetTitleSize(0.045);
  h_oddVsCent->GetYaxis()->SetTitleSize(0.05);
  h_oddVsCent->GetXaxis()->SetNdivisions(210);
  h_oddVsCent->SetMaximum(0.02);
  h_oddVsCent->SetMinimum(-0.06);
  h_oddVsCent->Draw();
  zeroLine->Draw("SAME");
  h_oddVsCent->Draw("EP SAME");
  //h_vn_pr->Draw("EP SAME");
  canvas->SaveAs("v3_oddVsCentrality.pdf");
  canvas->Clear();

  file->Close();
}
