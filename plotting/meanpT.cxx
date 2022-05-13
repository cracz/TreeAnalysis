#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TStyle.h"
#include "PlotUtils.h"

void meanpT(TString jobID, TString order_n_str)
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
  //canvas->SetLeftMargin(0.15);
  canvas->cd();

  gStyle->SetOptStat(0);


  TProfile *p_meanpT_vs_yCM_pp = (TProfile*)file->Get("p_meanpT_vs_yCM_pp");
  TProfile *p_meanpT_vs_yCM_pm = (TProfile*)file->Get("p_meanpT_vs_yCM_pm");
  TProfile *p_meanpT_vs_yCM_kp = (TProfile*)file->Get("p_meanpT_vs_yCM_kp");
  TProfile *p_meanpT_vs_yCM_km = (TProfile*)file->Get("p_meanpT_vs_yCM_km");
  TProfile *p_meanpT_vs_yCM_pr = (TProfile*)file->Get("p_meanpT_vs_yCM_pr");
  TProfile *p_meanpT_vs_yCM_pr_alt = (TProfile*)file->Get("p_meanpT_vs_yCM_pr_alt");
  TProfile *p_meanpT_vs_yCM_de = (TProfile*)file->Get("p_meanpT_vs_yCM_de");
  TProfile *p_meanpT_vs_yCM_tr = (TProfile*)file->Get("p_meanpT_vs_yCM_tr");

  p_meanpT_vs_yCM_kp->Rebin();
  p_meanpT_vs_yCM_km->Rebin();

  // Convert profiles to histograms
  TH1D *h_meanpT_vs_yCM_pp = p_meanpT_vs_yCM_pp->ProjectionX();
  TH1D *h_meanpT_vs_yCM_pm = p_meanpT_vs_yCM_pm->ProjectionX();
  TH1D *h_meanpT_vs_yCM_kp = p_meanpT_vs_yCM_kp->ProjectionX();
  TH1D *h_meanpT_vs_yCM_km = p_meanpT_vs_yCM_km->ProjectionX();
  TH1D *h_meanpT_vs_yCM_pr = p_meanpT_vs_yCM_pr->ProjectionX();
  TH1D *h_meanpT_vs_yCM_pr_alt = p_meanpT_vs_yCM_pr_alt->ProjectionX();
  TH1D *h_meanpT_vs_yCM_de = p_meanpT_vs_yCM_de->ProjectionX();
  TH1D *h_meanpT_vs_yCM_tr = p_meanpT_vs_yCM_tr->ProjectionX();
  
  h_meanpT_vs_yCM_pp->SetTitle("");
  h_meanpT_vs_yCM_pm->SetTitle("");
  h_meanpT_vs_yCM_kp->SetTitle("");
  h_meanpT_vs_yCM_km->SetTitle("");
  h_meanpT_vs_yCM_pr->SetTitle("");
  h_meanpT_vs_yCM_pr_alt->SetTitle("");
  h_meanpT_vs_yCM_de->SetTitle("");
  h_meanpT_vs_yCM_tr->SetTitle("");

  h_meanpT_vs_yCM_pp = PlotUtils::trimRapidityPlot(h_meanpT_vs_yCM_pp);
  h_meanpT_vs_yCM_pm = PlotUtils::trimRapidityPlot(h_meanpT_vs_yCM_pm);
  h_meanpT_vs_yCM_kp = PlotUtils::trimRapidityPlot(h_meanpT_vs_yCM_kp);
  h_meanpT_vs_yCM_km = PlotUtils::trimRapidityPlot(h_meanpT_vs_yCM_km);
  h_meanpT_vs_yCM_pr = PlotUtils::trimRapidityPlot(h_meanpT_vs_yCM_pr);
  h_meanpT_vs_yCM_pr_alt = PlotUtils::trimRapidityPlot(h_meanpT_vs_yCM_pr_alt);
  h_meanpT_vs_yCM_de = PlotUtils::trimRapidityPlot(h_meanpT_vs_yCM_de);
  h_meanpT_vs_yCM_tr = PlotUtils::trimRapidityPlot(h_meanpT_vs_yCM_tr);

  h_meanpT_vs_yCM_pp->SetMarkerStyle(20);
  h_meanpT_vs_yCM_pm->SetMarkerStyle(20);
  h_meanpT_vs_yCM_kp->SetMarkerStyle(20);
  h_meanpT_vs_yCM_km->SetMarkerStyle(20);
  h_meanpT_vs_yCM_pr->SetMarkerStyle(20);
  h_meanpT_vs_yCM_pr_alt->SetMarkerStyle(20);
  h_meanpT_vs_yCM_de->SetMarkerStyle(20);
  h_meanpT_vs_yCM_tr->SetMarkerStyle(20);

  h_meanpT_vs_yCM_pp->SetMarkerSize(2.5);
  h_meanpT_vs_yCM_pm->SetMarkerSize(2.5);
  h_meanpT_vs_yCM_kp->SetMarkerSize(2.5);
  h_meanpT_vs_yCM_km->SetMarkerSize(2.5);
  h_meanpT_vs_yCM_pr->SetMarkerSize(2.5);
  h_meanpT_vs_yCM_pr_alt->SetMarkerSize(2.5);
  h_meanpT_vs_yCM_de->SetMarkerSize(2.5);
  h_meanpT_vs_yCM_tr->SetMarkerSize(2.5);

  h_meanpT_vs_yCM_pp->SetLineColor(1);
  h_meanpT_vs_yCM_pm->SetLineColor(1);
  h_meanpT_vs_yCM_kp->SetLineColor(1);
  h_meanpT_vs_yCM_km->SetLineColor(1);
  h_meanpT_vs_yCM_pr->SetLineColor(1);
  h_meanpT_vs_yCM_pr_alt->SetLineColor(1);
  h_meanpT_vs_yCM_de->SetLineColor(1);
  h_meanpT_vs_yCM_tr->SetLineColor(1);

  h_meanpT_vs_yCM_pp->SetLineWidth(3);
  h_meanpT_vs_yCM_pm->SetLineWidth(3);
  h_meanpT_vs_yCM_kp->SetLineWidth(3);
  h_meanpT_vs_yCM_km->SetLineWidth(3);
  h_meanpT_vs_yCM_pr->SetLineWidth(3);
  h_meanpT_vs_yCM_pr_alt->SetLineWidth(3);
  h_meanpT_vs_yCM_de->SetLineWidth(3);
  h_meanpT_vs_yCM_tr->SetLineWidth(3);

  h_meanpT_vs_yCM_pp->GetYaxis()->SetTitleOffset(1.3);
  h_meanpT_vs_yCM_pm->GetYaxis()->SetTitleOffset(1.3);
  h_meanpT_vs_yCM_kp->GetYaxis()->SetTitleOffset(1.3);
  h_meanpT_vs_yCM_km->GetYaxis()->SetTitleOffset(1.3);
  h_meanpT_vs_yCM_pr->GetYaxis()->SetTitleOffset(1.3);
  h_meanpT_vs_yCM_pr_alt->GetYaxis()->SetTitleOffset(1.3);
  h_meanpT_vs_yCM_de->GetYaxis()->SetTitleOffset(1.3);
  h_meanpT_vs_yCM_tr->GetYaxis()->SetTitleOffset(1.3);


  Double_t min = 0.0;
  Double_t max = 2.2;
  
  h_meanpT_vs_yCM_pp->SetMinimum(min);
  h_meanpT_vs_yCM_pm->SetMinimum(min);
  h_meanpT_vs_yCM_kp->SetMinimum(min);
  h_meanpT_vs_yCM_km->SetMinimum(min);
  h_meanpT_vs_yCM_pr->SetMinimum(min);
  h_meanpT_vs_yCM_pr_alt->SetMinimum(min);
  h_meanpT_vs_yCM_de->SetMinimum(min);
  h_meanpT_vs_yCM_tr->SetMinimum(min);

  h_meanpT_vs_yCM_pp->SetMaximum(max);
  h_meanpT_vs_yCM_pm->SetMaximum(max);
  h_meanpT_vs_yCM_kp->SetMaximum(max);
  h_meanpT_vs_yCM_km->SetMaximum(max);
  h_meanpT_vs_yCM_pr->SetMaximum(max);
  h_meanpT_vs_yCM_pr_alt->SetMaximum(max);
  h_meanpT_vs_yCM_de->SetMaximum(max);
  h_meanpT_vs_yCM_tr->SetMaximum(max);


  TPaveText *ppText = new TPaveText(0.3, 1.35, 0.4, 1.8, "NB");
  ppText->AddText("#pi^{+}");
  ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  ppText->AddText("0.18 < p_{T} < 1.6 GeV");
  ppText->SetFillColorAlpha(0,0);
  ppText->SetLineColorAlpha(0,0);
  ppText->SetTextSize(.03);

  TPaveText *pmText = new TPaveText(0.3, 1.35, 0.4, 1.8, "NB");
  pmText->AddText("#pi^{-}");
  pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  pmText->AddText("0.18 < p_{T} < 1.6 GeV");
  pmText->SetFillColorAlpha(0,0);
  pmText->SetLineColorAlpha(0,0);
  pmText->SetTextSize(.03);

  TPaveText *kpText = new TPaveText(0.3, 1.35, 0.4, 1.8, "NB");
  kpText->AddText("K^{+}");
  kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  kpText->AddText("0.18 < p_{T} < 1.6 GeV");
  kpText->SetFillColorAlpha(0,0);
  kpText->SetLineColorAlpha(0,0);
  kpText->SetTextSize(.03);

  TPaveText *kmText = new TPaveText(0.3, 1.35, 0.4, 1.8, "NB");
  kmText->AddText("K^{-}");
  kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  kmText->AddText("0.18 < p_{T} < 1.6 GeV");
  kmText->SetFillColorAlpha(0,0);
  kmText->SetLineColorAlpha(0,0);
  kmText->SetTextSize(.03);

  TPaveText *prText = new TPaveText(0.3, 1.35, 0.4, 1.8, "NB");
  prText->AddText("Proton");
  prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  prText->AddText("0.4 < p_{T} < 2.0 GeV");
  prText->SetFillColorAlpha(0,0);
  prText->SetLineColorAlpha(0,0);
  prText->SetTextSize(.03);

  TPaveText *prAltText = new TPaveText(0.3, 1.35, 0.4, 1.8, "NB");
  prAltText->AddText("Proton");
  prAltText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  prAltText->AddText("1.0 < p_{T} < 2.5 GeV");
  prAltText->SetFillColorAlpha(0,0);
  prAltText->SetLineColorAlpha(0,0);
  prAltText->SetTextSize(.03);

  TPaveText *deText = new TPaveText(0.3, 0.2, 0.4, 0.65, "NB");
  deText->AddText("Deuteron");
  deText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  deText->AddText("0.4 < p_{T} < 2.0 GeV");
  deText->SetFillColorAlpha(0,0);
  deText->SetLineColorAlpha(0,0);
  deText->SetTextSize(.03);
  
  TPaveText *trText = new TPaveText(0.3, 0.2, 0.4, 0.65, "NB");
  trText->AddText("Triton");
  trText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  trText->AddText("0.4 < p_{T} < 2.0 GeV");
  trText->SetFillColorAlpha(0,0);
  trText->SetLineColorAlpha(0,0);
  trText->SetTextSize(.03);

  h_meanpT_vs_yCM_pp->Draw();
  ppText->Draw();
  canvas->SaveAs(jobID + "_meanpT_pp.png");
  canvas->Clear();

  h_meanpT_vs_yCM_pm->Draw();
  pmText->Draw();
  canvas->SaveAs(jobID + "_meanpT_pm.png");
  canvas->Clear();

  h_meanpT_vs_yCM_kp->Draw();
  kpText->Draw();
  canvas->SaveAs(jobID + "_meanpT_kp.png");
  canvas->Clear();

  h_meanpT_vs_yCM_km->Draw();
  kmText->Draw();
  canvas->SaveAs(jobID + "_meanpT_km.png");
  canvas->Clear();

  h_meanpT_vs_yCM_pr->Draw();
  prText->Draw();
  canvas->SaveAs(jobID + "_meanpT_pr.png");
  canvas->Clear();

  h_meanpT_vs_yCM_pr_alt->Draw();
  prAltText->Draw();
  canvas->SaveAs(jobID + "_meanpT_pr_alt.png");
  canvas->Clear();

  h_meanpT_vs_yCM_de->Draw();
  deText->Draw();
  canvas->SaveAs(jobID + "_meanpT_de.png");
  canvas->Clear();

  h_meanpT_vs_yCM_tr->Draw();
  trText->Draw();
  canvas->SaveAs(jobID + "_meanpT_tr.png");
  canvas->Clear();

}
