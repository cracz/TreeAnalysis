#include "PlotUtils.h"

void prelimCentralityPlots(TString jobID, Bool_t useSystematics = false, Double_t sqrt_s_NN = 3.0)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {std::cout << "Wrong file!" << std::endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1000, 1000);//1200, 1200);
  canvas->SetGridx(0);
  canvas->SetGridy(0);
  canvas->SetLogy(0);
  canvas->SetTicks();
  canvas->SetTopMargin(0.04);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetLeftMargin(0.17);
  canvas->cd();
  
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetLineWidth(3);
  gStyle->SetOptDate(0);


  TProfile *p_vn_pp = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_kp = (TProfile*)file->Get("p_vn_kp");
  TProfile *p_vn_km = (TProfile*)file->Get("p_vn_km");
  TProfile *p_vn_pr = (TProfile*)file->Get("p_vn_pr"); 
  p_vn_kp->Rebin();
  p_vn_km->Rebin();

  // Convert profiles to histograms
  TH1D *h_vn_pp = p_vn_pp->ProjectionX();
  TH1D *h_vn_pm = p_vn_pm->ProjectionX();
  TH1D *h_vn_kp = p_vn_kp->ProjectionX();
  TH1D *h_vn_km = p_vn_km->ProjectionX();
  TH1D *h_vn_pr = p_vn_pr->ProjectionX();

  // Flip centrality plots
  h_vn_pp = PlotUtils::flipHisto(h_vn_pp);
  h_vn_pm = PlotUtils::flipHisto(h_vn_pm);
  h_vn_kp = PlotUtils::flipHisto(h_vn_kp);
  h_vn_km = PlotUtils::flipHisto(h_vn_km);
  h_vn_pr = PlotUtils::flipHisto(h_vn_pr);

  // Trim and clean up x-axis
  h_vn_pp = PlotUtils::trimCentralityPlot(h_vn_pp);
  h_vn_pm = PlotUtils::trimCentralityPlot(h_vn_pm);
  h_vn_kp = PlotUtils::trimCentralityPlot(h_vn_kp);
  h_vn_km = PlotUtils::trimCentralityPlot(h_vn_km);
  h_vn_pr = PlotUtils::trimCentralityPlot(h_vn_pr);

  // Retrieve systematic uncertainties
  TFile* systematicFile;
  TGraphErrors* sys_pp;
  TGraphErrors* sys_pm;
  TGraphErrors* sys_kp;
  TGraphErrors* sys_km;
  TGraphErrors* sys_pr;

  if (useSystematics)
    {
      if (sqrt_s_NN == 3.0)
	systematicFile = TFile::Open("systematicErrors_3p0GeV.root", "READ");
      else if (sqrt_s_NN == 3.2)
	systematicFile = TFile::Open("systematicErrors_3p2GeV_SL23d.root", "READ");
      else if (sqrt_s_NN == 3.5)
	systematicFile = TFile::Open("systematicErrors_3p5GeV_SL23d.root", "READ");
      else if (sqrt_s_NN == 3.9)
	systematicFile = TFile::Open("systematicErrors_3p9GeV_SL23d.root", "READ");
      else if (sqrt_s_NN == 4.5)
	systematicFile = TFile::Open("systematicErrors_4p5GeV_SL23d.root", "READ");
      jobID = "Normal";
      sys_pp = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pp_"+jobID+"_flip");
      sys_pm = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pm_"+jobID+"_flip");
      sys_kp = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_kp_"+jobID+"_flip");
      sys_km = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_km_"+jobID+"_flip");
      sys_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pr_"+jobID+"_flip");

      sys_pp->SetMarkerColor(1);
      sys_pp->SetLineColor(1);

      sys_pm->SetMarkerColor(4);
      sys_pm->SetLineColor(4);

      sys_kp->SetMarkerColor(2);
      sys_kp->SetLineColor(2);

      sys_km->SetMarkerColor(4);
      sys_km->SetLineColor(4);

      //sys_pr should be fine as it is
    }
  ////
  
  // Set various aesthetics
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

  h_vn_kp->SetMarkerStyle(20);
  h_vn_kp->SetMarkerSize(2.5);
  h_vn_kp->SetMarkerColor(2);
  h_vn_kp->SetLineColor(2);
  h_vn_kp->SetLineWidth(3);
  h_vn_kp->GetYaxis()->SetTitleOffset(1.7);

  h_vn_km->SetMarkerStyle(24);
  h_vn_km->SetMarkerSize(3);
  h_vn_km->SetMarkerColor(4);
  h_vn_km->SetLineColor(4);
  h_vn_km->SetLineWidth(3);
  h_vn_km->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pr->SetMarkerStyle(20);
  h_vn_pr->SetMarkerSize(2.5);
  h_vn_pr->SetMarkerColor(kRed-4);
  h_vn_pr->SetLineColor(kRed-4);
  h_vn_pr->SetLineWidth(3);
  h_vn_pr->GetYaxis()->SetTitleOffset(1.7);
  ////

  // Remove points that have unreliable resolutions
  if (sqrt_s_NN == 3.5 || sqrt_s_NN == 3.9)
    {
      for (Int_t i = 9; i <= 12; i++)
	{
	  h_vn_pp->SetBinContent(i, 0.0);
	  h_vn_pp->SetBinError(i, 0.0);
	  h_vn_pm->SetBinContent(i, 0.0);
	  h_vn_pm->SetBinError(i, 0.0);
	  h_vn_pr->SetBinContent(i, 0.0);
	  h_vn_pr->SetBinError(i, 0.0);
	}
      for (Int_t i = 5; i <= 6; i++)
	{
	  h_vn_kp->SetBinContent(i, 0.0);
	  h_vn_kp->SetBinError(i, 0.0);
	  h_vn_km->SetBinContent(i, 0.0);
	  h_vn_km->SetBinError(i, 0.0);
	}
    }

  THStack *allCentralityStack = new THStack("allCentralityStack", ";Centrality (%);v_{3} {#Psi_{1}}");
  allCentralityStack->Add(h_vn_pp);
  allCentralityStack->Add(h_vn_pm);
  allCentralityStack->Add(h_vn_pr);

  THStack *kaCentralityStack = new THStack("kaCentralityStack", ";Centrality (%);v_{3} {#Psi_{1}}");
  kaCentralityStack->Add(h_vn_km);
  //kaCentralityStack->Add(h_vn_kp);
  
  // Make text boxes, legends, and line at zero
  TPaveText* prelimText = new TPaveText(15, 0.01, 45, 0.014, "NB");
  prelimText->AddText("STAR Preliminary");
  prelimText->SetTextColor(kRed);
  prelimText->SetFillColorAlpha(0,0);
  prelimText->SetTextSize(0.04);

   
  TPaveText *allText = new TPaveText(15, 0.007, 45, 0.018, "NB");
  if (sqrt_s_NN == 3.0)
    allText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");// (year 2018)");
  else if (sqrt_s_NN == 3.2 || sqrt_s_NN == 3.22)
    allText->AddText("Au+Au #sqrt{s_{NN}} = 3.2 GeV FXT");// (year 2019)");
  else if (sqrt_s_NN == 3.5)
    allText->AddText("Au+Au #sqrt{s_{NN}} = 3.5 GeV FXT");// (year 2020)");
  else if (sqrt_s_NN == 3.9)
    allText->AddText("Au+Au #sqrt{s_{NN}} = 3.9 GeV FXT");// (years 2019, 2020)");
  else if (sqrt_s_NN == 4.5)
    allText->AddText("Au+Au #sqrt{s_{NN}} = 4.5 GeV FXT");// (year 2020)");
  allText->AddText("0 < y_{c.m.} < 0.5");
  allText->AddText("STAR");
  allText->SetFillColorAlpha(0,0);
  allText->SetLineColorAlpha(0,0);
  allText->SetTextFont(23);
  allText->SetTextSize(45);

  TPaveText *kaText = new TPaveText(15, 0.038, 45, 0.095, "NB");//(15, 0.038, 48, 0.095, "NB");
  if (sqrt_s_NN == 3.0)
    kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  else if (sqrt_s_NN == 3.2 || sqrt_s_NN == 3.22)
    kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.2 GeV FXT");
  else if (sqrt_s_NN == 3.5)
    kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.5 GeV FXT");
  else if (sqrt_s_NN == 3.9)
    kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.9 GeV FXT");
  else if (sqrt_s_NN == 4.5)
    kaText->AddText("Au+Au #sqrt{s_{NN}} = 4.5 GeV FXT");
  kaText->AddText("0 < y_{c.m.} < 0.5");
  kaText->AddText("0.4 < p_{T} < 1.6 GeV/c");
  kaText->AddText("STAR");
  kaText->SetFillColorAlpha(0,0);
  kaText->SetLineColorAlpha(0,0);
  kaText->SetTextFont(23);
  kaText->SetTextSize(45);

  TLegend *allLegend = new TLegend(0.25, 0.14, 0.6, 0.3);
  allLegend->AddEntry(h_vn_pp,"#pi^{+}, 0.18 #leq p_{T} #leq 1.6 GeV/c");
  allLegend->AddEntry(h_vn_pm,"#pi^{-}, 0.18 #leq p_{T} #leq 1.6 GeV/c");
  allLegend->AddEntry(h_vn_pr,"p, 0.4 #leq p_{T} #leq 2.0 GeV/c");
  allLegend->SetFillColorAlpha(0,0);
  allLegend->SetLineColorAlpha(0,0);
  allLegend->SetTextFont(23);
  allLegend->SetTextSize(45);

  TLegend *kaLegend = new TLegend(0.5, 0.14, 0.7, 0.25);
  kaLegend->AddEntry(h_vn_kp,"K^{+}");
  kaLegend->AddEntry(h_vn_km,"K^{-}");
  kaLegend->SetFillColorAlpha(0,0);
  kaLegend->SetLineColorAlpha(0,0);
  kaLegend->SetTextFont(23);
  kaLegend->SetTextSize(45);


  TLine *zeroLine = new TLine(0, 0, 60, 0);
  zeroLine->SetLineStyle(9);
  zeroLine->SetLineWidth(3);
  ////

  allCentralityStack->Draw();
  allCentralityStack->GetXaxis()->SetLabelFont(133);
  allCentralityStack->GetYaxis()->SetLabelFont(133);
  allCentralityStack->GetXaxis()->SetLabelSize(55);
  allCentralityStack->GetYaxis()->SetLabelSize(55);
  allCentralityStack->GetXaxis()->SetTitleOffset(1.0);
  allCentralityStack->GetYaxis()->SetTitleOffset(1.5);
  allCentralityStack->GetXaxis()->SetTitleFont(133);
  allCentralityStack->GetYaxis()->SetTitleFont(133);
  allCentralityStack->GetXaxis()->SetTitleSize(55);
  allCentralityStack->GetYaxis()->SetTitleSize(55);
  allCentralityStack->GetXaxis()->SetNdivisions(210);
  allCentralityStack->GetYaxis()->SetNdivisions(505);
  allCentralityStack->SetMaximum(0.02);
  allCentralityStack->SetMinimum(-0.03);
  allCentralityStack->Draw("NOSTACK X0 EP");
  zeroLine->Draw("SAME");
  allCentralityStack->Draw("NOSTACK EP X0 SAME");
  if (useSystematics)
    {
      for (int i = 0; i < sys_pp->GetN(); i++)
	{
	  sys_pp->SetPointError(i, 0.0, sys_pp->GetErrorY(i));
	  sys_pm->SetPointError(i, 0.0, sys_pm->GetErrorY(i));
	  sys_pr->SetPointError(i, 0.0, sys_pr->GetErrorY(i));
	}
      sys_pp->Draw("[]");
      sys_pm->Draw("[]");
      sys_pr->Draw("[]");
    }
  allLegend->Draw();
  allText->Draw();
  //prelimText->Draw();
  canvas->SaveAs("v3_allCentralityStack.pdf");
  canvas->SaveAs("v3_allCentralityStack.png");
  canvas->Clear();


  TGraphErrors* g_kp;
  if (useSystematics)
    {
      g_kp = new TGraphErrors(h_vn_kp);
      for (int i = 0; i < g_kp->GetN(); i++)
	{
	  g_kp->SetPointError(i, 0.0, g_kp->GetErrorY(i));
	  sys_kp->SetPointError(i, 0.0, sys_kp->GetErrorY(i));
	  sys_km->SetPointError(i, 0.0, sys_km->GetErrorY(i));
	}
      PlotUtils::shiftGraphX(g_kp, -1.5);
      PlotUtils::shiftGraphX(sys_kp, -1.5);
    }
  
  kaCentralityStack->Draw();
  kaCentralityStack->GetXaxis()->SetLabelFont(133);
  kaCentralityStack->GetYaxis()->SetLabelFont(133);
  kaCentralityStack->GetXaxis()->SetLabelSize(55);
  kaCentralityStack->GetYaxis()->SetLabelSize(55);
  kaCentralityStack->GetXaxis()->SetTitleOffset(1.0);
  kaCentralityStack->GetYaxis()->SetTitleOffset(1.5);
  kaCentralityStack->GetXaxis()->SetTitleFont(133);
  kaCentralityStack->GetYaxis()->SetTitleFont(133);
  kaCentralityStack->GetXaxis()->SetTitleSize(55);
  kaCentralityStack->GetYaxis()->SetTitleSize(55);
  kaCentralityStack->GetXaxis()->SetNdivisions(210);
  kaCentralityStack->SetMaximum(0.1);
  kaCentralityStack->SetMinimum(-0.1);
  kaCentralityStack->Draw("NOSTACK EP X0");
  zeroLine->Draw("SAME");
  kaCentralityStack->Draw("NOSTACK EP X0 SAME");
  if (useSystematics)
    {
      g_kp->Draw("PZ");
      sys_kp->Draw("[]");
      sys_km->Draw("[]");
    }
  kaLegend->Draw();
  kaText->Draw();
  //prelimText->Draw();
  canvas->SaveAs("v3_kaCentralityStack.pdf");
  canvas->SaveAs("v3_kaCentralityStack.png");
  canvas->Clear();

  if (useSystematics)
    systematicFile->Close();

  file->Close();
}
