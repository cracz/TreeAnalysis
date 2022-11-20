#include "PlotUtils.h"

void prelimPdtPlots(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {std::cout << "Wrong file!" << std::endl; return;}

  TFile *dAndt_JAM_file = TFile::Open("dAndt_JAM.root");
  if(!dAndt_JAM_file) {std::cout << "No JAM file!" << std::endl; return;}

  TProfile *p_vn_pr_alt_JAM = (TProfile*)dAndt_JAM_file->Get("p_vn_pr_alt");
  TProfile *p_vn_de_JAM = (TProfile*)dAndt_JAM_file->Get("p_vn_de");
  TProfile *p_vn_tr_JAM = (TProfile*)dAndt_JAM_file->Get("p_vn_tr");
  p_vn_pr_alt_JAM->SetName("p_vn_pr_alt_JAM");
  p_vn_de_JAM->SetName("p_vn_de_JAM");
  p_vn_tr_JAM->SetName("p_vn_tr_JAM");

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
  gStyle->SetOptDate();


  TProfile *p_vn_de   = (TProfile*)file->Get("p_vn_de");
  TProfile *p_vn_tr   = (TProfile*)file->Get("p_vn_tr");
  TProfile *p_vn_pr_alt   = (TProfile*)file->Get("p_vn_pr_alt");

  // Convert profiles to histograms
  TH1D *h_vn_de = p_vn_de->ProjectionX();
  TH1D *h_vn_tr = p_vn_tr->ProjectionX();
  TH1D *h_vn_pr_alt = p_vn_pr_alt->ProjectionX();
  TH1D *h_vn_pr_alt_JAM = p_vn_pr_alt_JAM->ProjectionX();
  TH1D *h_vn_de_JAM = p_vn_de_JAM->ProjectionX();
  TH1D *h_vn_tr_JAM = p_vn_tr_JAM->ProjectionX();

  // Flip centrality plots
  h_vn_de = PlotUtils::flipHisto(h_vn_de);
  h_vn_tr = PlotUtils::flipHisto(h_vn_tr);
  h_vn_pr_alt = PlotUtils::flipHisto(h_vn_pr_alt);
  h_vn_pr_alt_JAM = PlotUtils::flipHisto(h_vn_pr_alt_JAM);
  h_vn_de_JAM = PlotUtils::flipHisto(h_vn_de_JAM);
  h_vn_tr_JAM = PlotUtils::flipHisto(h_vn_tr_JAM);

  // Trim and clean up x-axis
  h_vn_de = PlotUtils::trimCentralityPlot(h_vn_de);
  h_vn_tr = PlotUtils::trimCentralityPlot(h_vn_tr);
  h_vn_pr_alt = PlotUtils::trimCentralityPlot(h_vn_pr_alt);
  h_vn_pr_alt_JAM = PlotUtils::trimCentralityPlot(h_vn_pr_alt_JAM);
  h_vn_de_JAM = PlotUtils::trimCentralityPlot(h_vn_de_JAM);
  h_vn_tr_JAM = PlotUtils::trimCentralityPlot(h_vn_tr_JAM);
  
  // Retrieve systematic uncertainties
  TFile* systematicFile = TFile::Open("systematicErrors.root", "READ");
  //jobID = "Normal_afterDuplication_piKefficiencies";
  TGraphErrors* sys_de = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_de_"+jobID+"_flip");
  TGraphErrors* sys_tr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_tr_"+jobID+"_flip");
  TGraphErrors* sys_pr_alt = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pr_alt_"+jobID+"_flip");
  ////

  // Scale the errors
  //sys_de->Scale(1.0/2.0, "y"); // Function doesn't exist for unknown reason??
  //sys_tr->Scale(1.0/3.0, "y");

  for (int i=0; i < sys_de->GetN(); i++)
    {
      Double_t yValue = sys_de->GetPointY(i);
      Double_t yError = sys_de->GetErrorY(i);
      
      sys_de->SetPointY(i, yValue/2.0);
      sys_de->SetPointError(i, 0.0, yError/2.0);
    }

  for (int i=0; i < sys_tr->GetN(); i++)
    {
      Double_t yValue = sys_tr->GetPointY(i);
      Double_t yError = sys_tr->GetErrorY(i);
      
      sys_tr->SetPointY(i, yValue/3.0);
      sys_tr->SetPointError(i, 0.0, yError/3.0);
    }
  ////

  // Set various aesthetics
  sys_de->SetMarkerColor(1);
  sys_de->SetLineColor(1);

  sys_tr->SetMarkerColor(4);
  sys_tr->SetLineColor(4);

  //sys_pr_alt should be fine as it is
  
  h_vn_de->SetMarkerStyle(21);
  h_vn_de->SetMarkerSize(2.5);
  h_vn_de->SetMarkerColor(1);
  h_vn_de->SetLineColor(1);
  h_vn_de->SetLineWidth(3);
  h_vn_de->GetYaxis()->SetTitleOffset(1.7);

  h_vn_tr->SetMarkerStyle(22);
  h_vn_tr->SetMarkerSize(3);
  h_vn_tr->SetMarkerColor(4);
  h_vn_tr->SetLineColor(4);
  h_vn_tr->SetLineWidth(3);
  h_vn_tr->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pr_alt->SetMarkerStyle(20);
  h_vn_pr_alt->SetMarkerSize(2.5);
  h_vn_pr_alt->SetMarkerColor(kRed-4);
  h_vn_pr_alt->SetLineColor(kRed-4);
  h_vn_pr_alt->SetLineWidth(3);
  h_vn_pr_alt->GetYaxis()->SetTitleOffset(1.7);

  //h_vn_pr_alt->SetMarkerStyle(20);
  h_vn_pr_alt->SetMarkerSize(2.5);
  h_vn_pr_alt->SetMarkerColor(kRed-4);
  h_vn_pr_alt_JAM->SetFillColor(kRed-4);
  h_vn_pr_alt_JAM->SetFillStyle(3244);
  h_vn_pr_alt->SetLineColor(kRed-4);
  h_vn_pr_alt->SetLineWidth(3);
  h_vn_pr_alt->GetYaxis()->SetTitleOffset(1.7);

  //h_vn_de_JAM->SetMarkerStyle(21);
  h_vn_de_JAM->SetMarkerSize(2.5);
  h_vn_de_JAM->SetMarkerColor(1);
  h_vn_de_JAM->SetFillColor(1);
  h_vn_de_JAM->SetFillStyle(3244);
  h_vn_de_JAM->SetLineColor(1);
  h_vn_de_JAM->SetLineWidth(3);
  h_vn_de_JAM->GetYaxis()->SetTitleOffset(1.7);

  //h_vn_tr_JAM->SetMarkerStyle(22);
  h_vn_tr_JAM->SetMarkerSize(3);
  h_vn_tr_JAM->SetMarkerColor(4);
  h_vn_tr_JAM->SetFillColor(4);
  h_vn_tr_JAM->SetFillStyle(3244);
  h_vn_tr_JAM->SetLineColor(4);
  h_vn_tr_JAM->SetLineWidth(3);
  h_vn_tr_JAM->GetYaxis()->SetTitleOffset(1.7);
  ////


  // Scale flow from deuteron and triton.
  TH1D *h_vn_de_scaled = (TH1D*)h_vn_de->Clone();
  h_vn_de_scaled->Scale(1.0/2.0);

  TH1D *h_vn_tr_scaled = (TH1D*)h_vn_tr->Clone();
  h_vn_tr_scaled->Scale(1.0/3.0);

  TH1D *h_vn_de_JAM_scaled = (TH1D*)h_vn_de_JAM->Clone();
  h_vn_de_JAM_scaled->Scale(1.0/2.0);

  TH1D *h_vn_tr_JAM_scaled = (TH1D*)h_vn_tr_JAM->Clone();
  h_vn_tr_JAM_scaled->Scale(1.0/3.0);

  
  THStack *pdtCentralityStack = new THStack("pdtCentralityStack", ";Centrality (%);v_{3} {#Psi_{1}} / A");
  pdtCentralityStack->Add(h_vn_de_scaled);
  pdtCentralityStack->Add(h_vn_tr_scaled);
  pdtCentralityStack->Add(h_vn_pr_alt);
  
  // Make text boxes, legends, and line at zero
  TPaveText* prelimText = new TPaveText(15, -0.05, 45, -0.04, "NB");
  prelimText->AddText("STAR Preliminary");
  prelimText->SetTextColor(kRed);
  prelimText->SetFillColorAlpha(0,0);
  prelimText->SetTextSize(0.04);

  TPaveText *pdtText = new TPaveText(15, 0.01, 45, 0.028, "NB");
  pdtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT (year 2018)");
  pdtText->AddText("0 < y_{CM} < 1.0");
  pdtText->AddText("0.04 #leq (m_{T}-m_{0})/A #leq 0.4 GeV");
  pdtText->SetFillColorAlpha(0,0);
  pdtText->SetLineColorAlpha(0,0);
  pdtText->SetTextSize(0.045);

  TLegend *pdtLegend = new TLegend(0.2, 0.14, 0.45, 0.35);
  pdtLegend->AddEntry(h_vn_pr_alt,"p");
  pdtLegend->AddEntry(h_vn_de_scaled,"d");
  pdtLegend->AddEntry(h_vn_tr_scaled,"t");
  pdtLegend->SetFillColorAlpha(0,0);
  pdtLegend->SetLineColorAlpha(0,0);
  pdtLegend->SetTextSize(0.04);

  TLine *zeroLine = new TLine(0, 0, 60, 0);
  zeroLine->SetLineStyle(9);
  zeroLine->SetLineWidth(3);
  ////

  h_vn_pr_alt_JAM->GetYaxis()->SetRangeUser(-0.06, 0.03);
  h_vn_de_JAM_scaled->GetYaxis()->SetRangeUser(-0.06, 0.03);
  h_vn_tr_JAM_scaled->GetYaxis()->SetRangeUser(-0.06, 0.03);
  
  pdtCentralityStack->Draw();
  pdtCentralityStack->GetXaxis()->SetLabelSize(0.045);
  pdtCentralityStack->GetYaxis()->SetLabelSize(0.045);
  pdtCentralityStack->GetXaxis()->SetTitleOffset(1.0);
  pdtCentralityStack->GetYaxis()->SetTitleOffset(1.4);
  pdtCentralityStack->GetXaxis()->SetTitleSize(0.045);
  pdtCentralityStack->GetYaxis()->SetTitleSize(0.05);
  pdtCentralityStack->GetXaxis()->SetNdivisions(210);
  pdtCentralityStack->SetMaximum(0.03);
  pdtCentralityStack->SetMinimum(-0.06);
  pdtCentralityStack->Draw("NOSTACK E1P");
  zeroLine->Draw("SAME");
  //h_vn_pr_alt_JAM->Draw("E3 SAME");
  //h_vn_de_JAM_scaled->Draw("E3 SAME");
  //h_vn_tr_JAM_scaled->Draw("E3 SAME");
  pdtCentralityStack->Draw("NOSTACK E1P SAME");
  sys_de->Draw("[]");
  sys_tr->Draw("[]");
  sys_pr_alt->Draw("[]");
  pdtLegend->Draw();
  pdtText->Draw();
  
  //prelimText->Draw();
  canvas->SaveAs("v3_pdtCentralityStack.pdf");
  canvas->Clear();

  systematicFile->Close();
  file->Close();
}
