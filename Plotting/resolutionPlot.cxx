#include "PlotUtils.h"

void resolutionPlot()
{
  TFile* file = TFile::Open("eventPlaneSystematics_3p9GeV_1to5_max10_SL23d.root","READ");

  TH1D* h_resolutionsWithStats = (TH1D*)file->Get("h_resolutionsWithStats_flip");
  TH1D* h_resolutionsWithSysts = (TH1D*)file->Get("h_resolutionsWithSysts_flip");
  //h_resolutionsWithSysts = PlotUtils::trimCentralityPlot(h_resolutionsWithSysts);
  
  h_resolutionsWithStats->SetTitle("");
  h_resolutionsWithStats->GetYaxis()->SetTitle("R_{31}");
  h_resolutionsWithStats->GetYaxis()->SetTitleSize(0.055);
  h_resolutionsWithStats->GetYaxis()->SetTitleOffset(1.0);
  h_resolutionsWithStats->GetYaxis()->SetLabelSize(0.045);
  h_resolutionsWithStats->GetYaxis()->SetRangeUser(0.0, 0.25);// 1.0);
  h_resolutionsWithStats->GetXaxis()->SetTitle("Centrality (%)");
  h_resolutionsWithStats->GetXaxis()->SetTitleFont(132);
  h_resolutionsWithStats->GetYaxis()->SetTitleFont(132);

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
  h_resolutionsWithStats->SetFillColorAlpha(kBlue-4, 0.3);

  TGraphErrors* finalGraph = new TGraphErrors(h_resolutionsWithSysts);
  finalGraph->SetTitle("");
  finalGraph->GetYaxis()->SetTitle("R_{31}");
  finalGraph->GetYaxis()->SetTitleSize(0.055);
  finalGraph->GetYaxis()->SetTitleOffset(1.0);
  finalGraph->GetYaxis()->SetLabelSize(0.045);
  finalGraph->GetYaxis()->SetRangeUser(0.0, 0.25);// 1.0);
  finalGraph->GetXaxis()->SetTitle("Centrality (%)");
  finalGraph->GetXaxis()->SetTitleFont(132);
  finalGraph->GetYaxis()->SetTitleFont(132);
  
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
  finalGraph->SetFillColorAlpha(kBlue-4, 0.3);

  // Remove x errors
  //for (int i = 0; i < finalGraph->GetN(); i++)
  //finalGraph->SetPointError(i, 0.0, finalGraph->GetErrorY(i));



  /*
  ////// v1 results
  TFile* file_v1 = TFile::Open("eventPlaneSystematics_3p0GeV_v1.root","READ");
  TH1D* h_resolutionsWithStats_v1 = (TH1D*)file_v1->Get("h_resolutionsWithStats_flip");
  TH1D* h_resolutionsWithSysts_v1 = (TH1D*)file_v1->Get("h_resolutionsWithSysts_flip");
  //h_resolutionsWithSysts_v1 = PlotUtils::trimCentralityPlot(h_resolutionsWithSysts_v1);
  
  h_resolutionsWithStats_v1->SetTitle("");
  h_resolutionsWithStats_v1->GetYaxis()->SetTitle("");
  h_resolutionsWithStats_v1->GetYaxis()->SetTitleSize(0.055);
  h_resolutionsWithStats_v1->GetYaxis()->SetTitleOffset(1.0);
  h_resolutionsWithStats_v1->GetYaxis()->SetLabelSize(0.045);
  h_resolutionsWithStats_v1->GetYaxis()->SetRangeUser(0.0, 1.0);// 0.3);
  h_resolutionsWithStats_v1->GetXaxis()->SetTitle("Centrality (%)");
  h_resolutionsWithStats_v1->GetXaxis()->SetTitleFont(132);
  h_resolutionsWithStats_v1->GetYaxis()->SetTitleFont(132);

  //h_resolutionsWithStats_v1->GetXaxis()->SetRangeUser(0.0, 40.0);
  h_resolutionsWithStats_v1->GetXaxis()->SetRangeUser(0.0, 60.0);
  
  h_resolutionsWithStats_v1->GetXaxis()->SetTitleSize(0.045);
  h_resolutionsWithStats_v1->GetXaxis()->SetTitleOffset(1.1);
  h_resolutionsWithStats_v1->GetXaxis()->SetLabelSize(0.05);
  h_resolutionsWithStats_v1->SetLineWidth(2);
  h_resolutionsWithStats_v1->SetLineColor(kBlack);
  h_resolutionsWithStats_v1->SetMarkerStyle(21);
  h_resolutionsWithStats_v1->SetMarkerSize(2);
  h_resolutionsWithStats_v1->SetMarkerColor(kBlack);
  h_resolutionsWithStats_v1->SetFillColorAlpha(kBlack, 0.3);

  TGraphErrors* finalGraph_v1 = new TGraphErrors(h_resolutionsWithSysts_v1);
  finalGraph_v1->SetTitle("");
  finalGraph_v1->GetYaxis()->SetTitle("");
  finalGraph_v1->GetYaxis()->SetTitleSize(0.055);
  finalGraph_v1->GetYaxis()->SetTitleOffset(1.0);
  finalGraph_v1->GetYaxis()->SetLabelSize(0.045);
  finalGraph_v1->GetYaxis()->SetRangeUser(0.0, 1.0);// 0.3);
  finalGraph_v1->GetXaxis()->SetTitle("Centrality (%)");
  finalGraph_v1->GetXaxis()->SetTitleFont(132);
  finalGraph_v1->GetYaxis()->SetTitleFont(132);
  
  //finalGraph_v1->GetXaxis()->SetRangeUser(0.0, 40.0);
  finalGraph_v1->GetXaxis()->SetRangeUser(0.0, 60.0);

  finalGraph_v1->GetXaxis()->SetTitleSize(0.045);
  finalGraph_v1->GetXaxis()->SetTitleOffset(1.1);
  finalGraph_v1->GetXaxis()->SetLabelSize(0.05);
  finalGraph_v1->SetLineWidth(2);
  finalGraph_v1->SetLineColor(kBlack);
  finalGraph_v1->SetMarkerStyle(21);
  finalGraph_v1->SetMarkerSize(2);
  finalGraph_v1->SetMarkerColor(kBlack);
  finalGraph_v1->SetFillColorAlpha(kBlack, 0.3);

  // Remove x errors
  //for (int i = 0; i < finalGraph_v1->GetN(); i++)
  //finalGraph_v1->SetPointError(i, 0.0, finalGraph_v1->GetErrorY(i));
  //////////////
  */


  
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 1000);
  //canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetTopMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.13);
  canvas->cd();

  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetErrorX(0);

  TPaveText *text = new TPaveText(2, 0.85, 30, 0.95, "NB");
  //TPaveText *text = new TPaveText(2, 0.265, 30, 0.29, "NB");
  text->AddText("STAR Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  //text->AddText("E_{beam} = 3.85 GeV");
  //text->AddText("STAR");
  text->SetFillColorAlpha(0,0);
  text->SetLineColorAlpha(0,0);
  text->SetTextSize(0.045);
  text->SetTextAlign(13);
  text->SetTextFont(22);

  TLegend *legend = new TLegend(0.77, 0.39, 0.98, 0.55);
  //legend->AddEntry(h_resolutionsWithStats_v1, "n = 1");
  legend->AddEntry(h_resolutionsWithStats, "n = 3");
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0,0);
  legend->SetTextSize(0.04);
  legend->SetTextFont(132);
  

  h_resolutionsWithStats->Draw("EP");
  canvas->Update();
  finalGraph->Draw("2");
  /*
  h_resolutionsWithStats_v1->Draw("EP SAME");
  canvas->Update();
  finalGraph_v1->Draw("2");
  */
  //text->Draw();
  //legend->Draw();
  canvas->SaveAs("resolutionPlot.pdf");
  canvas->SaveAs("resolutionPlot.png");

  delete finalGraph;
  delete canvas;
}
