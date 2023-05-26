void resolutionPlot()
{
  TFile* file = TFile::Open("resolutionsWithSystematics.root","READ");

  TH1D* h_resolutions = (TH1D*)file->Get("h_resolutionsWithSystematics");

  TGraphErrors* finalGraph = new TGraphErrors(h_resolutions);
  finalGraph->SetTitle("");
  finalGraph->GetYaxis()->SetTitle("R_{31}");
  finalGraph->GetYaxis()->SetTitleSize(0.055);
  finalGraph->GetYaxis()->SetTitleOffset(1.15);
  finalGraph->GetYaxis()->SetLabelSize(0.045);
  finalGraph->GetYaxis()->SetRangeUser(0.0, 0.25);
  finalGraph->GetXaxis()->SetTitle("Centrality (%)");
  finalGraph->GetXaxis()->SetRangeUser(0.0, 60.0);
  finalGraph->GetXaxis()->SetTitleSize(0.045);
  finalGraph->GetXaxis()->SetTitleOffset(1.1);
  finalGraph->GetXaxis()->SetLabelSize(0.05);
  finalGraph->SetLineWidth(2);
  finalGraph->SetLineColor(kBlack);
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

  finalGraph->Draw();
  finalGraph->Draw("AP[]");
  canvas->SaveAs("resolutionPlot.pdf");

  delete finalGraph;
  delete canvas;
}
