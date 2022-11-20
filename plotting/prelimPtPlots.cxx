void prelimPtPlots(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {std::cout << "Wrong file!" << std::endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1000, 1000);
  //canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetLeftMargin(0.15);
  canvas->SetTopMargin(0.04);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.1);
  canvas->cd();

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetLineWidth(3);
  gStyle->SetOptDate();

  TProfile2D *p2_vn_pT_cent_pr = (TProfile2D*)file->Get("p2_vn_pT_cent_pr");
  TProfile *p_vn_pT_00to10_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_00to10_pr", 15, 16);
  TProfile *p_vn_pT_10to40_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_10to40_pr", 9, 14);
  TProfile *p_vn_pT_40to60_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_40to60_pr", 5, 8);

  TH1D *h_vn_pT_00to10_pr = new TH1D("h_vn_pT_00to10_pr", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pr = new TH1D("h_vn_pT_10to40_pr", ";p_{T} (GeV);v_{3}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pr = new TH1D("h_vn_pT_40to60_pr", ";p_{T} (GeV);v_{3}", 10, 0, 2);

  // Convert profiles to histograms
  h_vn_pT_00to10_pr = p_vn_pT_00to10_pr->ProjectionX();
  h_vn_pT_10to40_pr = p_vn_pT_10to40_pr->ProjectionX();
  h_vn_pT_40to60_pr = p_vn_pT_40to60_pr->ProjectionX();


  // Retrieve systematic uncertainties
  TFile* systematicFile = TFile::Open("systematicErrors.root", "READ");
  TGraphErrors* sys_pT_00to10_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pT_00to10_pr_px");
  TGraphErrors* sys_pT_10to40_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pT_10to40_pr_px");
  TGraphErrors* sys_pT_40to60_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pT_40to60_pr_px");
  ////


  // Set various aesthetics
  h_vn_pT_00to10_pr->SetMarkerStyle(21);
  h_vn_pT_10to40_pr->SetMarkerStyle(20);
  h_vn_pT_40to60_pr->SetMarkerStyle(22);
  h_vn_pT_00to10_pr->SetMarkerColor(kRed-4);
  h_vn_pT_10to40_pr->SetMarkerColor(kBlue-4);
  h_vn_pT_40to60_pr->SetMarkerColor(kGreen+1);
  h_vn_pT_00to10_pr->SetMarkerSize(2);
  h_vn_pT_10to40_pr->SetMarkerSize(2);
  h_vn_pT_40to60_pr->SetMarkerSize(2.5);
  h_vn_pT_00to10_pr->SetLineColor(kRed-4);
  h_vn_pT_10to40_pr->SetLineColor(kBlue-4);
  h_vn_pT_40to60_pr->SetLineColor(kGreen+1);
  h_vn_pT_00to10_pr->SetLineWidth(3);
  h_vn_pT_10to40_pr->SetLineWidth(3);
  h_vn_pT_40to60_pr->SetLineWidth(3);

  sys_pT_00to10_pr->SetMarkerColor(kRed-4);
  sys_pT_10to40_pr->SetMarkerColor(kBlue-4);
  sys_pT_40to60_pr->SetMarkerColor(kGreen+1);
  sys_pT_00to10_pr->SetLineColor(kRed-4);
  sys_pT_10to40_pr->SetLineColor(kBlue-4);
  sys_pT_40to60_pr->SetLineColor(kGreen+1);
  ////
  
  THStack *prPtStack   = new THStack("prPtStack", ";p_{T} (GeV/c);v_{3} {#Psi_{1}}");
  prPtStack->Add(h_vn_pT_00to10_pr);
  prPtStack->Add(h_vn_pT_10to40_pr);
  prPtStack->Add(h_vn_pT_40to60_pr);


  // Make text boxes, legends, and line at zero
  TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
  prLegend->AddEntry(h_vn_pT_00to10_pr, "0 - 10%");
  prLegend->AddEntry(h_vn_pT_10to40_pr, "10 - 40%");
  prLegend->AddEntry(h_vn_pT_40to60_pr, "40 - 60%");
  prLegend->SetBorderSize(0);
  prLegend->SetFillColorAlpha(0,0);
  prLegend->SetTextSize(0.04);

  //TPaveText *prText = new TPaveText(0.28, 0.03, 1.28, 0.055, "NB");
  TPaveText *prText = new TPaveText(0.0, 0.02, 1.8, 0.055, "NB");
  prText->SetTextAlign(12);
  prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT (year 2018)");
  prText->AddText("Proton");
  prText->AddText("0 < y_{CM} < 0.5");
  prText->SetFillColorAlpha(0,0);
  prText->SetLineColorAlpha(0,0);
  prText->SetTextSize(0.045);

  //TPaveText* prelimText = new TPaveText(0.3, 0.02, 1.28, 0.03, "NB");
  TPaveText* prelimText = new TPaveText(0.7, -0.08, 1.4, -0.07, "NB");
  prelimText->AddText("STAR Preliminary");
  prelimText->SetTextColor(kRed);
  prelimText->SetFillColorAlpha(0,0);
  prelimText->SetTextSize(0.04);

  TLine *zeroLine_pt = new TLine(0, 0, 2, 0);
  zeroLine_pt->SetLineStyle(9);
  ////

  prPtStack->Draw();
  prPtStack->GetYaxis()->SetLabelSize(0.043);
  prPtStack->GetXaxis()->SetLabelSize(0.043);
  prPtStack->GetYaxis()->SetTitleOffset(1.4);
  prPtStack->GetXaxis()->SetTitleOffset(1.0);
  prPtStack->GetXaxis()->SetNdivisions(210);
  prPtStack->GetXaxis()->SetTitleSize(0.045);
  prPtStack->GetYaxis()->SetTitleSize(0.05);
  prPtStack->Draw();
  prPtStack->SetMaximum(0.06);
  prPtStack->SetMinimum(-0.12);
  prPtStack->Draw("NOSTACK E1P");
  zeroLine_pt->Draw("SAME");
  prPtStack->Draw("NOSTACK E1P SAME");
  sys_pT_00to10_pr->Draw("[]");
  sys_pT_10to40_pr->Draw("[]");
  sys_pT_40to60_pr->Draw("[]");
  prLegend->Draw();
  prText->Draw();
  //prelimText->Draw();
  canvas->SaveAs("v3_prPtStack.pdf");
  canvas->Clear();
  
  systematicFile->Close();
  file->Close();
}
