void prelimPtPlots(TString jobID, Bool_t useSystematics = false, Double_t sqrt_s_NN = 3.0)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {std::cout << "Wrong file!" << std::endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1000, 1000);
  canvas->SetGridx(0);
  canvas->SetGridy(0);
  canvas->SetTicks();
  canvas->SetLeftMargin(0.17);
  canvas->SetTopMargin(0.04);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->cd();

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetLineWidth(3);
  gStyle->SetOptDate(0);

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
  TFile* systematicFile;
  TGraphErrors* sys_pT_00to10_pr;
  TGraphErrors* sys_pT_10to40_pr;
  TGraphErrors* sys_pT_40to60_pr;

  if (useSystematics)
    {
      systematicFile = TFile::Open("systematicErrors_3p0GeV.root", "READ");
      sys_pT_00to10_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pT_00to10_pr_px");
      sys_pT_10to40_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pT_10to40_pr_px");
      sys_pT_40to60_pr = (TGraphErrors*)systematicFile->Get("Graph_from_p_vn_pT_40to60_pr_px");
    }
  ////


  // Set various aesthetics
  h_vn_pT_00to10_pr->SetMarkerStyle(21);
  h_vn_pT_10to40_pr->SetMarkerStyle(20);
  h_vn_pT_40to60_pr->SetMarkerStyle(33);
  h_vn_pT_00to10_pr->SetMarkerColor(kRed-4);
  h_vn_pT_10to40_pr->SetMarkerColor(kBlack);
  h_vn_pT_40to60_pr->SetMarkerColor(kBlue-4);
  h_vn_pT_00to10_pr->SetMarkerSize(2);
  h_vn_pT_10to40_pr->SetMarkerSize(2);
  h_vn_pT_40to60_pr->SetMarkerSize(3);
  h_vn_pT_00to10_pr->SetLineColor(kRed-4);
  h_vn_pT_10to40_pr->SetLineColor(kBlack);
  h_vn_pT_40to60_pr->SetLineColor(kBlue-4);
  h_vn_pT_00to10_pr->SetLineWidth(3);
  h_vn_pT_10to40_pr->SetLineWidth(3);
  h_vn_pT_40to60_pr->SetLineWidth(3);

  if (useSystematics)
    {
      sys_pT_00to10_pr->SetMarkerColor(kRed-4);
      sys_pT_10to40_pr->SetMarkerColor(kBlack);
      sys_pT_40to60_pr->SetMarkerColor(kBlue-4);
      sys_pT_00to10_pr->SetLineColor(kRed-4);
      sys_pT_10to40_pr->SetLineColor(kBlack);
      sys_pT_40to60_pr->SetLineColor(kBlue-4);
    }
  ////

  THStack *prPtStack   = new THStack("prPtStack", ";p_{T} (GeV/c);v_{3} {#Psi_{1}}");
  prPtStack->Add(h_vn_pT_00to10_pr);
  prPtStack->Add(h_vn_pT_10to40_pr);
  if (sqrt_s_NN < 3.5)
    prPtStack->Add(h_vn_pT_40to60_pr);


  // Make text boxes, legends, and line at zero
  TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
  prLegend->AddEntry(h_vn_pT_00to10_pr, "0 - 10%");
  prLegend->AddEntry(h_vn_pT_10to40_pr, "10 - 40%");
  if (sqrt_s_NN < 3.5)
    prLegend->AddEntry(h_vn_pT_40to60_pr, "40 - 60%");
  prLegend->SetBorderSize(0);
  prLegend->SetFillColorAlpha(0,0);
  prLegend->SetTextFont(23);
  prLegend->SetTextSize(45);

  //TPaveText *prText = new TPaveText(0.28, 0.03, 1.28, 0.055, "NB");
  TPaveText *prText = new TPaveText(0.0, 0.003, 1.8, 0.035, "NB");
  prText->SetTextAlign(12);
  if (sqrt_s_NN == 3.0)
    prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");// (year 2018)");
  else if (sqrt_s_NN == 3.2 || sqrt_s_NN == 3.22)
    prText->AddText("Au+Au #sqrt{s_{NN}} = 3.2 GeV FXT");// (year 2019)");
  else if (sqrt_s_NN == 3.5)
    prText->AddText("Au+Au #sqrt{s_{NN}} = 3.5 GeV FXT");// (year 2020)");
  else if (sqrt_s_NN == 3.9)
    prText->AddText("Au+Au #sqrt{s_{NN}} = 3.9 GeV FXT");// (years 2019, 2020)");
  else if (sqrt_s_NN == 4.5)
    prText->AddText("Au+Au #sqrt{s_{NN}} = 4.5 GeV FXT");// (year 2020)");
  prText->AddText("Proton");
  prText->AddText("0 < y_{c.m.} < 0.5");
  prText->SetFillColorAlpha(0,0);
  prText->SetLineColorAlpha(0,0);
  prText->SetTextFont(23);
  prText->SetTextSize(45);

  //TPaveText* prelimText = new TPaveText(0.3, 0.02, 1.28, 0.03, "NB");
  TPaveText* prelimText = new TPaveText(0.1, -0.06, 0.5, -0.04, "NB");
  prelimText->AddText("STAR");
  //prelimText->SetTextColor(kRed);
  prelimText->SetFillColorAlpha(0,0);
  prelimText->SetTextFont(23);
  prelimText->SetTextSize(45);

  TLine *zeroLine_pt = new TLine(0, 0, 2, 0);
  zeroLine_pt->SetLineStyle(9);
  ////

  if (useSystematics)
    {    
      //Removing x error bars
      for (int i = 0; i < 10; i++)
	{
	  sys_pT_00to10_pr->SetPointError(i, 0.0, sys_pT_00to10_pr->GetErrorY(i));
	  sys_pT_10to40_pr->SetPointError(i, 0.0, sys_pT_10to40_pr->GetErrorY(i));
	  sys_pT_40to60_pr->SetPointError(i, 0.0, sys_pT_40to60_pr->GetErrorY(i));
	}
    }

  prPtStack->Draw();
  prPtStack->GetXaxis()->SetLabelFont(133);
  prPtStack->GetYaxis()->SetLabelFont(133);
  prPtStack->GetXaxis()->SetLabelSize(55);
  prPtStack->GetYaxis()->SetLabelSize(55);
  prPtStack->GetYaxis()->SetTitleOffset(1.6);
  prPtStack->GetXaxis()->SetTitleOffset(0.95);
  prPtStack->GetXaxis()->SetTitleFont(133);
  prPtStack->GetYaxis()->SetTitleFont(133);
  prPtStack->GetXaxis()->SetTitleSize(55);
  prPtStack->GetYaxis()->SetTitleSize(55);
  prPtStack->GetXaxis()->SetNdivisions(505);
  prPtStack->Draw();
  prPtStack->SetMaximum(0.04);
  prPtStack->SetMinimum(-0.12);
  prPtStack->Draw("NOSTACK EP");
  zeroLine_pt->Draw("SAME");
  if (useSystematics)
    {
      sys_pT_00to10_pr->Draw("[]");
      sys_pT_40to60_pr->Draw("[]");
      sys_pT_10to40_pr->Draw("[]");
    }
  prPtStack->Draw("NOSTACK EP SAME");
  prLegend->Draw();
  prText->Draw();
  prelimText->Draw();
  canvas->SaveAs("v3_prPtStack.pdf");
  canvas->SaveAs("v3_prPtStack.png");
  canvas->Clear();

  if (useSystematics)
    systematicFile->Close();
  
  file->Close();
}
