#include "PlotUtils.h"

void resolutionsEPDB(TString jobID, TString order_n_str)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";
  //TString fileName = jobID + ".root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TFile *resolutionInfo_INPUT = new TFile("resolutionInfo_INPUT.root", "RECREATE");
  /*
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 875, 675);
  canvas->SetGridx();
  canvas->SetGridy();
  //canvas->SetLeftMargin(0.15);
  canvas->cd();
  */
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 1000);
  canvas->SetTicks();
  canvas->SetGrid();
  canvas->SetTopMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.13);
  canvas->cd();
  
  TProfile *p_EpdAEpdB = (TProfile*)file->Get("p_EpdAEpdB");
  TProfile *p_TpcBEpdA = (TProfile*)file->Get("p_TpcBEpdA");
  TProfile *p_TpcBEpdB = (TProfile*)file->Get("p_TpcBEpdB");

  TH1D *h_EpdAEpdB = p_EpdAEpdB->ProjectionX();
  TH1D *h_TpcBEpdA = p_TpcBEpdA->ProjectionX();
  TH1D *h_TpcBEpdB = p_TpcBEpdB->ProjectionX();
  
  TH1D *h_EpdAEpdB_flip = PlotUtils::flipHisto(h_EpdAEpdB);
  TH1D *h_TpcBEpdA_flip = PlotUtils::flipHisto(h_TpcBEpdA);
  TH1D *h_TpcBEpdB_flip = PlotUtils::flipHisto(h_TpcBEpdB);

  Int_t centBins    = h_EpdAEpdB->GetNbinsX();
  Int_t firstCentID = h_EpdAEpdB->GetBinLowEdge(1);
  Int_t lastCentID  = h_EpdAEpdB->GetBinLowEdge(h_EpdAEpdB->GetNbinsX());

  // Make plots of EPD A resolutions.
  TH1D *h_resolEPDB = new TH1D("h_resolEPDB","EPD B vs EPD A and TPC B;Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);
  TH1D *h_resolutions = new TH1D("h_resolutions","EPD B Resolutions;Centrality;R_{"+order_n_str+"1}",centBins,0,centBins);
  
  
  Double_t EpdAEpdB;
  Double_t TpcBEpdA;
  Double_t TpcBEpdB;

  Double_t dEpdAEpdB;
  Double_t dTpcBEpdA;
  Double_t dTpcBEpdB;

  Double_t R_AvsB;
  Double_t R_BvsA;
  Double_t R_TpcB;
  Double_t dR_AvsB;
  Double_t dR_BvsA;
  Double_t dR_TpcB;

  Double_t EpdAEpdB_save;
  Double_t TpcBEpdA_save;
  Double_t TpcBEpdB_save;

  Double_t dEpdAEpdB_save;
  Double_t dTpcBEpdA_save;
  Double_t dTpcBEpdB_save;

  Double_t R_BvsA_save;
  Double_t dR_BvsA_save;
      
  // Fill resolution plots
  for (int i = 1; i <= centBins; i++)
    {
      EpdAEpdB_save = h_EpdAEpdB->GetBinContent(i);  //Don't use the flipped values here in the saved histogram!      
      TpcBEpdA_save = h_TpcBEpdA->GetBinContent(i);  // We need the centrality ID's in order, not the centrality percentages.
      TpcBEpdB_save = h_TpcBEpdB->GetBinContent(i);

      dEpdAEpdB_save = h_EpdAEpdB->GetBinError(i);
      dTpcBEpdA_save = h_TpcBEpdA->GetBinError(i);
      dTpcBEpdB_save = h_TpcBEpdB->GetBinError(i);

      R_BvsA_save  = TMath::Sqrt( (EpdAEpdB_save * TpcBEpdB_save) / TpcBEpdA_save );
      dR_BvsA_save = R_BvsA_save * TMath::Sqrt((dEpdAEpdB_save/(2*EpdAEpdB_save))*(dEpdAEpdB_save/(2*EpdAEpdB_save)) +
					       (dTpcBEpdA_save/(2*TpcBEpdA_save))*(dTpcBEpdA_save/(2*TpcBEpdA_save)) +
					       (dTpcBEpdB_save/(2*TpcBEpdB_save))*(dTpcBEpdB_save/(2*TpcBEpdB_save)));
      
      EpdAEpdB = h_EpdAEpdB_flip->GetBinContent(i);      
      TpcBEpdA = h_TpcBEpdA_flip->GetBinContent(i);
      TpcBEpdB = h_TpcBEpdB_flip->GetBinContent(i);

      dEpdAEpdB = h_EpdAEpdB_flip->GetBinError(i);
      dTpcBEpdA = h_TpcBEpdA_flip->GetBinError(i);
      dTpcBEpdB = h_TpcBEpdB_flip->GetBinError(i);

      R_AvsB = TMath::Sqrt( (EpdAEpdB * TpcBEpdA) / TpcBEpdB );
      R_BvsA = TMath::Sqrt( (EpdAEpdB * TpcBEpdB) / TpcBEpdA );
      R_TpcB = TMath::Sqrt( (TpcBEpdA * TpcBEpdB) / EpdAEpdB );

      dR_AvsB = R_AvsB * TMath::Sqrt((dEpdAEpdB/(2*EpdAEpdB))*(dEpdAEpdB/(2*EpdAEpdB)) +
				     (dTpcBEpdA/(2*TpcBEpdA))*(dTpcBEpdA/(2*TpcBEpdA)) +
				     (dTpcBEpdB/(2*TpcBEpdB))*(dTpcBEpdB/(2*TpcBEpdB)));

      dR_BvsA = R_BvsA * TMath::Sqrt((dEpdAEpdB/(2*EpdAEpdB))*(dEpdAEpdB/(2*EpdAEpdB)) +
				     (dTpcBEpdA/(2*TpcBEpdA))*(dTpcBEpdA/(2*TpcBEpdA)) +
				     (dTpcBEpdB/(2*TpcBEpdB))*(dTpcBEpdB/(2*TpcBEpdB)));

      dR_TpcB = R_TpcB * TMath::Sqrt((dTpcBEpdA/(2*TpcBEpdA))*(dTpcBEpdA/(2*TpcBEpdA)) +
				     (dTpcBEpdB/(2*TpcBEpdB))*(dTpcBEpdB/(2*TpcBEpdB)) +
				     (dEpdAEpdB/(2*EpdAEpdB))*(dEpdAEpdB/(2*EpdAEpdB)));


      if(TMath::IsNaN(R_AvsB)) { R_AvsB = 0; dR_AvsB = 0; }
      if(TMath::IsNaN(R_BvsA)) { R_BvsA = 0; dR_BvsA = 0; }
      if(TMath::IsNaN(R_TpcB)) { R_TpcB = 0; dR_TpcB = 0; }

      h_resolEPDB->SetBinContent(i, R_BvsA);
      h_resolEPDB->SetBinError(i, dR_BvsA);

      if(!TMath::IsNaN(R_BvsA_save))
	{
	  h_resolutions->SetBinContent(i, R_BvsA_save);
	  h_resolutions->SetBinError(i, dR_BvsA_save);
	}
    }

  h_resolutions->Write();
  
  
  gStyle->SetOptStat(0);

  h_resolEPDB->SetMarkerStyle(20);
  h_resolEPDB->SetMarkerSize(1.5);
  h_resolEPDB->SetMarkerColor(kBlue-7);
  h_resolEPDB->SetLineColor(kBlue-7);
  
  TLegend *legend2 = new TLegend(0.65, 0.82, 0.96, 0.96);
  legend2->AddEntry(h_resolEPDB,"Inner EPD #psi_{1}");
  legend2->SetTextSize(0.05);

  TPaveText *text_extra = new TPaveText(2, 0.25, 28, 0.3);
  //text_extra->AddText("#sqrt{s_{NN}} = 3.0 GeV FXT Au+Au");
  text_extra->AddText("#sqrt{s_{NN}} = 3.2 GeV FXT Au+Au");
  text_extra->AddText("Collisions at RHIC");
  //text_extra->SetFillColorAlpha(0,0);

  TPaveText* prelimText = new TPaveText(5, 0.2, 25, 0.23, "NB");
  prelimText->AddText("STAR Preliminary");
  prelimText->SetTextColor(kRed);
  prelimText->SetFillColorAlpha(0,0);
  prelimText->SetTextSize(0.04);


  canvas->SetTicks();
  canvas->SetLogy(0);
  canvas->SetTopMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.13);

  gStyle->SetErrorX(0);

  TPaveText *text = new TPaveText(2, 0.215, 30, 0.24, "NB");
  text->AddText("STAR Au+Au #sqrt{s_{NN}} = 4.5 GeV FXT");
  text->SetFillColorAlpha(0,0);
  text->SetLineColorAlpha(0,0);
  text->SetTextSize(0.045);
  text->SetTextAlign(13);
  text->SetTextFont(22);

  h_resolEPDB = PlotUtils::trimCentralityPlot(h_resolEPDB);
  /*
  h_resolEPDB->SetBinContent(9,0.0);
  h_resolEPDB->SetBinError(9,0.0);
  h_resolEPDB->SetBinContent(10,0.0);
  h_resolEPDB->SetBinError(10,0.0);
  h_resolEPDB->SetBinContent(11,0.0);
  h_resolEPDB->SetBinError(11,0.0);
  h_resolEPDB->SetBinContent(12,0.0);
  h_resolEPDB->SetBinError(12,0.0);
  */
  h_resolEPDB->SetLineWidth(2);
  h_resolEPDB->SetLineColor(kBlack);
  h_resolEPDB->SetMarkerStyle(20);
  h_resolEPDB->SetMarkerSize(2);
  h_resolEPDB->SetMarkerColor(kBlue);
  h_resolEPDB->SetFillColorAlpha(kBlue-4, 0.3);

  h_resolEPDB->GetXaxis()->SetLabelSize(0.05);
  h_resolEPDB->GetYaxis()->SetLabelSize(0.045);
  h_resolEPDB->GetXaxis()->SetTitleOffset(1.1);
  h_resolEPDB->GetYaxis()->SetTitleOffset(1.1);
  h_resolEPDB->GetXaxis()->SetTitleSize(0.045);
  h_resolEPDB->GetYaxis()->SetTitleSize(0.055);
  h_resolEPDB->GetYaxis()->SetLabelSize(0.045);
  h_resolEPDB->GetXaxis()->SetTitleFont(132);
  h_resolEPDB->GetYaxis()->SetTitleFont(132);
  //h_resolEPDB->SetMaximum(0.25);
  //h_resolEPDB->SetMaximum(1.0);
  h_resolEPDB->SetMinimum(0.0);
  h_resolEPDB->SetTitle("");
  //h_resolEPDB->SetMarkerColor(1);
  //h_resolEPDB->SetLineColor(1);
  h_resolEPDB->Draw("E1P");
  //legend2->Draw();
  //   text->Draw();
  //prelimText->Draw();
  canvas->SaveAs(jobID + "_resolutionBonly.pdf");
  canvas->Clear();

  resolutionInfo_INPUT->Close();
  file->Close();
}
