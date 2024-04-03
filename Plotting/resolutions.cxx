#include "PlotUtils.h"

void resolutions(TString jobID, TString order_n_str)
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
  TH1D *h_resolEPDA = new TH1D("h_resolEPDA","EPD A vs EPD B and TPC B;Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);
  TH1D *h_resolutions = new TH1D("h_resolutions","EPD A Resolutions;Centrality;R_{"+order_n_str+"1}",centBins,0,centBins);
  
  
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

  Double_t R_AvsB_save;
  Double_t dR_AvsB_save;
      
  // Fill resolution plots
  for (int i = 1; i <= centBins; i++)
    {
      EpdAEpdB_save = h_EpdAEpdB->GetBinContent(i);  //Don't use the flipped values here in the saved histogram!      
      TpcBEpdA_save = h_TpcBEpdA->GetBinContent(i);  // We need the centrality ID's in order, not the centrality percentages.
      TpcBEpdB_save = h_TpcBEpdB->GetBinContent(i);

      dEpdAEpdB_save = h_EpdAEpdB->GetBinError(i);
      dTpcBEpdA_save = h_TpcBEpdA->GetBinError(i);
      dTpcBEpdB_save = h_TpcBEpdB->GetBinError(i);

      R_AvsB_save  = TMath::Sqrt( (EpdAEpdB_save * TpcBEpdA_save) / TpcBEpdB_save );
      dR_AvsB_save = R_AvsB_save * TMath::Sqrt((dEpdAEpdB_save/(2*EpdAEpdB_save))*(dEpdAEpdB_save/(2*EpdAEpdB_save)) +
					       (dTpcBEpdA_save/(2*TpcBEpdA_save))*(dTpcBEpdA_save/(2*TpcBEpdA_save)) +
					       (dTpcBEpdB_save/(2*TpcBEpdB_save))*(dTpcBEpdB_save/(2*TpcBEpdB_save)));
      /*
      if (i == 5)
	{
	  std::cout << "( (EpdAEpdB_save * TpcBEpdA_save) / TpcBEpdB_save )^2 = "
		    << "(" << EpdAEpdB_save << " * " << TpcBEpdA_save << ") / " << TpcBEpdB_save
		    << std::endl << std::endl;

	  std::cout << "dEpdAEpdB_save = " << dEpdAEpdB_save << std::endl;
	  std::cout << "dTpcBEpdA_save = " << dTpcBEpdA_save << std::endl;
	  std::cout << "dTpcBEpdB_save = " << dTpcBEpdB_save << std::endl;
	}
      */
      
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

      h_resolEPDA->SetBinContent(i, R_AvsB);
      h_resolEPDA->SetBinError(i, dR_AvsB);

      if(!TMath::IsNaN(R_AvsB_save))
	{
	  h_resolutions->SetBinContent(i, R_AvsB_save);
	  h_resolutions->SetBinError(i, dR_AvsB_save);
	}
    }

  h_resolutions->Write();
  
  
  gStyle->SetOptStat(0);

  h_resolEPDA->SetMarkerStyle(20);
  h_resolEPDA->SetMarkerSize(1.5);
  h_resolEPDA->SetMarkerColor(kBlue-7);
  h_resolEPDA->SetLineColor(kBlue-7);
  
  TLegend *legend2 = new TLegend(0.65, 0.82, 0.96, 0.96);
  legend2->AddEntry(h_resolEPDA,"Inner EPD #psi_{1}");
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

  h_resolEPDA = PlotUtils::trimCentralityPlot(h_resolEPDA);

  h_resolEPDA->SetBinContent(9,0.0);
  h_resolEPDA->SetBinError(9,0.0);
  h_resolEPDA->SetBinContent(10,0.0);
  h_resolEPDA->SetBinError(10,0.0);
  h_resolEPDA->SetBinContent(11,0.0);
  h_resolEPDA->SetBinError(11,0.0);

  h_resolEPDA->SetLineWidth(2);
  h_resolEPDA->SetLineColor(kBlack);
  h_resolEPDA->SetMarkerStyle(20);
  h_resolEPDA->SetMarkerSize(2);
  h_resolEPDA->SetMarkerColor(kBlue);
  h_resolEPDA->SetFillColorAlpha(kBlue-4, 0.3);

  h_resolEPDA->GetXaxis()->SetLabelSize(0.05);
  h_resolEPDA->GetYaxis()->SetLabelSize(0.045);
  h_resolEPDA->GetXaxis()->SetTitleOffset(1.1);
  h_resolEPDA->GetYaxis()->SetTitleOffset(1.1);
  h_resolEPDA->GetXaxis()->SetTitleSize(0.045);
  h_resolEPDA->GetYaxis()->SetTitleSize(0.055);
  h_resolEPDA->GetYaxis()->SetLabelSize(0.045);
  h_resolEPDA->GetXaxis()->SetTitleFont(132);
  h_resolEPDA->GetYaxis()->SetTitleFont(132);
  h_resolEPDA->SetMaximum(0.25);
  //h_resolEPDA->SetMaximum(1.0);
  h_resolEPDA->SetMinimum(0.0);
  h_resolEPDA->SetTitle("");
  //h_resolEPDA->SetMarkerColor(1);
  //h_resolEPDA->SetLineColor(1);
  h_resolEPDA->Draw("E1P");
  //legend2->Draw();
  text->Draw();
  //prelimText->Draw();
  canvas->SaveAs(jobID + "_resolutionAonly.pdf");
  canvas->Clear();

  resolutionInfo_INPUT->Close();
  file->Close();
}
