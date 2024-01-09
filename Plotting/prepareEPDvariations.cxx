#include "PlotUtils.h"

void prepareEPDvariations()
{
  TFile* sourceFile = TFile::Open("eventPlaneSystematics_3p0GeV.root","READ");
  //TFile* sourceFile = TFile::Open("eventPlaneSystematics_3p2GeV_1to6_max13_SL23d.root","READ");
  //TFile* sourceFile = TFile::Open("eventPlaneSystematics_3p5GeV_1to6_max11.root","READ");
  //TFile* sourceFile = TFile::Open("eventPlaneSystematics_3p9GeV_1to5_max10.root","READ");

  TH1D* h_resolutionsCombinedError = (TH1D*)sourceFile->Get("h_resolutionsCombinedError");

  // EPD low settings
  TH1D* h_resolutions = (TH1D*)h_resolutionsCombinedError->Clone("h_resolutions");

  for (int bin = 1; bin <= h_resolutions->GetNbinsX(); bin++)
    {
      h_resolutions->SetBinContent(bin, h_resolutions->GetBinContent(bin) - h_resolutions->GetBinError(bin));
    }

  TFile* lowSettingOutputFile = TFile::Open("resolutionInfo_INPUT_epd_low.root","RECREATE");
  lowSettingOutputFile->cd();
  h_resolutions->Write();
  lowSettingOutputFile->Close();
  ////
  
  // EPD high settings
  h_resolutions = (TH1D*)h_resolutionsCombinedError->Clone("h_resolutions");

  for (int bin = 1; bin <= h_resolutions->GetNbinsX(); bin++)
    {
      h_resolutions->SetBinContent(bin, h_resolutions->GetBinContent(bin) + h_resolutions->GetBinError(bin));
    }

  TFile* highSettingOutputFile = TFile::Open("resolutionInfo_INPUT_epd_high.root","RECREATE");
  highSettingOutputFile->cd();
  h_resolutions->Write();
  highSettingOutputFile->Close();
  ////

  sourceFile->Close();

  std::cout << "Created resolutionInfo_INPUT_epd_low.root and resolutionInfo_INPUT_epd_high.root." << std::endl;
}
