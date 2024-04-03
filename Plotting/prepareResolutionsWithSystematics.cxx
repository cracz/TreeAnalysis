#include "PlotUtils.h"

double maxValue(double d1, double d2, double d3)
{
  double max = d1;
  if (d2 > max) max = d2;
  if (d3 > max) max = d3;
  return max;
}

void prepareResolutionsWithSystematics()
{
  /*
  // 3.0 GeV
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p0GeV_EPDB9to16.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p0GeV_EPDB10to16.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p0GeV_9to16scaled.root", "READ");
  */
  /*
  // 3.0 GeV v1
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v1_EPDB9to16.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v1_EPDB10to16.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v1_EPDB13to16.root", "READ");
  */
  /*
  // 3.0 GeV v2
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v2_EPDB9to16.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v2_EPDB10to16.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v2_EPDB13to16.root", "READ");
  */
  /*
  // 3.2 GeV
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p2GeV_EPDA1to6_EPDB7to13_withEff_SL23d.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p2GeV_EPDA1to6_EPDB8to13_withEff_SL23d.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p2GeV_EPDA1to6_EPDB9to13_withEff_SL23d.root", "READ");
  */
  /*  
  // 3.5 GeV 
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p5GeV_EPDA1to6_EPDB7to11_withEff_SL23d.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p5GeV_EPDA1to6_EPDB8to11_withEff_SL23d.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p5GeV_EPDA1to6_EPDB9to11_withEff_SL23d.root", "READ");
  */
  /*
  // 3.9 GeV
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p9GeV_EPDA1to5_EPDB6to10_withEff_SL23e.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p9GeV_EPDA1to5_EPDB7to10_withEff_SL23e.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p9GeV_EPDA1to5_EPDB8to10_withEff_SL23e.root", "READ");
  */

  // 4.5 GeV
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_4p5GeV_EPDA1to3_EPDB4to9_SL23e.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_4p5GeV_EPDA1to3_EPDB5to9_SL23e.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_4p5GeV_EPDA1to3_EPDB6to9_SL23e.root", "READ");

  
  if (!resVarFile1) std::cout << "No resVarFile1" << std::endl;
  if (!resVarFile2) std::cout << "No resVarFile2" << std::endl;
  if (!resVarFile3) std::cout << "No resVarFile3" << std::endl;

  TH1D* h_resVar1 = (TH1D*)resVarFile1->Get("h_resolutions");
  h_resVar1->SetName("h_resVar1");
  TH1D* h_resVar2 = (TH1D*)resVarFile2->Get("h_resolutions");
  h_resVar2->SetName("h_resVar2");
  TH1D* h_resVar3 = (TH1D*)resVarFile3->Get("h_resolutions");
  h_resVar3->SetName("h_resVar3");

  // Calculate the average resolution
  TH1D* h_avgRes = (TH1D*)h_resVar1->Clone("h_resolutions");
  h_avgRes->Add(h_resVar2);
  h_avgRes->Add(h_resVar3);
  h_avgRes->Scale(1.0/3.0);
  ////

  // Make a copy of h_avgRes to use with plotting
  TH1D* h_resolutionsWithStats = (TH1D*)h_avgRes->Clone("h_resolutionsWithStats");
  h_resolutionsWithStats->SetName("h_resolutionsWithStats");
  ////

  // Make a copy of h_avgRes that has systematic uncertainties instead of statistical uncertainties for plotting
  TH1D* h_resolutionsWithSysts = (TH1D*)h_avgRes->Clone("h_resolutionsWithSysts");
  h_resolutionsWithSysts->SetName("h_resolutionsWithSysts");
  
  for (int ithBin = 1; ithBin < h_resolutionsWithSysts->GetNbinsX(); ithBin++)
    {
      double diff1 = TMath::Abs(h_resolutionsWithSysts->GetBinContent(ithBin) - h_resVar1->GetBinContent(ithBin));
      double diff2 = TMath::Abs(h_resolutionsWithSysts->GetBinContent(ithBin) - h_resVar2->GetBinContent(ithBin));
      double diff3 = TMath::Abs(h_resolutionsWithSysts->GetBinContent(ithBin) - h_resVar3->GetBinContent(ithBin));
      double maxDiff = maxValue(diff1, diff2, diff3);
      h_resolutionsWithSysts->SetBinError(ithBin, maxDiff);
    }
  ////


  // Make a copy of h_avgRes that has statistical and systematic uncertainties combined in quadrature.
  TH1D* h_resolutionsCombinedError = (TH1D*)h_avgRes->Clone("h_resolutionsCombinedError");
  h_resolutionsCombinedError->SetName("h_resolutionsCombinedError");
  
  for (int ithBin = 1; ithBin < h_resolutionsCombinedError->GetNbinsX(); ithBin++)
    {
      Double_t statError = h_resolutionsWithStats->GetBinError(ithBin);
      Double_t systError = h_resolutionsWithSysts->GetBinError(ithBin);
      Double_t newError = TMath::Sqrt( TMath::Power(statError, 2.0) + TMath::Power(systError, 2.0) );
      h_resolutionsCombinedError->SetBinError(ithBin, newError);
    }
  ////


  // Flip and trim axis on some plots
  h_resolutionsWithStats = PlotUtils::flipHisto(h_resolutionsWithStats);
  //h_resolutionsWithStats = PlotUtils::trimCentralityPlot(h_resolutionsWithStats);      // Out to 60%
  h_resolutionsWithStats = PlotUtils::trimCentralityPlotStrict(h_resolutionsWithStats);  // Out to 40%

  h_resolutionsWithSysts = PlotUtils::flipHisto(h_resolutionsWithSysts);
  //h_resolutionsWithSysts = PlotUtils::trimCentralityPlot(h_resolutionsWithSysts);     // Out to 60%
  h_resolutionsWithSysts = PlotUtils::trimCentralityPlotStrict(h_resolutionsWithSysts); // Out to 40%
  ////

  TFile* newFile = TFile::Open("eventPlaneSystematics.root", "RECREATE");
  newFile->cd();
  h_avgRes->Write();
  h_resolutionsWithStats->Write();
  h_resolutionsWithSysts->Write();
  h_resolutionsCombinedError->Write();
  h_resVar1->Write();
  h_resVar2->Write();
  h_resVar3->Write();
  newFile->Close();

  resVarFile1->Close();
  resVarFile2->Close();
  resVarFile3->Close();

  std::cout << "Results saved in eventPlaneSystematics.root" << std::endl;
}
