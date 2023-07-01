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
  TFile* resVarFile1 = TFile::Open("v3_EPDref9to16_resolutionInfo.root", "READ");
  TFile* resVarFile2 = TFile::Open("v3_EPDref10to16_resolutionInfo.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_scaledByRes.root", "READ");
  */
  /*
  // 3.0 GeV v1
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v1_EPDB9to16.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v1_EPDB10to16.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v1_EPDB13to16.root", "READ");
  */

  // 3.0 GeV v2
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v2_EPDB9to16.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v2_EPDB10to16.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p0GeV_v2_EPDB13to16.root", "READ");

  /*
  // 3.2 GeV
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p22GeV_EPDA1to6_EPDB7to13_TPCB1p1to0.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p22GeV_EPDA1to6_EPDB8to13_TPCB1p1to0.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p22GeV_EPDA1to6_EPDB9to13_TPCB1p1to0.root", "READ");
  */
  /*
  // 3.5 GeV
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p5GeV_EPDA1to6_EPDref7to11_TPCB1p2to0.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p5GeV_EPDA1to6_EPDref8to11_TPCB1p2to0.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p5GeV_EPDA1to6_EPDref9to11_TPCB1p2to0.root", "READ");
  */
  /*
  // 3.9 GeV
  TFile* resVarFile1 = TFile::Open("resolutionInfo_INPUT_3p9GeV_EPDA1to5_EPDref6to10_TPCB1p32to0.root", "READ");
  TFile* resVarFile2 = TFile::Open("resolutionInfo_INPUT_3p9GeV_EPDA1to5_EPDref7to10_TPCB1p32to0.root", "READ");
  TFile* resVarFile3 = TFile::Open("resolutionInfo_INPUT_3p9GeV_EPDA1to5_EPDref8to10_TPCB1p32to0.root", "READ");
  */
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
  h_resolutionsWithStats = PlotUtils::flipHisto(h_resolutionsWithStats);
  h_resolutionsWithStats = PlotUtils::trimCentralityPlot(h_resolutionsWithStats);      // Out to 60%
  //h_resolutionsWithStats = PlotUtils::trimCentralityPlotStrict(h_resolutionsWithStats);  // Out to 40%
  h_resolutionsWithStats->SetName("h_resolutionsWithStats");
  ////

  // Make a copy of h_avgRes that has systematic uncertainties instead of statistical uncertainties for plotting
  TH1D* h_resolutionsWithSysts = (TH1D*)h_avgRes->Clone("h_resolutionsWithSysts");

  for (int ithBin = 1; ithBin < h_resolutionsWithSysts->GetNbinsX(); ithBin++)
    {
      double diff1 = TMath::Abs(h_resolutionsWithSysts->GetBinContent(ithBin) - h_resVar1->GetBinContent(ithBin));
      double diff2 = TMath::Abs(h_resolutionsWithSysts->GetBinContent(ithBin) - h_resVar2->GetBinContent(ithBin));
      double diff3 = TMath::Abs(h_resolutionsWithSysts->GetBinContent(ithBin) - h_resVar3->GetBinContent(ithBin));
      double maxDiff = maxValue(diff1, diff2, diff3);

      h_resolutionsWithSysts->SetBinError(ithBin, maxDiff);
    }
  h_resolutionsWithSysts = PlotUtils::flipHisto(h_resolutionsWithSysts);
  h_resolutionsWithSysts = PlotUtils::trimCentralityPlot(h_resolutionsWithSysts);       // Out to 60%
  //h_resolutionsWithSysts = PlotUtils::trimCentralityPlotStrict(h_resolutionsWithSysts); // Out to 40%
  h_resolutionsWithSysts->SetName("h_resolutionsWithSysts");
  ////

  TFile* newFile = TFile::Open("eventPlaneSystematics.root", "RECREATE");
  newFile->cd();
  h_avgRes->Write();
  h_resolutionsWithStats->Write();
  h_resolutionsWithSysts->Write();
  h_resVar1->Write();
  h_resVar2->Write();
  h_resVar3->Write();
  newFile->Close();

  resVarFile1->Close();
  resVarFile2->Close();
  resVarFile3->Close();

  std::cout << "Results saved in eventPlaneSystematics.root" << std::endl;
}
