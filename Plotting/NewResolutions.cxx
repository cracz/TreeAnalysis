// Bessel function resolution with m = 1 and k = 1
Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
    Double_t chi = x_val[0];
    Double_t arg = chi*chi/4.0;
    Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

    return norm * chi * TMath::Exp(-1.0*arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));
}
/*
// Bessel function resolution with m = 1 and k = 2
Double_t Resolution_12(Double_t *x_val, Double_t *par)
{
    Double_t y;
    Double_t chi = x_val[0];
    Double_t arg = chi*chi/4.0;
    Double_t norm = TMath::Sqrt(TMath::Pi())/(2.0*TMath::Sqrt(2.0));
    Double_t besselOneHalf = TMath::Sqrt(2.0*arg/TMath::Pi()) * TMath::SinH(arg)/arg;
    Double_t besselThreeHalf = TMath::Sqrt(2.0*arg/TMath::Pi()) * (TMath::CosH(arg)/arg - TMath::SinH(arg)/(arg*arg) );

    y = norm * chi * TMath::Exp(-1.0*arg) * (besselOneHalf + besselThreeHalf);

    return y;
}
*/
// Bessel function resolution with m = 1 and k = 3
Double_t Resolution_13(Double_t *x_val, Double_t *par)
{
  Double_t chi = x_val[0];
  Double_t arg = chi*chi/4.0;
  Double_t norm = TMath::Sqrt(TMath::Pi())/(2.0*TMath::Sqrt(2.0));

  return norm * chi * TMath::Exp(-1.0*arg) * (TMath::BesselI1(arg) + TMath::BesselI(2, arg));
}



void NewResolutions(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";
  //TString fileName = jobID + ".root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}
  
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 1000);
  canvas->SetTicks();
  canvas->SetGrid();
  canvas->SetTopMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.13);
  canvas->cd();
  
  TProfile *p_EpdAEpdB_R11 = (TProfile*)file->Get("p_EpdAEpdB_R11");
  TProfile *p_TpcBEpdA_R11 = (TProfile*)file->Get("p_TpcBEpdA_R11");
  TProfile *p_TpcBEpdB_R11 = (TProfile*)file->Get("p_TpcBEpdB_R11");
  /*
  TH1D *h_EpdAEpdB = p_EpdAEpdB_R11->ProjectionX();
  TH1D *h_TpcBEpdA = p_TpcBEpdA_R11->ProjectionX();
  TH1D *h_TpcBEpdB = p_TpcBEpdB_R11->ProjectionX();
  
  TH1D *h_EpdAEpdB_flip = PlotUtils::flipHisto(h_EpdAEpdB);
  TH1D *h_TpcBEpdA_flip = PlotUtils::flipHisto(h_TpcBEpdA);
  TH1D *h_TpcBEpdB_flip = PlotUtils::flipHisto(h_TpcBEpdB);
  */
  Int_t centBins    = p_EpdAEpdB_R11->GetNbinsX();
  Int_t firstCentID = p_EpdAEpdB_R11->GetBinLowEdge(1);
  Int_t lastCentID  = p_EpdAEpdB_R11->GetBinLowEdge(p_EpdAEpdB_R11->GetNbinsX());

  //// Calculate R11 in each centrality using the 3 sub-event equation
  std::vector<Double_t> v_R11;       // resolutions
  std::vector<Double_t> v_dR11;      // statistical uncertainties
  std::vector<Double_t> v_R11plus;   // resolutions + statistical uncertainties
  std::vector<Double_t> v_R11minus;  // resolutions - statistical uncertainties

  Double_t R11_epdA;
  Double_t dR11_epdA;

  Double_t EpdAEpdB;
  Double_t TpcBEpdA;
  Double_t TpcBEpdB;

  Double_t dEpdAEpdB;
  Double_t dTpcBEpdA;
  Double_t dTpcBEpdB;

  for (int i = 1; i <= centBins; i++)
    {
      EpdAEpdB = p_EpdAEpdB_R11->GetBinContent(i);
      TpcBEpdA = p_TpcBEpdA_R11->GetBinContent(i);
      TpcBEpdB = p_TpcBEpdB_R11->GetBinContent(i);

      dEpdAEpdB = p_EpdAEpdB_R11->GetBinError(i);
      dTpcBEpdA = p_TpcBEpdA_R11->GetBinError(i);
      dTpcBEpdB = p_TpcBEpdB_R11->GetBinError(i);

      R11_epdA  = TMath::Sqrt( (EpdAEpdB * TpcBEpdA) / TpcBEpdB );
      dR11_epdA = R11_epdA * TMath::Sqrt((dEpdAEpdB/(2*EpdAEpdB))*(dEpdAEpdB/(2*EpdAEpdB)) +
					 (dTpcBEpdA/(2*TpcBEpdA))*(dTpcBEpdA/(2*TpcBEpdA)) +
					 (dTpcBEpdB/(2*TpcBEpdB))*(dTpcBEpdB/(2*TpcBEpdB)));

      if(!TMath::IsNaN(R11_epdA))
	{
	  v_R11.push_back(R11_epdA);
	  v_dR11.push_back(dR11_epdA);
	  v_R11plus.push_back(R11_epdA + dR11_epdA);
	  v_R11minus.push_back(R11_epdA - dR11_epdA);
	}
      else
	{
	  v_R11.push_back(0.0);
	  v_dR11.push_back(0.0);
	  v_R11plus.push_back(0.0);
	  v_R11minus.push_back(0.0);
	}
    }
  //////////

  //   Extract Chi_1 for each centrality using the R11 values, then
  // plug those Chi_1 values into the R13 equation.
  std::vector<Double_t> v_R13;
  std::vector<Double_t> v_dR13plus;
  std::vector<Double_t> v_dR13minus;
  TF1 *f_BesselRes11 = new TF1("f_BesselRes11", Resolution_Full, 0, 10, 0);
  TF1 *f_BesselRes13 = new TF1("f_BesselRes13", Resolution_13, 0, 10, 0);

  for (int i = 0; i < v_R11.size(); i++)
    {
      if (v_R11.at(i) == 0.0)
	{
	  v_R13.push_back(0.0);
	  v_dR13plus.push_back(0.0);
	  v_dR13minus.push_back(0.0);
	  continue;
	}
      
      Double_t chi_1 = f_BesselRes11->GetX(v_R11.at(i));
      Double_t chi_1_plus  = f_BesselRes11->GetX(v_R11plus.at(i));
      Double_t chi_1_minus = f_BesselRes11->GetX(v_R11minus.at(i));

      Double_t R13   = f_BesselRes13->Eval(chi_1);
      Double_t R13plus  = f_BesselRes13->Eval(chi_1_plus);
      Double_t R13minus = f_BesselRes13->Eval(chi_1_minus);

      v_R13.push_back(R13);
      v_dR13plus.push_back(R13plus - R13);
      v_dR13minus.push_back(R13 - R13minus);
    }
  //////////

  /*
  for (int i = 0; i < v_R13.size(); i++)
    {
      std::cout << "Delta+ = " << v_R13plus.at(i) - v_R13.at(i) << std::endl;
      std::cout << "Delta- = " << v_R13.at(i) - v_R13minus.at(i) << std::endl;
      std::cout << std::endl;
    }
  */

  // Save R13 values for later use
  TH1D *h_resolutions = new TH1D("h_resolutions","EPD A Resolutions;Centrality;R_{31}",centBins,0,centBins);
  TH1D *h_resolutions_R11_3sub = new TH1D("h_resolutions_R11_3sub","EPD A 3-sub Resolutions;Centrality;R_{11}",centBins,0,centBins);

  for (int i = 1; i <= v_R13.size(); i++)
    {
      h_resolutions->SetBinContent(i, v_R13.at(i-1));
      h_resolutions->SetBinError(i, v_dR13plus.at(i-1));
      h_resolutions_R11_3sub->SetBinContent(i, v_R11.at(i-1));
      h_resolutions_R11_3sub->SetBinError(i, v_dR11.at(i-1));
    }

  TFile *resolutionInfo_INPUT = new TFile("resolutionInfo_INPUT.root", "RECREATE");
  h_resolutions->Write();
  h_resolutions_R11_3sub->Write();
  resolutionInfo_INPUT->Close();
  //////////
  

  // < C++11 initializations
  std::vector<Double_t> centralities;
  const double tmp_arr[] = {77.5, 72.5, 67.5, 62.5, 57.5, 52.5, 47.5, 42.5, 37.5, 32.5, 27.5, 22.5, 17.5, 12.5, 7.5, 2.5};
  for (unsigned int j = 0; j < (sizeof(tmp_arr) / sizeof(tmp_arr[0])); j++)
    centralities.push_back(tmp_arr[j]);

  // C++11 initializations
  //std::vector<Double_t> centralities = {77.5, 72.5, 67.5, 62.5, 57.5, 52.5, 47.5, 42.5, 37.5, 32.5, 27.5, 22.5, 17.5, 12.5, 7.5, 2.5};
  //std::vector<Double_t> centralities = {72.5, 67.5, 62.5, 57.5, 52.5, 47.5, 42.5, 37.5, 32.5, 27.5, 22.5, 17.5, 12.5, 7.5, 2.5};


  TGraphErrors *gr_R11 = new TGraphErrors(centBins);   // 3-sub
  TGraphErrors *gr_R13 = new TGraphErrors(centBins);   // full process

  Int_t graphPointIndex = 0;
  for (int i = v_R11.size()-1; i >= 0; i--)
    {
      gr_R11->SetPoint(graphPointIndex, centralities.at(i), v_R11.at(i));
      gr_R11->SetPointError(graphPointIndex, 0.0, v_dR11.at(i));
      graphPointIndex++;
    }

  graphPointIndex = 0;
  for (int i = v_R13.size()-1; i >= 0; i--)
    {
      gr_R13->SetPoint(graphPointIndex, centralities.at(i), v_R13.at(i));
      gr_R13->SetPointError(graphPointIndex, 0.0, v_dR13plus.at(i));
      graphPointIndex++;
    }

  
  //gr_R11->SetLineWidth(2);
  //gr_R11->SetLineColor(kBlack);
  gr_R11->SetMarkerStyle(20);
  gr_R11->SetMarkerSize(2);
  gr_R11->SetMarkerColor(kBlue);
  //gr_R11->SetFillColorAlpha(kBlue-4, 0.3);

  gr_R11->GetXaxis()->SetLabelSize(0.05);
  gr_R11->GetYaxis()->SetLabelSize(0.045);
  gr_R11->GetXaxis()->SetTitleOffset(1.1);
  gr_R11->GetYaxis()->SetTitleOffset(1.1);
  gr_R11->GetXaxis()->SetTitleSize(0.045);
  gr_R11->GetYaxis()->SetTitleSize(0.055);
  gr_R11->GetYaxis()->SetLabelSize(0.045);
  gr_R11->GetXaxis()->SetTitleFont(132);
  gr_R11->GetYaxis()->SetTitleFont(132);
  gr_R11->GetXaxis()->SetTitle("Centrality(%)");
  gr_R11->GetYaxis()->SetTitle("R_{11}");


  //gr_R13->SetLineWidth(2);
  //gr_R13->SetLineColor(kBlack);
  gr_R13->SetMarkerStyle(20);
  gr_R13->SetMarkerSize(2);
  gr_R13->SetMarkerColor(kBlue);
  //gr_R13->SetFillColorAlpha(kBlue-4, 0.3);

  gr_R13->GetXaxis()->SetLabelSize(0.05);
  gr_R13->GetYaxis()->SetLabelSize(0.045);
  gr_R13->GetXaxis()->SetTitleOffset(1.1);
  gr_R13->GetYaxis()->SetTitleOffset(1.1);
  gr_R13->GetXaxis()->SetTitleSize(0.045);
  gr_R13->GetYaxis()->SetTitleSize(0.055);
  gr_R13->GetYaxis()->SetLabelSize(0.045);
  gr_R13->GetXaxis()->SetTitleFont(132);
  gr_R13->GetYaxis()->SetTitleFont(132);
  gr_R13->GetXaxis()->SetTitle("Centrality(%)");
  gr_R13->GetYaxis()->SetTitle("R_{13}");

  //gr_R11->SetMaximum(0.25);
  //gr_R11->SetMaximum(1.0);
  //gr_R11->SetMinimum(0.0);
  gr_R11->SetTitle("");
  gr_R11->Draw("ACP");

  canvas->SaveAs(jobID + "_resolutionAonly_R11.pdf");
  canvas->Clear();

  gr_R13->SetTitle("");
  gr_R13->Draw("ACP");

  canvas->SaveAs(jobID + "_resolutionAonly_R13.pdf");
  canvas->Clear();
  
}
