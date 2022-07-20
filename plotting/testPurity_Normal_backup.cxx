double getIntersection(double amp1, double mean1, double std1, double amp2, double mean2, double std2)
{
  // First is background, second is signal!
  
  double a = std::pow(std2, 2) - std::pow(std1, 2);
  double b = 2 * mean2 * std::pow(std1, 2) - 2 * mean1 * std::pow(std2, 2);
  double c = std::pow(mean1 * std2, 2) - std::pow(mean2 * std1, 2) - 2 * std::pow(std1 * std2, 2) * std::log(amp1/amp2);
  /*
  std::cout << std::endl
	    << "a = " << a << std::endl
	    << "b = " << b << std::endl
	    << "c = " << c << std::endl
	    << std::endl;
  */
  double disc = std::pow(b, 2) - 4 * a * c;
  if (disc < 0) return -999.0;

  double posSolution =  (-b + std::sqrt(disc)) / (2 * a);
  double negSolution =  (-b - std::sqrt(disc)) / (2 * a);

  return (mean2 > mean1) ? posSolution : negSolution;
}

void printHistos(TH2D *h2, TString species)
{
  TCanvas* canvas = new TCanvas("canvas","canvas",1400,2400);
  canvas->Divide(3, 6);
  canvas->SaveAs("allPlots_"+species+".pdf[");
  int iPad = 0;
  
  TH1D *h_nSig;
  int nBins = h2->GetNbinsX();
  int nPlots = 0;
  Double_t pLow;
  Double_t pHigh;
  TString pLowStr;
  TString pHighStr;
  TString pBinStr;
  int pBinStart = (species == "de") ? 5/*9*/ : 11;

  TPaveText *text;
  TLine *lineLow;
  TLine *lineHigh;

  // fullFunction_de = low bkg + signal + high bkg
  TF1 *fullFunction_de =
    new TF1("fullFunction_de", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2) + [6]*exp(-0.5*((x-[7])/[8])**2)", -0.7, 0.7);

  TF1 *fullFunction_tr =
    new TF1("fullFunction_tr", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2) + [6]*exp(-0.5*((x-[7])/[8])**2)", -0.7, 0.7);

  /*
  TF1 *fullFunction_tr =
    new TF1("fullFunction_tr", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2) + [6]*exp(-0.5*((x-[7])/[8])**2) + [9]*exp(-0.5*((x-[10])/[11])**2)", -0.7, 0.7);
  */
  /*
  TF1 *fullFunction_tr =
    new TF1("fullFunction_tr", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2) + [6]*exp(-0.5*((x-[7])/[8])**2) + [9]*exp(-0.5*((x-[10])/[11])**2) + [12]*exp(-0.5*((x-[13])/[14])**2)", -0.7, 0.7);
  */
  // Low bkg for tritons
  TF1 *bkg0 = new TF1("bkg0", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)", -0.7, 0.7);
  bkg0->SetLineColor(4);
  bkg0->SetLineStyle(2);

  // Low single bkg
  TF1 *bkg1 = new TF1("bkg1", "[0]*exp(-0.5*((x-[1])/[2])**2)", -0.7, 0.7);
  bkg1->SetLineColor(4);
  bkg1->SetLineStyle(2);

  TF1 *signal = new TF1("signal", "[0]*exp(-0.5*((x-[1])/[2])**2)", -0.7, 0.7);
  signal->SetLineColor(8);
  signal->SetLineStyle(2);

  // High single bkg
  TF1 *bkg2 = new TF1("bkg2", "[0]*exp(-0.5*((x-[1])/[2])**2)", -0.7, 0.7);
  bkg2->SetLineColor(6);
  bkg2->SetLineStyle(2);

  // High double bkg
  TF1 *bkg3 = new TF1("bkg3", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)", -0.7, 0.7);
  bkg3->SetLineColor(6);
  bkg3->SetLineStyle(2);


  Double_t signalIntegral = 0.0;
  Double_t fullIntegral = 0.0;
  Double_t purity = 0.0;
  TString purityStr;

  Double_t lowBound = -0.2; // Default values
  Double_t highBound = 0.2;

  std::vector<int> binsUsed;
  std::vector<double> lowBounds;
  std::vector<double> highBounds;
  std::vector<double> purities;
  
  for (int pBin = /*16*/pBinStart; pBin <= nBins; pBin++)
    {
      pBinStr.Form("%d", pBin);
      h_nSig = h2->ProjectionY("h_"+species+"_bin"+pBinStr, pBin, pBin);
      h_nSig->SetMinimum(1.0);

      pLow  = h2->GetXaxis()->GetBinLowEdge(pBin);
      pHigh = pLow + h2->GetXaxis()->GetBinWidth(pBin);
      pLowStr.Form("%1.1f", pLow);
      pHighStr.Form("%1.1f", pHigh);

      h_nSig->SetTitle(pLowStr+" #leq |p| < "+pHighStr);
      h_nSig->GetYaxis()->SetTitle("Counts");
      h_nSig->GetYaxis()->SetTitleOffset(0.9);
      h_nSig->SetTitleSize(0.06, "x");
      h_nSig->SetTitleSize(0.06, "y");
      h_nSig->SetLabelSize(0.06, "x");
      h_nSig->SetLabelSize(0.06, "y");

      gStyle->SetTitleSize(0.09, "t");
      gStyle->SetStatW(0.25);
      gStyle->SetStatH(0.25);
      
      if (h_nSig->GetEntries() == 0) continue;
      
      canvas->cd(++iPad);
      //gPad->SetTopMargin(0);
      gPad->SetBottomMargin(0.2);
      gPad->SetLogy();

      text = new TPaveText(0.65, 0.85, 0.89, 0.9, "NB NDC");
      text->SetFillColorAlpha(0,0);

      //lineLow = new TLine(lowBound, 0.0, lowBound, h_nSig->GetMaximum());    
      //lineHigh = new TLine(highBound, 0.0, highBound, h_nSig->GetMaximum());

      // Only set params for first fit, for others just let
      //the previous values carry over.
      if (species == "de" && pBin == 5/*9*/)
	{
	  fullFunction_de->SetParameter(0, 2e+8); // Amps
	  fullFunction_de->SetParameter(3, 2.3e+06);
	  fullFunction_de->SetParameter(6, 4.88e+05);

	  fullFunction_de->SetParameter(2, 0.5);   // Std devs
	  //fullFunction_de->SetParLimits(2, 0.1/(2*TMath::Sqrt(2*TMath::Log(2))), 0.55/(2*TMath::Sqrt(2*TMath::Log(2))));
	  fullFunction_de->SetParameter(5, 0.1);
	  //fullFunction_de->SetParLimits(5, 0.1/(2*TMath::Sqrt(2*TMath::Log(2))), 0.55/(2*TMath::Sqrt(2*TMath::Log(2))));
	  fullFunction_de->SetParameter(8, 0.1);
	  //fullFunction_de->SetParLimits(8, 0.1/(2*TMath::Sqrt(2*TMath::Log(2))), 0.55/(2*TMath::Sqrt(2*TMath::Log(2))));
	  
	  fullFunction_de->SetParameter(1, -1.8);    // Means
	  fullFunction_de->SetParameter(4, 0.0);
	  fullFunction_de->SetParLimits(4, -0.1, 0.1);
	  fullFunction_de->SetParameter(7, 0.628);
	}

      else if (species == "tr" && pBin == 11)
	{
	  fullFunction_tr->SetParameter(0, 2e+6); // Amps
	  fullFunction_tr->SetParameter(3, 4e+5);
	  fullFunction_tr->SetParameter(6, 2.4e+05);

	  fullFunction_tr->SetParameter(2, 0.1);   // Std devs
	  fullFunction_tr->SetParameter(5, 0.1);
	  //fullFunction_tr->SetParLimits(2, 0.1/(2*TMath::Sqrt(2*TMath::Log(2))), 0.55/(2*TMath::Sqrt(2*TMath::Log(2))));
	  fullFunction_tr->SetParameter(8, 0.1);
	  
	  fullFunction_tr->SetParameter(1, -0.6);    // Means
	  fullFunction_tr->SetParameter(4, 0.0);
	  fullFunction_tr->SetParLimits(4, -0.1, 0.1);
	  fullFunction_tr->SetParameter(7, 0.4);
	}
      else if (species == "tr" && pBin == 16)
	{
	  fullFunction_tr->SetParameter(0, 1.8e+6); // Amps
	  fullFunction_tr->SetParameter(3, 2.4e+5);
	  fullFunction_tr->SetParameter(6, 2e+04);

	  fullFunction_tr->SetParameter(2, 0.1);   // Std devs
	  fullFunction_tr->SetParameter(5, 0.1);
	  //fullFunction_tr->SetParLimits(2, 0.1/(2*TMath::Sqrt(2*TMath::Log(2))), 0.55/(2*TMath::Sqrt(2*TMath::Log(2))));
	  fullFunction_tr->SetParameter(8, 0.1);
	  fullFunction_tr->SetParLimits(8, 0.05, 0.12);
	  
	  fullFunction_tr->SetParameter(1, -0.48);    // Means
	  fullFunction_tr->SetParameter(4, 0.0);
	  fullFunction_tr->SetParLimits(4, -0.1, 0.1);
	  fullFunction_tr->SetParameter(7, 0.38);
	  fullFunction_tr->SetParLimits(7, 0.3, 0.42);
	}
      else if (species == "tr" && pBin == 19)
	{
	  delete fullFunction_tr;
	  fullFunction_tr = new TF1("fullFunction_tr", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2) + [6]*exp(-0.5*((x-[7])/[8])**2) + [9]*exp(-0.5*((x-[10])/[11])**2)", -0.7, 0.7);
	  
	  fullFunction_tr->SetParameter(0, 1e+07); // Amps
	  fullFunction_tr->SetParameter(3, 1.5e+06);
	  fullFunction_tr->SetParameter(6, 2e+05);
	  fullFunction_tr->SetParameter(9, 1e+03);

	  fullFunction_tr->SetParameter(2, 0.1);   // Std devs
	  fullFunction_tr->SetParameter(5, 0.1);
	  fullFunction_tr->SetParameter(8, 0.1);
	  fullFunction_tr->SetParameter(11, 0.1);
	  fullFunction_tr->SetParLimits(11, 0.05, 0.12);
	  
	  fullFunction_tr->SetParameter(1, -0.8);    // Means
	  fullFunction_tr->SetParameter(4, -0.4);
	  fullFunction_tr->SetParameter(7, 0.0);
	  fullFunction_tr->SetParLimits(7, -0.1, 0.1);
	  fullFunction_tr->SetParameter(10, 0.33);
	  fullFunction_tr->SetParLimits(10, 0.3, 0.42);
	}
      else if (species == "tr" && pBin == 34)
	{
	  fullFunction_tr->ReleaseParameter(10);	  
	}

      
      
      if (species == "de")
	{
	  h_nSig->Fit(fullFunction_de);
		
	  bkg1->SetParameter(0, fullFunction_de->GetParameter(0));
	  bkg1->SetParameter(1, fullFunction_de->GetParameter(1));
	  bkg1->SetParameter(2, fullFunction_de->GetParameter(2));

	  signal->SetParameter(0, fullFunction_de->GetParameter(3));
	  signal->SetParameter(1, fullFunction_de->GetParameter(4));
	  signal->SetParameter(2, fullFunction_de->GetParameter(5));

	  bkg2->SetParameter(0, fullFunction_de->GetParameter(6));
	  bkg2->SetParameter(1, fullFunction_de->GetParameter(7));
	  bkg2->SetParameter(2, fullFunction_de->GetParameter(8));

	  double sig_bkg1_cross = getIntersection(bkg1->GetParameter(0),bkg1->GetParameter(1),bkg1->GetParameter(2),
						  signal->GetParameter(0),signal->GetParameter(1),signal->GetParameter(2));

	  double sig_bkg2_cross = getIntersection(bkg2->GetParameter(0),bkg2->GetParameter(1),bkg2->GetParameter(2),
						  signal->GetParameter(0),signal->GetParameter(1),signal->GetParameter(2));

	  lowBound = sig_bkg1_cross;
	  highBound = sig_bkg2_cross;

	  // Extra cut tweaks to increase purity
	  if (pBin >= 20 && pBin < 26) lowBound += 0.05;
	  else if (pBin >= 26)
	    {
	      lowBound += 0.05;
	      highBound -= 0.03;
	    }
	  
	  lineLow = new TLine(lowBound, 0.0, lowBound, h_nSig->GetMaximum());    
	  lineHigh = new TLine(highBound, 0.0, highBound, h_nSig->GetMaximum());
	  
	  signalIntegral = signal->Integral(lowBound, highBound);
	  fullIntegral = fullFunction_de->Integral(lowBound, highBound);
	  purity = signalIntegral / fullIntegral;
	  purityStr.Form("%1.3f", purity);
	  text->AddText("Purity = " + purityStr);
      

	  h_nSig->Draw();
	  bkg1->DrawCopy("SAME");
	  bkg2->DrawCopy("SAME");
	  signal->DrawCopy("SAME");
	  text->Draw("SAME");
	  lineLow->Draw("SAME");
	  lineHigh->Draw("SAME");

	  binsUsed.push_back(pBin);
	  lowBounds.push_back(lowBound);
	  highBounds.push_back(highBound);
	  purities.push_back(purity);
	}
      else if (species == "tr" && pBin < 19)
	{
	  h_nSig->Fit(fullFunction_tr);

	  bkg1->SetParameter(0, fullFunction_tr->GetParameter(0));
	  bkg1->SetParameter(1, fullFunction_tr->GetParameter(1));
	  bkg1->SetParameter(2, fullFunction_tr->GetParameter(2));
	  
	  signal->SetParameter(0, fullFunction_tr->GetParameter(3));
	  signal->SetParameter(1, fullFunction_tr->GetParameter(4));
	  signal->SetParameter(2, fullFunction_tr->GetParameter(5));

	  bkg2->SetParameter(0, fullFunction_tr->GetParameter(6));
	  bkg2->SetParameter(1, fullFunction_tr->GetParameter(7));
	  bkg2->SetParameter(2, fullFunction_tr->GetParameter(8));

	  double sig_bkg1_cross = getIntersection(bkg1->GetParameter(0),bkg1->GetParameter(1),bkg1->GetParameter(2),
						  signal->GetParameter(0),signal->GetParameter(1),signal->GetParameter(2));

	  double sig_bkg2_cross = getIntersection(bkg2->GetParameter(0),bkg2->GetParameter(1),bkg2->GetParameter(2),
						  signal->GetParameter(0),signal->GetParameter(1),signal->GetParameter(2));

	  lowBound = sig_bkg1_cross;
	  highBound = sig_bkg2_cross;

	  lineLow = new TLine(lowBound, 0.0, lowBound, h_nSig->GetMaximum());    
	  lineHigh = new TLine(highBound, 0.0, highBound, h_nSig->GetMaximum());

	  signalIntegral = signal->Integral(lowBound, highBound);
	  fullIntegral = fullFunction_tr->Integral(lowBound, highBound);
	  purity = signalIntegral / fullIntegral;
	  purityStr.Form("%1.3f", purity);
	  text->AddText("Purity = " + purityStr);
	  
	  h_nSig->Draw();
	  signal->DrawCopy("SAME");
	  bkg1->DrawCopy("SAME");
	  bkg2->DrawCopy("SAME");
	  text->Draw("SAME");
	  lineLow->Draw("SAME");
	  lineHigh->Draw("SAME");

	  binsUsed.push_back(pBin);
	  lowBounds.push_back(lowBound);
	  highBounds.push_back(highBound);
	  purities.push_back(purity);
	}
      else if (species == "tr" && pBin >= 19)
	{
	  h_nSig->Fit(fullFunction_tr);

	  bkg0->SetParameter(0, fullFunction_tr->GetParameter(0));
	  bkg0->SetParameter(1, fullFunction_tr->GetParameter(1));
	  bkg0->SetParameter(2, fullFunction_tr->GetParameter(2));
	  bkg0->SetParameter(3, fullFunction_tr->GetParameter(3));
	  bkg0->SetParameter(4, fullFunction_tr->GetParameter(4));
	  bkg0->SetParameter(5, fullFunction_tr->GetParameter(5));
	  
	  signal->SetParameter(0, fullFunction_tr->GetParameter(6));
	  signal->SetParameter(1, fullFunction_tr->GetParameter(7));
	  signal->SetParameter(2, fullFunction_tr->GetParameter(8));

	  bkg2->SetParameter(0, fullFunction_tr->GetParameter(9));
	  bkg2->SetParameter(1, fullFunction_tr->GetParameter(10));
	  bkg2->SetParameter(2, fullFunction_tr->GetParameter(11));

	  double sig_bkg1_cross = getIntersection(bkg0->GetParameter(3),bkg0->GetParameter(4),bkg0->GetParameter(5),
	  					  signal->GetParameter(0),signal->GetParameter(1),signal->GetParameter(2));

	  double sig_bkg2_cross = getIntersection(bkg2->GetParameter(0),bkg2->GetParameter(1),bkg2->GetParameter(2),
						  signal->GetParameter(0),signal->GetParameter(1),signal->GetParameter(2));

	  lowBound = sig_bkg1_cross;
	  highBound = sig_bkg2_cross;

	  // Extra cut tweaks to increase purity
	  if (pBin >= 21 && pBin < 24) lowBound += 0.05;
	  else if (pBin >= 24 && pBin < 34)
	    {
	      lowBound += 0.05;
	      highBound -= 0.03;
	    }
	  else if (pBin == 34)
	    {
	      lowBound += 0.08;
	      highBound -= 0.05;
	    }
	  else if (pBin == 35)
	    {
	      lowBound += 0.1;
	      highBound -= 0.03;
	    }
	  else if (pBin == 36)
	    {
	      lowBound +=0.15;
	      highBound -= 0.05;
	    }
	  else if (pBin == 37)
	    {
	      lowBound += 0.09;
	      highBound -= 0.05;
	    }
	  else if (pBin > 37)
	    {
	      lowBound += 0.08;
	      highBound -= 0.05;
	    }

	  lineLow = new TLine(lowBound, 0.0, lowBound, h_nSig->GetMaximum());    
	  lineHigh = new TLine(highBound, 0.0, highBound, h_nSig->GetMaximum());

	  signalIntegral = signal->Integral(lowBound, highBound);
	  fullIntegral = fullFunction_tr->Integral(lowBound, highBound);
	  purity = signalIntegral / fullIntegral;
	  purityStr.Form("%1.3f", purity);
	  text->AddText("Purity = " + purityStr);

	  h_nSig->Draw();
	  signal->DrawCopy("SAME");
	  bkg0->DrawCopy("SAME");
	  bkg2->DrawCopy("SAME");
	  text->Draw("SAME");
	  lineLow->Draw("SAME");
	  lineHigh->Draw("SAME");

	  binsUsed.push_back(pBin);
	  lowBounds.push_back(lowBound);
	  highBounds.push_back(highBound);
	  purities.push_back(purity);
	}
      
      
      nPlots++;
      if (nPlots % 18 == 0)
	{
	  canvas->SaveAs("allPlots_"+species+".pdf"); // save page
	  iPad = 0; // reset to top left pad
	  canvas->Clear(); // next page is copy of previous, so clear/reset it first
	  canvas->Divide(3, 6);
	}
    }

  canvas->SaveAs("allPlots_"+species+".pdf");
  canvas->SaveAs("allPlots_"+species+".pdf]");
  delete canvas;

  // Show where the ideal cuts are
  for (int i = 0; i < binsUsed.size(); i++)
    {
      int bin = binsUsed.at(i);
      
      std::cout << h2->GetXaxis()->GetBinLowEdge(bin) << " <= p < " << h2->GetXaxis()->GetBinLowEdge(bin) + h2->GetXaxis()->GetBinWidth(bin) << ": "
		<< "PURITY = " << purities.at(i) << std::endl
		<< "   " << lowBounds.at(i) << " < z < " << highBounds.at(i)
		<< std::endl
		<< std::endl;
    }
}


void testPurity(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TH2D *h2_nSigp_vs_mom = (TH2D*)file->Get("h2_nSigp_vs_mom");
  TH2D *h2_zd_vs_mom = (TH2D*)file->Get("h2_zd_vs_mom");
  TH2D *h2_zt_vs_mom = (TH2D*)file->Get("h2_zt_vs_mom");

  gStyle->SetOptStat(0);
  
  //printHistos(h2_nSigp_vs_mom, "pr");
  printHistos(h2_zd_vs_mom, "de");
  printHistos(h2_zt_vs_mom, "tr");
  
  file->Close();
}
