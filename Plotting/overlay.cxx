void usage();

void overlay(TString jobID, TString fitType)
{
  if (jobID.EqualTo("")) { usage(); return; }
  if (!fitType.EqualTo("E") && !fitType.EqualTo("L") && !fitType.EqualTo("B")) { usage(); return; }
  
  TString rawFileName = jobID + ".picoDst.result.combined.root";
  TString intFileName = jobID + ".picoDst.result."+fitType+"yields.root";

  TFile *rawFile = TFile::Open(rawFileName, "read");
  if (!rawFile) { std::cout << "File with raw spectra not found!" << std::endl; return; }

  TFile *intFile = TFile::Open(intFileName, "read");
  if (!intFile) { std::cout << "File with integrated spectra not found!" << std::endl; return; }


  TH1D *pp_raw = (TH1D*)rawFile->Get("h_pp_dndy");
  TH1D *pm_raw = (TH1D*)rawFile->Get("h_pm_dndy");
  TH1D *kp_raw = (TH1D*)rawFile->Get("h_kp_dndy");
  TH1D *km_raw = (TH1D*)rawFile->Get("h_km_dndy");
  TH1D *pr_raw = (TH1D*)rawFile->Get("h_pr_dndy");

  TH1D *pp_int = (TH1D*)intFile->Get("pp_yield");
  TH1D *pm_int = (TH1D*)intFile->Get("pm_yield");
  TH1D *kp_int = (TH1D*)intFile->Get("kp_yield");
  TH1D *km_int = (TH1D*)intFile->Get("km_yield");
  TH1D *pr_int = (TH1D*)intFile->Get("pr_yield");

  
  TObjArray *rawArray = new TObjArray();
  rawArray->Add(pp_raw);
  rawArray->Add(pm_raw);
  rawArray->Add(kp_raw);
  rawArray->Add(km_raw);
  rawArray->Add(pr_raw);
  rawArray->Compress();

  TObjArray *intArray = new TObjArray();
  intArray->Add(pp_int);
  intArray->Add(pm_int);
  intArray->Add(kp_int);
  intArray->Add(km_int);
  intArray->Add(pr_int);
  intArray->Compress();


  TCanvas *canvas = new TCanvas("canvas", "Title", 875, 675);
  canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetLeftMargin(0.13);
  gStyle->SetOptStat(0);
  
  TLegend *leg = new TLegend(0.72,0.72,0.9,0.9);


  for (Int_t i = 0; i < rawArray->GetEntriesFast(); i++)
    {
      TH1D *rawYield = (TH1D*)rawArray->At(i);
      TH1D *intYield = (TH1D*)intArray->At(i);

      rawYield->SetMarkerStyle(27);
      rawYield->SetMarkerColor(4);
      rawYield->SetMarkerSize(2.5);
      rawYield->SetLineColor(4);

      intYield->SetMarkerStyle(24);
      intYield->SetMarkerColor(2);
      intYield->SetMarkerSize(1.5);
      intYield->SetLineColor(2);
      intYield->SetTitleOffset(1.4, "y");

      intYield->DrawClone("E1");
      rawYield->DrawClone("E1 SAME");

      intYield->SetMarkerStyle(1);
      rawYield->SetMarkerStyle(1);
      intYield->Draw("E1 SAME");
      rawYield->Draw("E1 SAME");

      leg->Clear();
      leg->AddEntry(rawYield->GetName(), "Raw", "ep");
      leg->AddEntry(intYield->GetName(), "Integrated", "ep");
      leg->Draw();

      canvas->Update();
      canvas->SaveAs(jobID+"."+rawYield->GetName()+"_fit"+fitType+"_overlay.png");
      canvas->Clear();
    }

  canvas->Close();
  rawFile->Close();
  intFile->Close();
  delete rawArray;
  delete intArray;
  delete rawFile;
  delete intFile;
}//End overlay()




void usage()
{
  std::cout << "Usage: root -l -b -q 'overlay.cxx(\"[job ID]\",\"[fit type]\")'" << std::endl
	    << "Fit types: E = exponential, L = Levy, B = Blast-wave" << std::endl;
}

