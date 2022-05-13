TH1D* flipHist(TH1D *hist)
{
  TH1D *h_flipped = (TH1D*)hist->Clone();//(TString)hist->GetName() + "_flip");
  Int_t bins = hist->GetNbinsX();
  
  Int_t j = 1;
  for (int i = bins; i >= 1; i--)
    {
      h_flipped->SetBinContent(j, hist->GetBinContent(i));
      h_flipped->SetBinError(j, hist->GetBinError(i));
      j++;
    }

  return h_flipped;
}
