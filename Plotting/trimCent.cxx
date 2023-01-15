TH1D* trimCent(TH1D *hist)
{
  TH1D *h_trimmed = (TH1D*)hist->Clone();
  Int_t oldBins = hist->GetNbinsX();
  Int_t newBins = (oldBins == 16) ? 12 : 6; // Usually 12, 6 for rebinned kaons.
  h_trimmed->SetBins(newBins, 0, 60);
  h_trimmed->GetXaxis()->SetTitle("Centrality (%)");

  for (int i = 1; i < newBins; i++)
    {
      h_trimmed->SetBinContent(i, hist->GetBinContent(i));
      h_trimmed->SetBinError(i, hist->GetBinError(i));
    }

  return h_trimmed;
}
