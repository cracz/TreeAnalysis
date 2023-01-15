void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);
 
void acceptanceCombinedPDT(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";
  //TString fileName = "Normal_September1.picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TH2D *h2_KToverA_vs_yCM_pr = (TH2D*)file->Get("h2_KToverA_vs_yCM_pr");
  TH2D *h2_KToverA_vs_yCM_de = (TH2D*)file->Get("h2_KToverA_vs_yCM_de");
  TH2D *h2_KToverA_vs_yCM_tr = (TH2D*)file->Get("h2_KToverA_vs_yCM_tr");
  h2_KToverA_vs_yCM_pr->SetTitle("");
  h2_KToverA_vs_yCM_de->SetTitle("");
  h2_KToverA_vs_yCM_tr->SetTitle("");

  Double_t maxScaleKT = h2_KToverA_vs_yCM_pr->GetMaximum();
  h2_KToverA_vs_yCM_de->SetMaximum(maxScaleKT);
  h2_KToverA_vs_yCM_tr->SetMaximum(maxScaleKT);

  TPaveText *text_pr = new TPaveText(-0.6, 0.95, -0.4, 1.35);
  text_pr->SetFillColorAlpha(0, 0);
  text_pr->AddText("p");

  TPaveText *text_de = new TPaveText(-0.9, 0.95, -0.6, 1.45);
  text_de->SetFillColorAlpha(0, 0);
  text_de->AddText("d");

  TPaveText *text_tr = new TPaveText(-0.9, 1.0, -0.6, 1.4);
  text_tr->SetFillColorAlpha(0, 0);
  text_tr->AddText("t");

  TPaveText *text_extra = new TPaveText(-1.0, 2.0, 1.0, 2.5);
  //text_extra->SetFillColorAlpha(0,0);
  text_extra->SetFillColor(0);
  text_extra->SetBorderSize(1);
  text_extra->AddText("#sqrt{s_{NN}} = 3.0 GeV FXT Au+Au");
  text_extra->AddText("Collisions at RHIC");

  TLine *y_mid = new TLine(0, 0, 0, 2.5);
  y_mid->SetLineColor(kRed);
  y_mid->SetLineWidth(2);

  TLine *y_target = new TLine(1.05, 0, 1.05, 2.5);
  y_target->SetLineStyle(9);
  y_target->SetLineColor(kRed);
  y_target->SetLineWidth(2);

  Double_t yCM_low_pr  = 0.0;
  Double_t yCM_high_pr = 0.5;
  Double_t pT_low_pr   = 0.4;
  Double_t pT_high_pr  = 2.0;

  Double_t yCM_low_de  = 0.0;
  Double_t yCM_high_de = 1.0;
  Double_t pT_low_de   = 0.2;
  Double_t pT_high_de  = 1.0;
  Double_t KT_low_de   = 0.04;
  Double_t KT_high_de  = 0.4;

  Double_t yCM_low_tr  = 0.0;
  Double_t yCM_high_tr = 1.0;
  Double_t pT_low_tr   = 0.2;
  Double_t pT_high_tr  = 1.0;
  Double_t KT_low_tr   = 0.04;
  Double_t KT_high_tr  = 0.4;

  TLine *left_pr = new TLine(yCM_low_pr, pT_low_pr, yCM_low_pr, pT_high_pr);
  TLine *right_pr = new TLine(yCM_high_pr, pT_low_pr, yCM_high_pr, pT_high_pr);
  TLine *top_pr = new TLine(yCM_low_pr, pT_high_pr, yCM_high_pr, pT_high_pr);
  TLine *bottom_pr = new TLine(yCM_low_pr, pT_low_pr, yCM_high_pr, pT_low_pr);
  left_pr->SetLineWidth(2);
  right_pr->SetLineWidth(2);
  top_pr->SetLineWidth(2);
  bottom_pr->SetLineWidth(2);

  TLine *left_de = new TLine(yCM_low_de, KT_low_de, yCM_low_de, KT_high_de);
  TLine *right_de = new TLine(yCM_high_de, KT_low_de, yCM_high_de, KT_high_de);
  TLine *top_de = new TLine(yCM_low_de, KT_high_de, yCM_high_de, KT_high_de);
  TLine *bottom_de = new TLine(yCM_low_de, KT_low_de, yCM_high_de, KT_low_de);
  left_de->SetLineWidth(2);
  right_de->SetLineWidth(2);
  top_de->SetLineWidth(2);
  bottom_de->SetLineWidth(2);
  
  TLine *left_tr = new TLine(yCM_low_tr, KT_low_tr, yCM_low_tr, KT_high_tr);
  TLine *right_tr = new TLine(yCM_high_tr, KT_low_tr, yCM_high_tr, KT_high_tr);
  TLine *top_tr = new TLine(yCM_low_tr, KT_high_tr, yCM_high_tr, KT_high_tr);
  TLine *bottom_tr = new TLine(yCM_low_tr, KT_low_tr, yCM_high_tr, KT_low_tr);
  left_tr->SetLineWidth(2);
  right_tr->SetLineWidth(2);
  top_tr->SetLineWidth(2);
  bottom_tr->SetLineWidth(2);

  
  gStyle->SetOptStat(0);
 
  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1024,448);
  C->SetFillStyle(4000);
 
  // Number of PADS
  const Int_t Nx = 3;
  const Int_t Ny = 1;
 
  // Margins
  Float_t lMargin = 0.09;
  Float_t rMargin = 0.09;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.03;
 
  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
  
  TH2D *h;
  TPad *pad[Nx][Ny];
 
  for (Int_t i=0;i<Nx;i++)
    {
      for (Int_t j=0;j<Ny;j++)
	{
	  char hname[16];
	  sprintf(hname,"h_%i_%i",i,j);
      
	  if (i == 0 && j == 0) h = (TH2D*)h2_KToverA_vs_yCM_pr->Clone(hname);
	  else if (i == 1 && j == 0) h = (TH2D*)h2_KToverA_vs_yCM_de->Clone(hname);
	  else if (i == 2 && j == 0) h = (TH2D*)h2_KToverA_vs_yCM_tr->Clone(hname);
	  
	  C->cd(0);
 
	  // Get the pads previously created.
	  char pname[16];
	  sprintf(pname,"pad_%i_%i",i,j);
	  pad[i][j] = (TPad*) gROOT->FindObject(pname);
	  pad[i][j]->Draw();
	  pad[i][j]->SetFillStyle(4000);
	  pad[i][j]->SetFrameFillStyle(4000);
	  pad[i][j]->SetLogz();
	  pad[i][j]->SetTicks();
	  pad[i][j]->cd();
      
	  // Size factors
	  Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
	  Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();

	  TH2D *hFrame = (TH2D*)h->Clone();
	  hFrame->Reset();
	  hFrame->Draw();
      
	  // y axis range
	  //h->GetYaxis()->SetRangeUser(0.0001,1.2*h->GetMaximum());
 
	  // Format for y axis
	  h->GetYaxis()->SetLabelFont(43);
	  h->GetYaxis()->SetLabelSize(20);
	  h->GetYaxis()->SetLabelOffset(0.02);
	  h->GetYaxis()->SetTitleFont(43);
	  h->GetYaxis()->SetTitleSize(29);
	  h->GetYaxis()->SetTitleOffset(1.5);
	  h->GetYaxis()->CenterTitle();
	  h->GetYaxis()->SetNdivisions(505);
 
	  // TICKS Y Axis
	  //h->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
 
	  // Format for x axis
	  h->GetXaxis()->SetLabelFont(43);
	  h->GetXaxis()->SetLabelSize(20);
	  h->GetXaxis()->SetLabelOffset(0.02);
	  h->GetXaxis()->SetTitleFont(43);
	  h->GetXaxis()->SetTitleSize(29);
	  h->GetXaxis()->SetTitleOffset(1.0);
	  h->GetXaxis()->CenterTitle();
	  h->GetXaxis()->SetNdivisions(505);
 
	  // TICKS X Axis
	  //h->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);

	  // TICKS Z Axis
	  h->GetZaxis()->SetTickLength(0.07);
	  h->GetZaxis()->SetLabelSize(0.08);

	  
	  h->SetTitle(";y-y_{mid};(m_{T}-m_{0}) / A (GeV/c^{2})");

	  
	  if (i == 0 && j == 0)
	    {
	      h->Draw("col");
	      y_mid->Draw("same");
	      y_target->Draw("same");
	      left_de->Draw("SAME");
	      right_de->Draw("SAME");
	      top_de->Draw("SAME");
	      bottom_de->Draw("SAME");
	      text_pr->Draw("SAME");
	    }
	  else if (i == 1 && j == 0)
	    {
	      h->Draw("col");
	      y_mid->Draw("same");
	      y_target->Draw("same");
	      left_de->Draw("SAME");
	      right_de->Draw("SAME");
	      top_de->Draw("SAME");
	      bottom_de->Draw("SAME");
	      text_extra->Draw("SAME");
	      text_de->Draw("SAME");
	    }
	  else if (i == 2 && j == 0)
	    {
	      h->Draw("colz");
	      y_mid->Draw("same");
	      y_target->Draw("same");
	      gPad->Update();
	      TPaletteAxis *palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
	      palette->SetX1NDC(0.78);
	      palette->SetX2NDC(0.85);
	      palette->Draw();
	      left_tr->Draw("SAME");
	      right_tr->Draw("SAME");
	      top_tr->Draw("SAME");
	      bottom_tr->Draw("SAME");
	      text_tr->Draw("SAME");
	    }
	}
    }
  C->SaveAs("acceptanceCombinedPDT.png");
  C->cd();
}
 
 
 
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
  if (!C) return;
 
  // Setup Pad layout:
  Float_t vSpacing = 0.0;
  Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
 
  Float_t hSpacing = 0.0;
  Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
 
  Float_t vposd,vposu,vmard,vmaru,vfactor;
  Float_t hposl,hposr,hmarl,hmarr,hfactor;
 
  for (Int_t i=0;i<Nx;i++) {
 
    if (i==0) {
      hposl = 0.0;
      hposr = lMargin + hStep;
      hfactor = hposr-hposl;
      hmarl = lMargin / hfactor;
      hmarr = 0.0;
    } else if (i == Nx-1) {
      hposl = hposr + hSpacing;
      hposr = hposl + hStep + rMargin;
      hfactor = hposr-hposl;
      hmarl = 0.0;
      hmarr = rMargin / (hposr-hposl);
    } else {
      hposl = hposr + hSpacing;
      hposr = hposl + hStep;
      hfactor = hposr-hposl;
      hmarl = 0.0;
      hmarr = 0.0;
    }
 
    for (Int_t j=0;j<Ny;j++) {
 
      if (j==0) {
	vposd = 0.0;
	vposu = bMargin + vStep;
	vfactor = vposu-vposd;
	vmard = bMargin / vfactor;
	//vmaru = 0.0;
	vmaru = 0.05;
      } else if (j == Ny-1) {
	vposd = vposu + vSpacing;
	vposu = vposd + vStep + tMargin;
	vfactor = vposu-vposd;
	vmard = 0.0;
	vmaru = tMargin / (vposu-vposd);
      } else {
	vposd = vposu + vSpacing;
	vposu = vposd + vStep;
	vfactor = vposu-vposd;
	vmard = 0.0;
	vmaru = 0.0;
      }
 
      C->cd(0);
 
      char name[16];
      sprintf(name,"pad_%i_%i",i,j);
      TPad *pad = (TPad*) gROOT->FindObject(name);
      if (pad) delete pad;
      pad = new TPad(name,"",hposl,vposd,hposr,vposu);
      pad->SetLeftMargin(hmarl);
      pad->SetRightMargin(hmarr);
      pad->SetBottomMargin(vmard);
      pad->SetTopMargin(vmaru);
 
      pad->SetFrameBorderMode(0);
      pad->SetBorderMode(0);
      pad->SetBorderSize(0);
 
      pad->Draw();
    }
  }
}
