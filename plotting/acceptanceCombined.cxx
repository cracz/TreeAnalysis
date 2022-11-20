void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);
 
void acceptanceCombined(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";
  //TString fileName = "Normal_nTracksCorrect_nHits15_pdtEfficiency.picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TProfile2D *p2_pp_vs_eta = (TProfile2D*)file->Get("p2_pp_vs_eta");
  TH2D *h2_phi_vs_eta_EPD = (TH2D*)file->Get("h2_phi_vs_eta_EPD");
  
  TH2D *h2_pT_vs_yCM_pp = (TH2D*)file->Get("h2_pT_vs_yCM_pp");
  TH2D *h2_pT_vs_yCM_pm = (TH2D*)file->Get("h2_pT_vs_yCM_pm");
  TH2D *h2_pT_vs_yCM_kp = (TH2D*)file->Get("h2_pT_vs_yCM_kp");
  TH2D *h2_pT_vs_yCM_km = (TH2D*)file->Get("h2_pT_vs_yCM_km");
  TH2D *h2_pT_vs_yCM_pr = (TH2D*)file->Get("h2_pT_vs_yCM_pr");
  //TH2D *h2_pT_vs_yCM_pr_alt = (TH2D*)file->Get("h2_pT_vs_yCM_pr_alt");
  h2_pT_vs_yCM_pp->SetTitle("");
  h2_pT_vs_yCM_pm->SetTitle("");
  h2_pT_vs_yCM_kp->SetTitle("");
  h2_pT_vs_yCM_km->SetTitle("");
  h2_pT_vs_yCM_pr->SetTitle("");

  Double_t maxScale = h2_pT_vs_yCM_pr->GetMaximum();
  h2_pT_vs_yCM_pp->SetMaximum(maxScale);
  h2_pT_vs_yCM_pm->SetMaximum(maxScale);
  h2_pT_vs_yCM_kp->SetMaximum(maxScale);
  h2_pT_vs_yCM_km->SetMaximum(maxScale);



  TPaveText *text_pp = new TPaveText(-0.9, 0.2, -0.6, 0.5);
  text_pp->SetFillColorAlpha(0, 0);
  text_pp->AddText("#pi^{+}");

  TPaveText *text_pm = new TPaveText(-0.9, 0.2, -0.6, 0.5);
  text_pm->SetFillColorAlpha(0, 0);
  text_pm->AddText("#pi^{-}");

  TPaveText *text_kp = new TPaveText(-0.9, 0.2, -0.6, 0.5);
  text_kp->SetFillColorAlpha(0, 0);
  text_kp->AddText("K^{+}");

  TPaveText *text_km = new TPaveText(-0.9, 0.2, -0.6, 0.5);
  text_km->SetFillColorAlpha(0, 0);
  text_km->AddText("K^{-}");

  TPaveText *text_pr = new TPaveText(-0.9, 0.2, -0.6, 0.5);
  text_pr->SetFillColorAlpha(0, 0);
  text_pr->AddText("p");

  TPaveText *text_extra = new TPaveText(-1.0, 1.5, 1.0, 2.0);
  text_extra->SetFillColorAlpha(0,0);
  text_extra->AddText("#sqrt{s_{NN}} = 3.0 GeV FXT Au+Au");
  text_extra->AddText("Collisions at RHIC");


  TLine *y_mid = new TLine(0, 0, 0, 2.5);
  y_mid->SetLineColor(kRed);
  y_mid->SetLineWidth(2);

  TLine *y_target = new TLine(1.05, 0, 1.05, 2.5);
  y_target->SetLineStyle(9);
  y_target->SetLineColor(kRed);
  y_target->SetLineWidth(2);


  Double_t yCM_low_pp  = 0.0;
  Double_t yCM_high_pp = 0.5;
  Double_t pT_low_pp   = 0.18;
  Double_t pT_high_pp  = 1.6;

  Double_t yCM_low_pm  = 0.0;
  Double_t yCM_high_pm = 0.5;
  Double_t pT_low_pm   = 0.18;
  Double_t pT_high_pm  = 1.6;

  Double_t yCM_low_kp  = 0.0;
  Double_t yCM_high_kp = 0.5;
  Double_t pT_low_kp   = 0.4;
  Double_t pT_high_kp  = 1.6;

  Double_t yCM_low_km  = 0.0;
  Double_t yCM_high_km = 0.5;
  Double_t pT_low_km   = 0.4;
  Double_t pT_high_km  = 1.6;

  Double_t yCM_low_pr  = 0.0;
  Double_t yCM_high_pr = 0.5;
  Double_t pT_low_pr   = 0.4;
  Double_t pT_high_pr  = 2.0;

  TLine *left_pp = new TLine(yCM_low_pp, pT_low_pp, yCM_low_pp, pT_high_pp);
  TLine *right_pp = new TLine(yCM_high_pp, pT_low_pp, yCM_high_pp, pT_high_pp);
  TLine *top_pp = new TLine(yCM_low_pp, pT_high_pp, yCM_high_pp, pT_high_pp);
  TLine *bottom_pp = new TLine(yCM_low_pp, pT_low_pp, yCM_high_pp, pT_low_pp);
  left_pp->SetLineWidth(2);
  right_pp->SetLineWidth(2);
  top_pp->SetLineWidth(2);
  bottom_pp->SetLineWidth(2);
  
  TLine *left_pm = new TLine(yCM_low_pm, pT_low_pm, yCM_low_pm, pT_high_pm);
  TLine *right_pm = new TLine(yCM_high_pm, pT_low_pm, yCM_high_pm, pT_high_pm);
  TLine *top_pm = new TLine(yCM_low_pm, pT_high_pm, yCM_high_pm, pT_high_pm);
  TLine *bottom_pm = new TLine(yCM_low_pm, pT_low_pm, yCM_high_pm, pT_low_pm);
  left_pm->SetLineWidth(2);
  right_pm->SetLineWidth(2);
  top_pm->SetLineWidth(2);
  bottom_pm->SetLineWidth(2);

  TLine *left_kp = new TLine(yCM_low_kp, pT_low_kp, yCM_low_kp, pT_high_kp);
  TLine *right_kp = new TLine(yCM_high_kp, pT_low_kp, yCM_high_kp, pT_high_kp);
  TLine *top_kp = new TLine(yCM_low_kp, pT_high_kp, yCM_high_kp, pT_high_kp);
  TLine *bottom_kp = new TLine(yCM_low_kp, pT_low_kp, yCM_high_kp, pT_low_kp);
  left_kp->SetLineWidth(2);
  right_kp->SetLineWidth(2);
  top_kp->SetLineWidth(2);
  bottom_kp->SetLineWidth(2);

  TLine *left_km = new TLine(yCM_low_km, pT_low_km, yCM_low_km, pT_high_km);
  TLine *right_km = new TLine(yCM_high_km, pT_low_km, yCM_high_km, pT_high_km);
  TLine *top_km = new TLine(yCM_low_km, pT_high_km, yCM_high_km, pT_high_km);
  TLine *bottom_km = new TLine(yCM_low_km, pT_low_km, yCM_high_km, pT_low_km);
  left_km->SetLineWidth(2);
  right_km->SetLineWidth(2);
  top_km->SetLineWidth(2);
  bottom_km->SetLineWidth(2);

  TLine *left_pr = new TLine(yCM_low_pr, pT_low_pr, yCM_low_pr, pT_high_pr);
  TLine *right_pr = new TLine(yCM_high_pr, pT_low_pr, yCM_high_pr, pT_high_pr);
  TLine *top_pr = new TLine(yCM_low_pr, pT_high_pr, yCM_high_pr, pT_high_pr);
  TLine *bottom_pr = new TLine(yCM_low_pr, pT_low_pr, yCM_high_pr, pT_low_pr);
  left_pr->SetLineWidth(2);
  right_pr->SetLineWidth(2);
  top_pr->SetLineWidth(2);
  bottom_pr->SetLineWidth(2);

  TLine *left_pr_sym = new TLine(-0.5,1.0,-0.5,2.5);
  TLine *right_pr_sym = new TLine(0.5,1.0,0.5,2.5);
  TLine *top_pr_sym = new TLine(-0.5,2.5,0.5,2.5);
  TLine *bottom_pr_sym = new TLine(-0.5,1.0,0.5,1.0);
  left_pr_sym->SetLineWidth(2);
  right_pr_sym->SetLineWidth(2);
  top_pr_sym->SetLineWidth(2);
  bottom_pr_sym->SetLineWidth(2);
  //left_pr_sym->SetLineColor(4);
  //right_pr_sym->SetLineColor(4);
  //top_pr_sym->SetLineColor(4);
  //bottom_pr_sym->SetLineColor(4);
  left_pr_sym->SetLineStyle(2);
  right_pr_sym->SetLineStyle(2);
  top_pr_sym->SetLineStyle(2);
  bottom_pr_sym->SetLineStyle(2);


  
 
  gStyle->SetOptStat(0);
 
  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1024,640);
  C->SetFillStyle(4000);
 
  // Number of PADS
  const Int_t Nx = 3;
  const Int_t Ny = 2;
 
  // Margins
  Float_t lMargin = 0.08;
  Float_t rMargin = 0.08;
  Float_t bMargin = 0.11;
  Float_t tMargin = 0.05;
 
  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
   
  /*
  // Dummy histogram.
  TH1F *h = (TH1F*) gROOT->FindObject("histo");
  if (h) delete h;
  h = new TH1F("histo","",100,-5.0,5.0);
  h->FillRandom("gaus",10000);
  h->GetXaxis()->SetTitle("x axis");
  h->GetYaxis()->SetTitle("y axis");
  */

  TH2D *h;
  TPad *pad[Nx][Ny];
 
  for (Int_t i=0;i<Nx;i++) {
    for (Int_t j=0;j<Ny;j++) {

      char hname[16];
      sprintf(hname,"h_%i_%i",i,j);
      
      if (i == 0 && j == 0) h = (TH2D*)h2_pT_vs_yCM_pm->Clone(hname);
      else if (i == 0 && j == 1) h = (TH2D*)h2_pT_vs_yCM_pp->Clone(hname);
      else if (i == 1 && j == 0) h = (TH2D*)h2_pT_vs_yCM_km->Clone(hname);
      else if (i == 1 && j == 1) h = (TH2D*)h2_pT_vs_yCM_kp->Clone(hname);
      else if (i == 2 && j == 0) h = (TH2D*)h2_pT_vs_yCM_pr->Clone(hname);
      //else if (i == 2 && j == 1) continue;

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
      if (i == 2 && j == 1)
	{
	  hFrame->GetXaxis()->SetTickLength(0.0);
	  hFrame->GetYaxis()->SetTickLength(0.0);
	}
      hFrame->Draw();
      
      // y axis range
      //h->GetYaxis()->SetRangeUser(0.0001,1.2*h->GetMaximum());
 
      // Format for y axis
      h->GetYaxis()->SetLabelFont(43);
      h->GetYaxis()->SetLabelSize(20);
      h->GetYaxis()->SetLabelOffset(0.02);
      h->GetYaxis()->SetTitleFont(43);
      h->GetYaxis()->SetTitleSize(29);
      h->GetYaxis()->SetTitleOffset(1.9);
 
      h->GetYaxis()->CenterTitle();
      //h->GetYaxis()->SetNdivisions(505);
 
      // TICKS Y Axis
      //h->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
 
      // Format for x axis
      h->GetXaxis()->SetLabelFont(43);
      h->GetXaxis()->SetLabelSize(20);
      h->GetXaxis()->SetLabelOffset(0.02);
      h->GetXaxis()->SetTitleFont(43);
      h->GetXaxis()->SetTitleSize(29);
      h->GetXaxis()->SetTitleOffset(1.9);
      h->GetXaxis()->CenterTitle();
      h->GetXaxis()->SetNdivisions(505);
 
      // TICKS X Axis
      //h->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);

      // TICKS Z Axis
      h->GetZaxis()->SetTickLength(0.07);
      h->GetZaxis()->SetLabelSize(0.08);

      if (i == 2 && j == 0)
	{
	  h->Draw("colz");
	  y_mid->Draw("same");
	  y_target->Draw("same");
	  gPad->Update();
	  TPaletteAxis *palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
	  palette->SetX2NDC(0.86);
	  palette->Draw();
	}
      else if ( !(i == 2 && j == 1) )
	{
	  h->Draw("col");
	  y_mid->Draw("same");
	  y_target->Draw("same");
	}

      // Draw text and lines
      if (i == 0 && j == 0)
	{
	  text_pm->Draw("same");
	  left_pm->Draw("SAME");
	  right_pm->Draw("SAME");
	  top_pm->Draw("SAME");
	  bottom_pm->Draw("SAME");
	}
      else if (i == 0 && j == 1)
	{
	  text_pp->Draw("same");
	  left_pp->Draw("SAME");
	  right_pp->Draw("SAME");
	  top_pp->Draw("SAME");
	  bottom_pp->Draw("SAME");
	}
      else if (i == 1 && j == 0)
	{
	  text_km->Draw("same");
	  left_km->Draw("SAME");
	  right_km->Draw("SAME");
	  top_km->Draw("SAME");
	  bottom_km->Draw("SAME");
	}
      else if (i == 1 && j == 1)
	{
	  text_kp->Draw("same");
	  left_kp->Draw("SAME");
	  right_kp->Draw("SAME");
	  top_kp->Draw("SAME");
	  bottom_kp->Draw("SAME");
	}
      else if (i == 2 && j == 0)
	{
	  text_pr->Draw("same");
	  left_pr->Draw("SAME");
	  right_pr->Draw("SAME");
	  top_pr->Draw("SAME");
	  bottom_pr->Draw("SAME");
	  left_pr_sym->Draw("SAME");
	  right_pr_sym->Draw("SAME");
	  top_pr_sym->Draw("SAME");
	  bottom_pr_sym->Draw("SAME");
	}
      else if (i == 2 && j == 1)
	{
	  text_extra->Draw("same");
	}
    }
  }
  C->SaveAs("acceptanceCombined.png");
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
	vmaru = 0.0;
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
