void vnVsKt(TString jobID, TString order_n_str)
{
  //TH1::SetDefaultSumw2();
  
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 800);
  //canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetLeftMargin(0.14);
  canvas->cd();
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);


  TProfile2D *p2_vn_KT_cent_pp = (TProfile2D*)file->Get("p2_vn_KT_cent_pp");
  TProfile2D *p2_vn_KT_cent_pm = (TProfile2D*)file->Get("p2_vn_KT_cent_pm");
  TProfile2D *p2_vn_KT_cent_kp = (TProfile2D*)file->Get("p2_vn_KT_cent_kp");
  TProfile2D *p2_vn_KT_cent_km = (TProfile2D*)file->Get("p2_vn_KT_cent_km");
  TProfile2D *p2_vn_KT_cent_pr = (TProfile2D*)file->Get("p2_vn_KT_cent_pr");
  TProfile2D *p2_vn_KT_cent_de = (TProfile2D*)file->Get("p2_vn_KT_cent_de");
  TProfile2D *p2_vn_KT_cent_tr = (TProfile2D*)file->Get("p2_vn_KT_cent_tr");

  p2_vn_KT_cent_kp->RebinY();
  p2_vn_KT_cent_km->RebinY();


  TProfile *p_vn_KT_00to10_pp = p2_vn_KT_cent_pp->ProfileY("p_vn_KT_00to10_pp", 15, 16);
  TProfile *p_vn_KT_10to40_pp = p2_vn_KT_cent_pp->ProfileY("p_vn_KT_10to40_pp", 9, 14);
  TProfile *p_vn_KT_40to60_pp = p2_vn_KT_cent_pp->ProfileY("p_vn_KT_40to60_pp", 5, 8);
  
  TProfile *p_vn_KT_00to10_pm = p2_vn_KT_cent_pm->ProfileY("p_vn_KT_00to10_pm", 15, 16);
  TProfile *p_vn_KT_10to40_pm = p2_vn_KT_cent_pm->ProfileY("p_vn_KT_10to40_pm", 9, 14);
  TProfile *p_vn_KT_40to60_pm = p2_vn_KT_cent_pm->ProfileY("p_vn_KT_40to60_pm", 5, 8);

  TProfile *p_vn_KT_00to10_kp = p2_vn_KT_cent_kp->ProfileY("p_vn_KT_00to10_kp", 15, 16);
  TProfile *p_vn_KT_10to40_kp = p2_vn_KT_cent_kp->ProfileY("p_vn_KT_10to40_kp", 9, 14);
  TProfile *p_vn_KT_40to60_kp = p2_vn_KT_cent_kp->ProfileY("p_vn_KT_40to60_kp", 5, 8);

  TProfile *p_vn_KT_00to10_km = p2_vn_KT_cent_km->ProfileY("p_vn_KT_00to10_km", 15, 16);
  TProfile *p_vn_KT_10to40_km = p2_vn_KT_cent_km->ProfileY("p_vn_KT_10to40_km", 9, 14);
  TProfile *p_vn_KT_40to60_km = p2_vn_KT_cent_km->ProfileY("p_vn_KT_40to60_km", 5, 8);

  TProfile *p_vn_KT_00to10_pr = p2_vn_KT_cent_pr->ProfileY("p_vn_KT_00to10_pr", 15, 16);
  TProfile *p_vn_KT_10to40_pr = p2_vn_KT_cent_pr->ProfileY("p_vn_KT_10to40_pr", 9, 14);
  TProfile *p_vn_KT_40to60_pr = p2_vn_KT_cent_pr->ProfileY("p_vn_KT_40to60_pr", 5, 8);

  TProfile *p_vn_KT_00to10_de = p2_vn_KT_cent_de->ProfileY("p_vn_KT_00to10_de", 15, 16);
  TProfile *p_vn_KT_10to40_de = p2_vn_KT_cent_de->ProfileY("p_vn_KT_10to40_de", 9, 14);
  TProfile *p_vn_KT_40to60_de = p2_vn_KT_cent_de->ProfileY("p_vn_KT_40to60_de", 5, 8);

  TProfile *p_vn_KT_00to10_tr = p2_vn_KT_cent_tr->ProfileY("p_vn_KT_00to10_tr", 15, 16);
  TProfile *p_vn_KT_10to40_tr = p2_vn_KT_cent_tr->ProfileY("p_vn_KT_10to40_tr", 9, 14);
  TProfile *p_vn_KT_40to60_tr = p2_vn_KT_cent_tr->ProfileY("p_vn_KT_40to60_tr", 5, 8);


  TH1D *h_vn_KT_00to10_pp = new TH1D("h_vn_KT_00to10_pp", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_10to40_pp = new TH1D("h_vn_KT_10to40_pp", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_40to60_pp = new TH1D("h_vn_KT_40to60_pp", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  TH1D *h_vn_KT_00to10_pm = new TH1D("h_vn_KT_00to10_pm", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_10to40_pm = new TH1D("h_vn_KT_10to40_pm", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_40to60_pm = new TH1D("h_vn_KT_40to60_pm", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  TH1D *h_vn_KT_00to10_kp = new TH1D("h_vn_KT_00to10_kp", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  TH1D *h_vn_KT_10to40_kp = new TH1D("h_vn_KT_10to40_kp", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  TH1D *h_vn_KT_40to60_kp = new TH1D("h_vn_KT_40to60_kp", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 5, 0, 2);

  TH1D *h_vn_KT_00to10_km = new TH1D("h_vn_KT_00to10_km", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  TH1D *h_vn_KT_10to40_km = new TH1D("h_vn_KT_10to40_km", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  TH1D *h_vn_KT_40to60_km = new TH1D("h_vn_KT_40to60_km", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 5, 0, 2);

  TH1D *h_vn_KT_00to10_pr = new TH1D("h_vn_KT_00to10_pr", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_10to40_pr = new TH1D("h_vn_KT_10to40_pr", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_40to60_pr = new TH1D("h_vn_KT_40to60_pr", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  TH1D *h_vn_KT_00to10_de = new TH1D("h_vn_KT_00to10_de", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_10to40_de = new TH1D("h_vn_KT_10to40_de", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_40to60_de = new TH1D("h_vn_KT_40to60_de", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  TH1D *h_vn_KT_00to10_tr = new TH1D("h_vn_KT_00to10_tr", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_10to40_tr = new TH1D("h_vn_KT_10to40_tr", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_KT_40to60_tr = new TH1D("h_vn_KT_40to60_tr", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}", 10, 0, 2);


  h_vn_KT_00to10_kp = p_vn_KT_00to10_kp->ProjectionX();
  h_vn_KT_10to40_kp = p_vn_KT_10to40_kp->ProjectionX();
  h_vn_KT_40to60_kp = p_vn_KT_40to60_kp->ProjectionX();
  h_vn_KT_00to10_km = p_vn_KT_00to10_km->ProjectionX();
  h_vn_KT_10to40_km = p_vn_KT_10to40_km->ProjectionX();
  h_vn_KT_40to60_km = p_vn_KT_40to60_km->ProjectionX();
  
  h_vn_KT_00to10_pp = p_vn_KT_00to10_pp->ProjectionX();
  h_vn_KT_10to40_pp = p_vn_KT_10to40_pp->ProjectionX();
  h_vn_KT_40to60_pp = p_vn_KT_40to60_pp->ProjectionX();

  h_vn_KT_00to10_pm = p_vn_KT_00to10_pm->ProjectionX();
  h_vn_KT_10to40_pm = p_vn_KT_10to40_pm->ProjectionX();
  h_vn_KT_40to60_pm = p_vn_KT_40to60_pm->ProjectionX();

  h_vn_KT_00to10_pr = p_vn_KT_00to10_pr->ProjectionX();
  h_vn_KT_10to40_pr = p_vn_KT_10to40_pr->ProjectionX();
  h_vn_KT_40to60_pr = p_vn_KT_40to60_pr->ProjectionX();

  h_vn_KT_00to10_de = p_vn_KT_00to10_de->ProjectionX();
  h_vn_KT_10to40_de = p_vn_KT_10to40_de->ProjectionX();
  h_vn_KT_40to60_de = p_vn_KT_40to60_de->ProjectionX();

  h_vn_KT_00to10_tr = p_vn_KT_00to10_tr->ProjectionX();
  h_vn_KT_10to40_tr = p_vn_KT_10to40_tr->ProjectionX();
  h_vn_KT_40to60_tr = p_vn_KT_40to60_tr->ProjectionX();


  THStack *ppKtStack   = new THStack("ppKtStack", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}");
  THStack *pmKtStack   = new THStack("pmKtStack", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}");
  THStack *kpKtStack   = new THStack("kpKtStack", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}");
  THStack *kmKtStack   = new THStack("kmKtStack", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}");
  THStack *prKtStack   = new THStack("prKtStack", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}");
  THStack *deKtStack   = new THStack("deKtStack", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}");
  THStack *trKtStack   = new THStack("trKtStack", ";m_{T}-m_{0} (GeV);v_{"+order_n_str+"}");


  h_vn_KT_00to10_pp->SetMarkerStyle(20);
  h_vn_KT_10to40_pp->SetMarkerStyle(20);
  h_vn_KT_40to60_pp->SetMarkerStyle(20);
  h_vn_KT_00to10_pp->SetMarkerColor(2);
  h_vn_KT_10to40_pp->SetMarkerColor(4);
  h_vn_KT_40to60_pp->SetMarkerColor(8);
  h_vn_KT_00to10_pp->SetMarkerSize(2);
  h_vn_KT_10to40_pp->SetMarkerSize(2);
  h_vn_KT_40to60_pp->SetMarkerSize(2);
  h_vn_KT_00to10_pp->SetLineColor(2);
  h_vn_KT_10to40_pp->SetLineColor(4);
  h_vn_KT_40to60_pp->SetLineColor(8);
  
  h_vn_KT_00to10_pm->SetMarkerStyle(20);
  h_vn_KT_10to40_pm->SetMarkerStyle(20);
  h_vn_KT_40to60_pm->SetMarkerStyle(20);
  h_vn_KT_00to10_pm->SetMarkerColor(2);
  h_vn_KT_10to40_pm->SetMarkerColor(4);
  h_vn_KT_40to60_pm->SetMarkerColor(8);
  h_vn_KT_00to10_pm->SetMarkerSize(2);
  h_vn_KT_10to40_pm->SetMarkerSize(2);
  h_vn_KT_40to60_pm->SetMarkerSize(2);
  h_vn_KT_00to10_pm->SetLineColor(2);
  h_vn_KT_10to40_pm->SetLineColor(4);
  h_vn_KT_40to60_pm->SetLineColor(8);
  
  h_vn_KT_00to10_kp->SetMarkerStyle(20);
  h_vn_KT_10to40_kp->SetMarkerStyle(20);
  h_vn_KT_40to60_kp->SetMarkerStyle(20);
  h_vn_KT_00to10_kp->SetMarkerColor(2);
  h_vn_KT_10to40_kp->SetMarkerColor(4);
  h_vn_KT_40to60_kp->SetMarkerColor(8);
  h_vn_KT_00to10_kp->SetMarkerSize(2);
  h_vn_KT_10to40_kp->SetMarkerSize(2);
  h_vn_KT_40to60_kp->SetMarkerSize(2);
  h_vn_KT_00to10_kp->SetLineColor(2);
  h_vn_KT_10to40_kp->SetLineColor(4);
  h_vn_KT_40to60_kp->SetLineColor(8);
  
  h_vn_KT_00to10_km->SetMarkerStyle(20);
  h_vn_KT_10to40_km->SetMarkerStyle(20);
  h_vn_KT_40to60_km->SetMarkerStyle(20);
  h_vn_KT_00to10_km->SetMarkerColor(2);
  h_vn_KT_10to40_km->SetMarkerColor(4);
  h_vn_KT_40to60_km->SetMarkerColor(8);
  h_vn_KT_00to10_km->SetMarkerSize(2);
  h_vn_KT_10to40_km->SetMarkerSize(2);
  h_vn_KT_40to60_km->SetMarkerSize(2);
  h_vn_KT_00to10_km->SetLineColor(2);
  h_vn_KT_10to40_km->SetLineColor(4);
  h_vn_KT_40to60_km->SetLineColor(8);
  
  h_vn_KT_00to10_pr->SetMarkerStyle(20);
  h_vn_KT_10to40_pr->SetMarkerStyle(20);
  h_vn_KT_40to60_pr->SetMarkerStyle(20);
  h_vn_KT_00to10_pr->SetMarkerColor(2);
  h_vn_KT_10to40_pr->SetMarkerColor(4);
  h_vn_KT_40to60_pr->SetMarkerColor(8);
  h_vn_KT_00to10_pr->SetMarkerSize(2);
  h_vn_KT_10to40_pr->SetMarkerSize(2);
  h_vn_KT_40to60_pr->SetMarkerSize(2);
  h_vn_KT_00to10_pr->SetLineColor(2);
  h_vn_KT_10to40_pr->SetLineColor(4);
  h_vn_KT_40to60_pr->SetLineColor(8);
  
  h_vn_KT_00to10_de->SetMarkerStyle(20);
  h_vn_KT_10to40_de->SetMarkerStyle(20);
  h_vn_KT_40to60_de->SetMarkerStyle(20);
  h_vn_KT_00to10_de->SetMarkerColor(2);
  h_vn_KT_10to40_de->SetMarkerColor(4);
  h_vn_KT_40to60_de->SetMarkerColor(8);
  h_vn_KT_00to10_de->SetMarkerSize(2);
  h_vn_KT_10to40_de->SetMarkerSize(2);
  h_vn_KT_40to60_de->SetMarkerSize(2);
  h_vn_KT_00to10_de->SetLineColor(2);
  h_vn_KT_10to40_de->SetLineColor(4);
  h_vn_KT_40to60_de->SetLineColor(8);

  h_vn_KT_00to10_tr->SetMarkerStyle(20);
  h_vn_KT_10to40_tr->SetMarkerStyle(20);
  h_vn_KT_40to60_tr->SetMarkerStyle(20);
  h_vn_KT_00to10_tr->SetMarkerColor(2);
  h_vn_KT_10to40_tr->SetMarkerColor(4);
  h_vn_KT_40to60_tr->SetMarkerColor(8);
  h_vn_KT_00to10_tr->SetMarkerSize(2);
  h_vn_KT_10to40_tr->SetMarkerSize(2);
  h_vn_KT_40to60_tr->SetMarkerSize(2);
  h_vn_KT_00to10_tr->SetLineColor(2);
  h_vn_KT_10to40_tr->SetLineColor(4);
  h_vn_KT_40to60_tr->SetLineColor(8);


  if (order_n_str == "2")
    {
      ppKtStack->Add(h_vn_KT_00to10_pp);
      ppKtStack->Add(h_vn_KT_10to40_pp);
      ppKtStack->Add(h_vn_KT_40to60_pp);

      pmKtStack->Add(h_vn_KT_00to10_pm);
      pmKtStack->Add(h_vn_KT_10to40_pm);
      pmKtStack->Add(h_vn_KT_40to60_pm);

      kpKtStack->Add(h_vn_KT_00to10_kp);
      kpKtStack->Add(h_vn_KT_10to40_kp);
      kpKtStack->Add(h_vn_KT_40to60_kp);

      kmKtStack->Add(h_vn_KT_10to40_km);

      prKtStack->Add(h_vn_KT_00to10_pr);
      prKtStack->Add(h_vn_KT_10to40_pr);
      prKtStack->Add(h_vn_KT_40to60_pr);


      ppKtStack->Draw();
      ppKtStack->GetYaxis()->SetTitleOffset(1.7);
      ppKtStack->GetXaxis()->SetNdivisions(210);
      ppKtStack->Draw();
      //ppKtStack->SetMaximum(0.03);
      //ppKtStack->SetMinimum(-0.05);
      ppKtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //ppLegend->Draw();
      //pKText->Draw();
      canvas->SaveAs(jobID + "_ppKtStack.png");
      canvas->Clear();

      pmKtStack->Draw();
      pmKtStack->GetYaxis()->SetTitleOffset(1.7);
      pmKtStack->GetXaxis()->SetNdivisions(210);
      pmKtStack->Draw();
      //pmKtStack->SetMaximum(0.03);
      //pmKtStack->SetMinimum(-0.04);
      pmKtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //pmLegend->Draw();
      //pmText->Draw();
      canvas->SaveAs(jobID + "_pmKtStack.png");
      canvas->Clear();

      kpKtStack->Draw();
      kpKtStack->GetYaxis()->SetTitleOffset(1.7);
      kpKtStack->GetXaxis()->SetNdivisions(210);
      kpKtStack->Draw();
      //kpKtStack->SetMaximum(0.08);
      //kpKtStack->SetMinimum(-0.13);
      kpKtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //kpLegend->Draw();
      //kKText->Draw();
      canvas->SaveAs(jobID + "_kpKtStack.png");
      canvas->Clear();

      kmKtStack->Draw();
      kmKtStack->GetYaxis()->SetTitleOffset(1.7);
      kmKtStack->GetXaxis()->SetNdivisions(210);
      kmKtStack->Draw();
      //kmKtStack->SetMaximum(0.02);
      //kmKtStack->SetMinimum(-0.06);
      kmKtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //kmLegend->Draw();
      //kmText->Draw();
      canvas->SaveAs(jobID + "_kmKtStack.png");
      canvas->Clear();

      prKtStack->Draw();
      prKtStack->GetYaxis()->SetTitleOffset(1.7);
      prKtStack->GetXaxis()->SetNdivisions(210);
      prKtStack->Draw();
      //prKtStack->SetMaximum(0.1);
      //prKtStack->SetMinimum(-0.08);
      prKtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //prLegend->Draw();
      //prText_y->Draw();
      canvas->SaveAs(jobID + "_prKtStack.png");
      canvas->Clear();
    }
  else if (order_n_str == "3")
    {
      TLine *zeroLine_Kt = new TLine(0, 0, 2, 0);
      zeroLine_Kt->SetLineStyle(9);

      
      ppKtStack->Add(h_vn_KT_00to10_pp);
      ppKtStack->Add(h_vn_KT_10to40_pp);
      ppKtStack->Add(h_vn_KT_40to60_pp);

      pmKtStack->Add(h_vn_KT_00to10_pm);
      pmKtStack->Add(h_vn_KT_10to40_pm);
      pmKtStack->Add(h_vn_KT_40to60_pm);

      kpKtStack->Add(h_vn_KT_00to10_kp);
      kpKtStack->Add(h_vn_KT_10to40_kp);
      kpKtStack->Add(h_vn_KT_40to60_kp);

      kmKtStack->Add(h_vn_KT_10to40_km);

      prKtStack->Add(h_vn_KT_00to10_pr);
      prKtStack->Add(h_vn_KT_10to40_pr);
      prKtStack->Add(h_vn_KT_40to60_pr);

      deKtStack->Add(h_vn_KT_00to10_de);
      deKtStack->Add(h_vn_KT_10to40_de);
      deKtStack->Add(h_vn_KT_40to60_de);

      trKtStack->Add(h_vn_KT_00to10_tr);
      trKtStack->Add(h_vn_KT_10to40_tr);
      trKtStack->Add(h_vn_KT_40to60_tr);


      TLegend *ppLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      ppLegend->AddEntry(h_vn_KT_00to10_pp, "0 - 10%");
      ppLegend->AddEntry(h_vn_KT_10to40_pp, "10 - 40%");
      ppLegend->AddEntry(h_vn_KT_40to60_pp, "40 - 60%");
      ppLegend->SetBorderSize(0);
      ppLegend->SetFillColorAlpha(0,0);

      TLegend *pmLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      pmLegend->AddEntry(h_vn_KT_00to10_pm, "0 - 10%");
      pmLegend->AddEntry(h_vn_KT_10to40_pm, "10 - 40%");
      pmLegend->AddEntry(h_vn_KT_40to60_pm, "40 - 60%");
      pmLegend->SetBorderSize(0);
      pmLegend->SetFillColorAlpha(0,0);

      TLegend *kpLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      kpLegend->AddEntry(h_vn_KT_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(h_vn_KT_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(h_vn_KT_40to60_kp, "40 - 60%");
      kpLegend->SetBorderSize(0);
      kpLegend->SetFillColorAlpha(0,0);

      TLegend *kmLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      //kmLegend->AddEntry(h_vn_KT_00to10_km, "0 - 10%");
      kmLegend->AddEntry(h_vn_KT_10to40_km, "10 - 40%");
      //kmLegend->AddEntry(h_vn_KT_40to60_km, "40 - 60%");
      kmLegend->SetBorderSize(0);
      kmLegend->SetFillColorAlpha(0,0);

      TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      prLegend->AddEntry(h_vn_KT_00to10_pr, "0 - 10%");
      prLegend->AddEntry(h_vn_KT_10to40_pr, "10 - 40%");
      prLegend->AddEntry(h_vn_KT_40to60_pr, "40 - 60%");
      prLegend->SetBorderSize(0);
      prLegend->SetFillColorAlpha(0,0);

      TLegend *deLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      deLegend->AddEntry(h_vn_KT_00to10_de, "0 - 10%");
      deLegend->AddEntry(h_vn_KT_10to40_de, "10 - 40%");
      deLegend->AddEntry(h_vn_KT_40to60_de, "40 - 60%");
      deLegend->SetBorderSize(0);
      deLegend->SetFillColorAlpha(0,0);

      TLegend *trLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      trLegend->AddEntry(h_vn_KT_00to10_tr, "0 - 10%");
      trLegend->AddEntry(h_vn_KT_10to40_tr, "10 - 40%");
      trLegend->AddEntry(h_vn_KT_40to60_tr, "40 - 60%");
      trLegend->SetBorderSize(0);
      trLegend->SetFillColorAlpha(0,0);

      

      TPaveText *ppText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      ppText->AddText("#pi^{+}");
      ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppText->AddText("0 < y_{CM} < 0.5 GeV");
      ppText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      ppText->SetFillColorAlpha(0,0);
      ppText->SetLineColorAlpha(0,0);

      TPaveText *pmText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      pmText->AddText("#pi^{-}");
      pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmText->AddText("0 < y_{CM} < 0.5 GeV");
      pmText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      pmText->SetFillColorAlpha(0,0);
      pmText->SetLineColorAlpha(0,0);

      TPaveText *kpText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      kpText->AddText("K^{+}");
      kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpText->AddText("0 < y_{CM} < 0.5 GeV");
      kpText->AddText("0.4 #leq p_{T} #leq 1.6 GeV");
      kpText->SetFillColorAlpha(0,0);
      kpText->SetLineColorAlpha(0,0);

      TPaveText *kmText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      kmText->AddText("K^{-}");
      kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmText->AddText("0 < y_{CM} < 0.5 GeV");
      kmText->AddText("0.4 #leq p_{T} #leq 1.6 GeV");
      kmText->SetFillColorAlpha(0,0);
      kmText->SetLineColorAlpha(0,0);

      TPaveText *prText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);

      TPaveText *deText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      deText->AddText("Deuteron");
      deText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      deText->AddText("0 < y_{CM} < 0.5 GeV");
      deText->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
      deText->SetFillColorAlpha(0,0);
      deText->SetLineColorAlpha(0,0);

      TPaveText *trText = new TPaveText(0.8, 0.07, 1.8, 0.18, "NB");
      trText->AddText("Triton");
      trText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      trText->AddText("0 < y_{CM} < 0.5 GeV");
      trText->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
      trText->SetFillColorAlpha(0,0);
      trText->SetLineColorAlpha(0,0);


      Double_t ptUpperBound = 0.25;
      Double_t ptLowerBound = -0.25;


      ppKtStack->Draw();
      ppKtStack->GetYaxis()->SetTitleOffset(1.7);
      ppKtStack->GetXaxis()->SetNdivisions(210);
      ppKtStack->Draw();
      ppKtStack->SetMaximum(ptUpperBound);
      ppKtStack->SetMinimum(ptLowerBound);
      ppKtStack->Draw("NOSTACK E1P");
      zeroLine_Kt->Draw("SAME");
      ppKtStack->Draw("NOSTACK E1P SAME");
      ppLegend->Draw();
      ppText->Draw();
      canvas->SaveAs(jobID + "_ppKtStack.png");
      canvas->Clear();

      pmKtStack->Draw();
      pmKtStack->GetYaxis()->SetTitleOffset(1.7);
      pmKtStack->GetXaxis()->SetNdivisions(210);
      pmKtStack->Draw();
      pmKtStack->SetMaximum(ptUpperBound);
      pmKtStack->SetMinimum(ptLowerBound);
      pmKtStack->Draw("NOSTACK E1P");
      zeroLine_Kt->Draw("SAME");
      pmKtStack->Draw("NOSTACK E1P SAME");
      pmLegend->Draw();
      pmText->Draw();
      canvas->SaveAs(jobID + "_pmKtStack.png");
      canvas->Clear();

      kpKtStack->Draw();
      kpKtStack->GetYaxis()->SetTitleOffset(1.7);
      kpKtStack->GetXaxis()->SetNdivisions(210);
      kpKtStack->Draw();
      kpKtStack->SetMaximum(ptUpperBound);
      kpKtStack->SetMinimum(ptLowerBound);
      kpKtStack->Draw("NOSTACK E1P");
      zeroLine_Kt->Draw("SAME");
      kpKtStack->Draw("NOSTACK E1P SAME");
      kpLegend->Draw();
      kpText->Draw();
      canvas->SaveAs(jobID + "_kpKtStack.png");
      canvas->Clear();
      
      kmKtStack->Draw();
      kmKtStack->GetYaxis()->SetTitleOffset(1.7);
      kmKtStack->GetXaxis()->SetNdivisions(210);
      kmKtStack->Draw();
      kmKtStack->SetMaximum(ptUpperBound);
      kmKtStack->SetMinimum(ptLowerBound);
      kmKtStack->Draw("NOSTACK E1P");
      zeroLine_Kt->Draw("SAME");
      kmKtStack->Draw("NOSTACK E1P SAME");
      kmLegend->Draw();
      kmText->Draw();
      canvas->SaveAs(jobID + "_kmKtStack.png");
      canvas->Clear();
      gStyle->SetErrorX(0);
      
      prKtStack->Draw();
      prKtStack->GetYaxis()->SetTitleOffset(1.7);
      prKtStack->GetXaxis()->SetNdivisions(210);
      prKtStack->Draw();
      prKtStack->SetMaximum(ptUpperBound);
      prKtStack->SetMinimum(ptLowerBound);
      prKtStack->Draw("NOSTACK E1P");
      zeroLine_Kt->Draw("SAME");
      prKtStack->Draw("NOSTACK E1P SAME");
      prLegend->Draw();
      prText->Draw();
      canvas->SaveAs(jobID + "_prKtStack.png");
      canvas->Clear();

      deKtStack->Draw();
      deKtStack->GetYaxis()->SetTitleOffset(1.7);
      deKtStack->GetXaxis()->SetNdivisions(210);
      deKtStack->Draw();
      deKtStack->SetMaximum(ptUpperBound);
      deKtStack->SetMinimum(ptLowerBound);
      deKtStack->Draw("NOSTACK E1P");
      zeroLine_Kt->Draw("SAME");
      deKtStack->Draw("NOSTACK E1P SAME");
      deLegend->Draw();
      deText->Draw();
      canvas->SaveAs(jobID + "_deKtStack.png");
      canvas->Clear();

      trKtStack->Draw();
      trKtStack->GetYaxis()->SetTitleOffset(1.7);
      trKtStack->GetXaxis()->SetNdivisions(210);
      trKtStack->Draw();
      trKtStack->SetMaximum(ptUpperBound);
      trKtStack->SetMinimum(ptLowerBound);
      trKtStack->Draw("NOSTACK E1P");
      zeroLine_Kt->Draw("SAME");
      trKtStack->Draw("NOSTACK E1P SAME");
      trLegend->Draw();
      trText->Draw();
      canvas->SaveAs(jobID + "_trKtStack.png");
      canvas->Clear();
    }

  file->Close();
  
}
