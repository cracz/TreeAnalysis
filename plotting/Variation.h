#ifndef VARIATION_H
#define VARIATION_H

#include <vector>
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"

class Variation
{
 public:
  Variation(TString prefix, TString order_n_str);
  ~Variation();
  TString ID;

  TH1D *h_vn_pp;
  TH1D *h_vn_pm;
  TH1D *h_vn_kp;
  TH1D *h_vn_km;
  TH1D *h_vn_pr;
  TH1D *h_vn_pr_alt;
  TH1D *h_vn_de;
  TH1D *h_vn_tr;
  
  TH1D *h_vn_pp_ext;
  TH1D *h_vn_pm_ext;
  TH1D *h_vn_kp_ext;
  TH1D *h_vn_km_ext;
  TH1D *h_vn_pr_ext;
  
  TH1D *h_vn_pr_for;
  
  TH1D *h_vn_yCM_00to10_pp;
  TH1D *h_vn_yCM_10to40_pp;
  TH1D *h_vn_yCM_40to60_pp;
  TH1D *h_vn_yCM_00to10_pm;
  TH1D *h_vn_yCM_10to40_pm;
  TH1D *h_vn_yCM_40to60_pm;
  TH1D *h_vn_yCM_00to10_kp;
  TH1D *h_vn_yCM_10to40_kp;
  TH1D *h_vn_yCM_40to60_kp;
  TH1D *h_vn_yCM_00to10_km;
  TH1D *h_vn_yCM_10to40_km;
  TH1D *h_vn_yCM_40to60_km;
  TH1D *h_vn_yCM_00to10_pr;
  TH1D *h_vn_yCM_10to40_pr;
  TH1D *h_vn_yCM_40to60_pr;
  TH1D *h_vn_yCM_00to10_pr_symm;
  TH1D *h_vn_yCM_10to40_pr_symm;
  TH1D *h_vn_yCM_40to60_pr_symm;

  TH1D *h_vn_pT_00to10_pp;
  TH1D *h_vn_pT_10to40_pp;
  TH1D *h_vn_pT_40to60_pp;
  TH1D *h_vn_pT_00to10_pm;
  TH1D *h_vn_pT_10to40_pm;
  TH1D *h_vn_pT_40to60_pm;
  TH1D *h_vn_pT_00to10_kp;
  TH1D *h_vn_pT_10to40_kp;
  TH1D *h_vn_pT_40to60_kp;
  TH1D *h_vn_pT_00to10_km;
  TH1D *h_vn_pT_10to40_km;
  TH1D *h_vn_pT_40to60_km;
  TH1D *h_vn_pT_00to10_pr;
  TH1D *h_vn_pT_10to40_pr;
  TH1D *h_vn_pT_40to60_pr;

  TH1D *h_vn_yCM_pr_1;
  TH1D *h_vn_yCM_pr_2;
  TH1D *h_vn_yCM_pr_3;

  
  
 private:
  TString fileName;
  TFile *file;
  void initialize(TString order_n_str);
  void fixAttributes(TString order_n_str);
};

#endif
