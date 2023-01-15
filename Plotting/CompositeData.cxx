#include <iostream>
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "PlotUtils.h"
#include "CompositeData.h"


// Constructor for one variation.
CompositeData::CompositeData(TString prefix, Variation* normalData, Variation* var1Data)
{
  onlyOneVariation = true;
  ID = prefix;
  initialize();
  combine(normalData, var1Data);
  saveDetails(normalData);
}


// Constructor for two related variations.
CompositeData::CompositeData(TString prefix, Variation* normalData, Variation* var1Data, Variation* var2Data)
{
  onlyOneVariation = false;
  ID = prefix;
  initialize();
  combine(normalData, var1Data, var2Data);
  saveDetails(normalData);
}


CompositeData::~CompositeData()
{
  delete barlow_vn_pp;
  delete barlow_vn_pm;
  delete barlow_vn_kp;
  delete barlow_vn_km;
  delete barlow_vn_pr;
  delete barlow_vn_pr_alt;
  delete barlow_vn_de;
  delete barlow_vn_tr;
  
  delete barlow_vn_yCM_00to10_pp;
  delete barlow_vn_yCM_10to40_pp;
  delete barlow_vn_yCM_40to60_pp;
  delete barlow_vn_yCM_00to10_pm;
  delete barlow_vn_yCM_10to40_pm;
  delete barlow_vn_yCM_40to60_pm;
  delete barlow_vn_yCM_00to10_kp;
  delete barlow_vn_yCM_10to40_kp;
  delete barlow_vn_yCM_40to60_kp;
  delete barlow_vn_yCM_00to10_km;
  delete barlow_vn_yCM_10to40_km;
  delete barlow_vn_yCM_40to60_km;
  delete barlow_vn_yCM_00to10_pr;
  delete barlow_vn_yCM_10to40_pr;
  delete barlow_vn_yCM_40to60_pr;
  delete barlow_vn_yCM_00to10_pr_symm;
  delete barlow_vn_yCM_10to40_pr_symm;
  delete barlow_vn_yCM_40to60_pr_symm;

  delete barlow_vn_pT_00to10_pp;
  delete barlow_vn_pT_10to40_pp;
  delete barlow_vn_pT_40to60_pp;
  delete barlow_vn_pT_00to10_pm;
  delete barlow_vn_pT_10to40_pm;
  delete barlow_vn_pT_40to60_pm;
  delete barlow_vn_pT_00to10_kp;
  delete barlow_vn_pT_10to40_kp;
  delete barlow_vn_pT_40to60_kp;
  delete barlow_vn_pT_00to10_km;
  delete barlow_vn_pT_10to40_km;
  delete barlow_vn_pT_40to60_km;
  delete barlow_vn_pT_00to10_pr;
  delete barlow_vn_pT_10to40_pr;
  delete barlow_vn_pT_40to60_pr;
}


void CompositeData::initialize()
{
  // Centrality
  barlow_vn_pp = new TH1D("barlow_vn_pp_"+ID, "pp vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);
  barlow_vn_pm = new TH1D("barlow_vn_pm_"+ID, "pm vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);
  barlow_vn_kp = new TH1D("barlow_vn_kp_"+ID, "kp vs cent;Centrality;#Delta/#sigma_{#Delta}", 6, 0, 60);
  barlow_vn_km = new TH1D("barlow_vn_km_"+ID, "km vs cent;Centrality;#Delta/#sigma_{#Delta}", 6, 0, 60);
  barlow_vn_pr = new TH1D("barlow_vn_pr_"+ID, "pr vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);
  barlow_vn_pr_alt = new TH1D("barlow_vn_pr_alt_"+ID, "Alternate pr vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);
  barlow_vn_de = new TH1D("barlow_vn_de_"+ID, "de vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);
  barlow_vn_tr = new TH1D("barlow_vn_tr_"+ID, "tr vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);

  // Pi+ y
  barlow_vn_yCM_00to10_pp = new TH1D("barlow_vn_yCM_00to10_pp_"+ID, "0-10% pp vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 10, 0, 1);
  barlow_vn_yCM_10to40_pp = new TH1D("barlow_vn_yCM_10to40_pp_"+ID, "10-40% pp vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 10, 0, 1);
  barlow_vn_yCM_40to60_pp = new TH1D("barlow_vn_yCM_40to60_pp_"+ID, "40-60% pp vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 10, 0, 1);

  // Pi- y
  barlow_vn_yCM_00to10_pm = new TH1D("barlow_vn_yCM_00to10_pm_"+ID, "0-10% pm vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 10, 0, 1);
  barlow_vn_yCM_10to40_pm = new TH1D("barlow_vn_yCM_10to40_pm_"+ID, "10-40% pm vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 10, 0, 1);
  barlow_vn_yCM_40to60_pm = new TH1D("barlow_vn_yCM_40to60_pm_"+ID, "40-60% pm vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 10, 0, 1);

  // K+ y
  barlow_vn_yCM_00to10_kp = new TH1D("barlow_vn_yCM_00to10_kp_"+ID, "0-10% kp vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 5, 0, 1);
  barlow_vn_yCM_10to40_kp = new TH1D("barlow_vn_yCM_10to40_kp_"+ID, "10-40% kp vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 5, 0, 1);
  barlow_vn_yCM_40to60_kp = new TH1D("barlow_vn_yCM_40to60_kp_"+ID, "40-60% kp vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 5, 0, 1);

  // K- y
  barlow_vn_yCM_00to10_km = new TH1D("barlow_vn_yCM_00to10_km_"+ID, "0-10% km vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 5, 0, 1);
  barlow_vn_yCM_10to40_km = new TH1D("barlow_vn_yCM_10to40_km_"+ID, "10-40% km vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 5, 0, 1);
  barlow_vn_yCM_40to60_km = new TH1D("barlow_vn_yCM_40to60_km_"+ID, "40-60% km vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 5, 0, 1);

  // Proton y
  barlow_vn_yCM_00to10_pr = new TH1D("barlow_vn_yCM_00to10_pr_"+ID, "0-10% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_10to40_pr = new TH1D("barlow_vn_yCM_10to40_pr_"+ID, "10-40% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_40to60_pr = new TH1D("barlow_vn_yCM_40to60_pr_"+ID, "40-60% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);

  // Proton y symmetric
  barlow_vn_yCM_00to10_pr_symm = new TH1D("barlow_vn_yCM_00to10_pr_symm_"+ID, "0-10% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_10to40_pr_symm = new TH1D("barlow_vn_yCM_10to40_pr_symm_"+ID, "10-40% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_40to60_pr_symm = new TH1D("barlow_vn_yCM_40to60_pr_symm_"+ID, "40-60% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  
  // Pi+ pT
  barlow_vn_pT_00to10_pp = new TH1D("barlow_vn_pT_00to10_pp_"+ID, "0-10% pp vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
  barlow_vn_pT_10to40_pp = new TH1D("barlow_vn_pT_10to40_pp_"+ID, "10-40% pp vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
  barlow_vn_pT_40to60_pp = new TH1D("barlow_vn_pT_40to60_pp_"+ID, "40-60% pp vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);

  // Pi- pT
  barlow_vn_pT_00to10_pm = new TH1D("barlow_vn_pT_00to10_pm_"+ID, "0-10% pm vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
  barlow_vn_pT_10to40_pm = new TH1D("barlow_vn_pT_10to40_pm_"+ID, "10-40% pm vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
  barlow_vn_pT_40to60_pm = new TH1D("barlow_vn_pT_40to60_pm_"+ID, "40-60% pm vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);

  // K+ pT
  barlow_vn_pT_00to10_kp = new TH1D("barlow_vn_pT_00to10_kp_"+ID, "0-10% kp vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 5, 0, 2);
  barlow_vn_pT_10to40_kp = new TH1D("barlow_vn_pT_10to40_kp_"+ID, "10-40% kp vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 5, 0, 2);
  barlow_vn_pT_40to60_kp = new TH1D("barlow_vn_pT_40to60_kp_"+ID, "40-60% kp vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 5, 0, 2);

  // K- pT
  barlow_vn_pT_00to10_km = new TH1D("barlow_vn_pT_00to10_km_"+ID, "0-10% km vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 5, 0, 2);
  barlow_vn_pT_10to40_km = new TH1D("barlow_vn_pT_10to40_km_"+ID, "10-40% km vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 5, 0, 2);
  barlow_vn_pT_40to60_km = new TH1D("barlow_vn_pT_40to60_km_"+ID, "40-60% km vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 5, 0, 2);

  // Proton pT
  barlow_vn_pT_00to10_pr = new TH1D("barlow_vn_pT_00to10_pr_"+ID, "0-10% pr vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
  barlow_vn_pT_10to40_pr = new TH1D("barlow_vn_pT_10to40_pr_"+ID, "10-40% pr vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
  barlow_vn_pT_40to60_pr = new TH1D("barlow_vn_pT_40to60_pr_"+ID, "40-60% pr vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
}


// Save normal and varied flow values for each point to the file supplied.
void CompositeData::addRawValuesToFile(TFile* file, TString histogramName, std::vector<DataPoint> vectorOfPoints)
{
  TH1D* temp;
  TString nameWithBinNo;
  file->cd();

  for (int i = 0; i < vectorOfPoints.size(); i++)
    {
      nameWithBinNo.Form(histogramName+"_bin%d", i+1);
      temp = new TH1D(nameWithBinNo, nameWithBinNo, 500, 1, 1);

      if (onlyOneVariation)
	{
	  temp->Fill(vectorOfPoints.at(i).var1Value);
	  temp->Fill(vectorOfPoints.at(i).normalValue);
	}
      else
	{
	  temp->Fill(vectorOfPoints.at(i).var1Value);
	  temp->Fill(vectorOfPoints.at(i).var2Value);
	  temp->Fill(vectorOfPoints.at(i).normalValue);
	}
      
      temp->Write();
      delete temp;
    }
}


void CompositeData::addBarlowValuesToFile(TFile* file, TH1D* barlowHistogram, std::vector<DataPoint> vectorOfPoints)
{
  for (int i = 1; i <= vectorOfPoints.size(); i++)
    { barlowHistogram->SetBinContent(i, vectorOfPoints.at(i-1).deltaByDeltaError); }

  file->cd();
  barlowHistogram->Write();
}


// Save each flow value along with the same value from one variation to see the effect on every point.
//Also save the values used for the Barlow check.
//Normal Variation is only used to get names of histograms.
void CompositeData::saveDetails(Variation* normalData)
{
  TFile *detailsFile = new TFile("Details_"+ID+".root", "RECREATE");

  // Save the flow values 
  addRawValuesToFile(detailsFile, normalData->h_vn_pp->GetName(), v_vn_pp);
  addRawValuesToFile(detailsFile, normalData->h_vn_pm->GetName(), v_vn_pm);
  addRawValuesToFile(detailsFile, normalData->h_vn_kp->GetName(), v_vn_kp);
  addRawValuesToFile(detailsFile, normalData->h_vn_km->GetName(), v_vn_km);
  addRawValuesToFile(detailsFile, normalData->h_vn_pr->GetName(), v_vn_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_pr_alt->GetName(), v_vn_pr_alt);
  addRawValuesToFile(detailsFile, normalData->h_vn_pr->GetName(), v_vn_de);
  addRawValuesToFile(detailsFile, normalData->h_vn_pr->GetName(), v_vn_tr);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to10_pp->GetName(), v_vn_yCM_00to10_pp);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to40_pp->GetName(), v_vn_yCM_10to40_pp);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to60_pp->GetName(), v_vn_yCM_40to60_pp);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to10_pm->GetName(), v_vn_yCM_00to10_pm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to40_pm->GetName(), v_vn_yCM_10to40_pm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to60_pm->GetName(), v_vn_yCM_40to60_pm);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to10_kp->GetName(), v_vn_yCM_00to10_kp);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to40_kp->GetName(), v_vn_yCM_10to40_kp);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to60_kp->GetName(), v_vn_yCM_40to60_kp);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to10_km->GetName(), v_vn_yCM_00to10_km);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to40_km->GetName(), v_vn_yCM_10to40_km);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to60_km->GetName(), v_vn_yCM_40to60_km);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to10_pr->GetName(), v_vn_yCM_00to10_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to40_pr->GetName(), v_vn_yCM_10to40_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to60_pr->GetName(), v_vn_yCM_40to60_pr);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to10_pr_symm->GetName(), v_vn_yCM_00to10_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to40_pr_symm->GetName(), v_vn_yCM_10to40_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to60_pr_symm->GetName(), v_vn_yCM_40to60_pr_symm);

  addRawValuesToFile(detailsFile, normalData->h_vn_pT_00to10_pp->GetName(), v_vn_pT_00to10_pp);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_10to40_pp->GetName(), v_vn_pT_10to40_pp);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_40to60_pp->GetName(), v_vn_pT_40to60_pp);

  addRawValuesToFile(detailsFile, normalData->h_vn_pT_00to10_pm->GetName(), v_vn_pT_00to10_pm);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_10to40_pm->GetName(), v_vn_pT_10to40_pm);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_40to60_pm->GetName(), v_vn_pT_40to60_pm);

  addRawValuesToFile(detailsFile, normalData->h_vn_pT_00to10_kp->GetName(), v_vn_pT_00to10_kp);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_10to40_kp->GetName(), v_vn_pT_10to40_kp);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_40to60_kp->GetName(), v_vn_pT_40to60_kp);

  addRawValuesToFile(detailsFile, normalData->h_vn_pT_00to10_km->GetName(), v_vn_pT_00to10_km);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_10to40_km->GetName(), v_vn_pT_10to40_km);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_40to60_km->GetName(), v_vn_pT_40to60_km);

  addRawValuesToFile(detailsFile, normalData->h_vn_pT_00to10_pr->GetName(), v_vn_pT_00to10_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_10to40_pr->GetName(), v_vn_pT_10to40_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_40to60_pr->GetName(), v_vn_pT_40to60_pr);
  ////

  // Save the significance values
  addBarlowValuesToFile(detailsFile, barlow_vn_pp, v_vn_pp);
  addBarlowValuesToFile(detailsFile, barlow_vn_pm, v_vn_pm);
  addBarlowValuesToFile(detailsFile, barlow_vn_kp, v_vn_kp);
  addBarlowValuesToFile(detailsFile, barlow_vn_km, v_vn_km);
  addBarlowValuesToFile(detailsFile, barlow_vn_pr, v_vn_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_pr_alt, v_vn_pr_alt);
  addBarlowValuesToFile(detailsFile, barlow_vn_de, v_vn_de);
  addBarlowValuesToFile(detailsFile, barlow_vn_tr, v_vn_tr);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to10_pp, v_vn_yCM_00to10_pp);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to40_pp, v_vn_yCM_10to40_pp);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to60_pp, v_vn_yCM_40to60_pp);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to10_pm, v_vn_yCM_00to10_pm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to40_pm, v_vn_yCM_10to40_pm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to60_pm, v_vn_yCM_40to60_pm);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to10_kp, v_vn_yCM_00to10_kp);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to40_kp, v_vn_yCM_10to40_kp);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to60_kp, v_vn_yCM_40to60_kp);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to10_km, v_vn_yCM_00to10_km);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to40_km, v_vn_yCM_10to40_km);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to60_km, v_vn_yCM_40to60_km);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  addBarlowValuesToFile(detailsFile, barlow_vn_pT_00to10_pp, v_vn_pT_00to10_pp);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_10to40_pp, v_vn_pT_10to40_pp);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_40to60_pp, v_vn_pT_40to60_pp);

  addBarlowValuesToFile(detailsFile, barlow_vn_pT_00to10_pm, v_vn_pT_00to10_pm);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_10to40_pm, v_vn_pT_10to40_pm);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_40to60_pm, v_vn_pT_40to60_pm);

  addBarlowValuesToFile(detailsFile, barlow_vn_pT_00to10_kp, v_vn_pT_00to10_kp);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_10to40_kp, v_vn_pT_10to40_kp);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_40to60_kp, v_vn_pT_40to60_kp);

  addBarlowValuesToFile(detailsFile, barlow_vn_pT_00to10_km, v_vn_pT_00to10_km);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_10to40_km, v_vn_pT_10to40_km);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_40to60_km, v_vn_pT_40to60_km);

  addBarlowValuesToFile(detailsFile, barlow_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_40to60_pr, v_vn_pT_40to60_pr);
  ////
  
  detailsFile->Close();
  delete detailsFile;
}


// Calculate necessary info for one histogram type with normal data and one variation
void CompositeData::mergePoints(TH1D* normalHisto, TH1D* var1Histo, std::vector<DataPoint>& vectorOfPoints)
{
  DataPoint point;

  for (int i = 1; i <= normalHisto->GetNbinsX(); i++)
    {
      point.normalValue = normalHisto->GetBinContent(i);
      point.normalError = normalHisto->GetBinError(i);
      
      point.var1Value = var1Histo->GetBinContent(i);
      point.var1Error = var1Histo->GetBinError(i);
      
      point.delta = TMath::Abs(point.var1Value - point.normalValue);
      point.deltaError = TMath::Sqrt(TMath::Abs(TMath::Power(point.var1Error, 2) - TMath::Power(point.normalError, 2)));
      point.deltaByDeltaError = (point.deltaError == 0.0)?0.0:point.delta/point.deltaError;
      point.stdDev = point.delta / TMath::Sqrt(12.0);
      point.variance = TMath::Power(point.delta/TMath::Sqrt(12.0), 2.0);
      vectorOfPoints.push_back(point);
    }
}


// Calculate necessary info for one histogram type with normal data and two variations
void CompositeData::mergePoints(TH1D* normalHisto, TH1D* var1Histo, TH1D* var2Histo, std::vector<DataPoint>& vectorOfPoints)
{
  DataPoint point;

  for (int i = 1; i <= normalHisto->GetNbinsX(); i++)
    {
      point.normalValue = normalHisto->GetBinContent(i);
      point.normalError = normalHisto->GetBinError(i);
      
      point.var1Value = var1Histo->GetBinContent(i);
      point.var1Error = var1Histo->GetBinError(i);

      point.var2Value = var2Histo->GetBinContent(i);
      point.var2Error = var2Histo->GetBinError(i);

      point.getMax();
      point.getMin();
      
      point.delta = TMath::Abs(point.maxValue - point.minValue);
      point.deltaError = TMath::Sqrt(TMath::Abs(TMath::Power(point.maxError, 2) - TMath::Power(point.minError, 2)));
      point.deltaByDeltaError = (point.deltaError == 0.0)?0.0:point.delta/point.deltaError;
      point.stdDev = point.delta / TMath::Sqrt(12.0);
      point.variance = TMath::Power(point.delta/TMath::Sqrt(12.0), 2.0);
      vectorOfPoints.push_back(point);
    }
}


// Combine one variation with the normal data to get the attributes for systematics
void CompositeData::combine(Variation* normalData, Variation* var1Data)
{
  mergePoints(normalData->h_vn_pp, var1Data->h_vn_pp, v_vn_pp);
  mergePoints(normalData->h_vn_pm, var1Data->h_vn_pm, v_vn_pm);
  mergePoints(normalData->h_vn_kp, var1Data->h_vn_kp, v_vn_kp);
  mergePoints(normalData->h_vn_km, var1Data->h_vn_km, v_vn_km);
  mergePoints(normalData->h_vn_pr, var1Data->h_vn_pr, v_vn_pr);
  mergePoints(normalData->h_vn_pr_alt, var1Data->h_vn_pr_alt, v_vn_pr_alt);
  mergePoints(normalData->h_vn_de, var1Data->h_vn_de, v_vn_de);
  mergePoints(normalData->h_vn_tr, var1Data->h_vn_tr, v_vn_tr);

  mergePoints(normalData->h_vn_yCM_00to10_pp, var1Data->h_vn_yCM_00to10_pp, v_vn_yCM_00to10_pp);
  mergePoints(normalData->h_vn_yCM_10to40_pp, var1Data->h_vn_yCM_10to40_pp, v_vn_yCM_10to40_pp);
  mergePoints(normalData->h_vn_yCM_40to60_pp, var1Data->h_vn_yCM_40to60_pp, v_vn_yCM_40to60_pp);

  mergePoints(normalData->h_vn_yCM_00to10_pm, var1Data->h_vn_yCM_00to10_pm, v_vn_yCM_00to10_pm);
  mergePoints(normalData->h_vn_yCM_10to40_pm, var1Data->h_vn_yCM_10to40_pm, v_vn_yCM_10to40_pm);
  mergePoints(normalData->h_vn_yCM_40to60_pm, var1Data->h_vn_yCM_40to60_pm, v_vn_yCM_40to60_pm);

  mergePoints(normalData->h_vn_yCM_00to10_kp, var1Data->h_vn_yCM_00to10_kp, v_vn_yCM_00to10_kp);
  mergePoints(normalData->h_vn_yCM_10to40_kp, var1Data->h_vn_yCM_10to40_kp, v_vn_yCM_10to40_kp);
  mergePoints(normalData->h_vn_yCM_40to60_kp, var1Data->h_vn_yCM_40to60_kp, v_vn_yCM_40to60_kp);

  mergePoints(normalData->h_vn_yCM_00to10_km, var1Data->h_vn_yCM_00to10_km, v_vn_yCM_00to10_km);
  mergePoints(normalData->h_vn_yCM_10to40_km, var1Data->h_vn_yCM_10to40_km, v_vn_yCM_10to40_km);
  mergePoints(normalData->h_vn_yCM_40to60_km, var1Data->h_vn_yCM_40to60_km, v_vn_yCM_40to60_km);

  mergePoints(normalData->h_vn_yCM_00to10_pr, var1Data->h_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  mergePoints(normalData->h_vn_yCM_10to40_pr, var1Data->h_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  mergePoints(normalData->h_vn_yCM_40to60_pr, var1Data->h_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  mergePoints(normalData->h_vn_yCM_00to10_pr_symm, var1Data->h_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to40_pr_symm, var1Data->h_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to60_pr_symm, var1Data->h_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  mergePoints(normalData->h_vn_pT_00to10_pp, var1Data->h_vn_pT_00to10_pp, v_vn_pT_00to10_pp);
  mergePoints(normalData->h_vn_pT_10to40_pp, var1Data->h_vn_pT_10to40_pp, v_vn_pT_10to40_pp);
  mergePoints(normalData->h_vn_pT_40to60_pp, var1Data->h_vn_pT_40to60_pp, v_vn_pT_40to60_pp);

  mergePoints(normalData->h_vn_pT_00to10_pm, var1Data->h_vn_pT_00to10_pm, v_vn_pT_00to10_pm);
  mergePoints(normalData->h_vn_pT_10to40_pm, var1Data->h_vn_pT_10to40_pm, v_vn_pT_10to40_pm);
  mergePoints(normalData->h_vn_pT_40to60_pm, var1Data->h_vn_pT_40to60_pm, v_vn_pT_40to60_pm);

  mergePoints(normalData->h_vn_pT_00to10_kp, var1Data->h_vn_pT_00to10_kp, v_vn_pT_00to10_kp);
  mergePoints(normalData->h_vn_pT_10to40_kp, var1Data->h_vn_pT_10to40_kp, v_vn_pT_10to40_kp);
  mergePoints(normalData->h_vn_pT_40to60_kp, var1Data->h_vn_pT_40to60_kp, v_vn_pT_40to60_kp);

  mergePoints(normalData->h_vn_pT_00to10_km, var1Data->h_vn_pT_00to10_km, v_vn_pT_00to10_km);
  mergePoints(normalData->h_vn_pT_10to40_km, var1Data->h_vn_pT_10to40_km, v_vn_pT_10to40_km);
  mergePoints(normalData->h_vn_pT_40to60_km, var1Data->h_vn_pT_40to60_km, v_vn_pT_40to60_km);

  mergePoints(normalData->h_vn_pT_00to10_pr, var1Data->h_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  mergePoints(normalData->h_vn_pT_10to40_pr, var1Data->h_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  mergePoints(normalData->h_vn_pT_40to60_pr, var1Data->h_vn_pT_40to60_pr, v_vn_pT_40to60_pr);
}



// Combine two variations with the normal data to get the attributes for systematics
void CompositeData::combine(Variation* normalData, Variation* var1Data, Variation* var2Data)
{
  mergePoints(normalData->h_vn_pp, var1Data->h_vn_pp, var2Data->h_vn_pp, v_vn_pp);
  mergePoints(normalData->h_vn_pm, var1Data->h_vn_pm, var2Data->h_vn_pm, v_vn_pm);
  mergePoints(normalData->h_vn_kp, var1Data->h_vn_kp, var2Data->h_vn_kp, v_vn_kp);
  mergePoints(normalData->h_vn_km, var1Data->h_vn_km, var2Data->h_vn_km, v_vn_km);
  mergePoints(normalData->h_vn_pr, var1Data->h_vn_pr, var2Data->h_vn_pr, v_vn_pr);
  mergePoints(normalData->h_vn_pr_alt, var1Data->h_vn_pr_alt, var2Data->h_vn_pr_alt, v_vn_pr_alt);
  mergePoints(normalData->h_vn_de, var1Data->h_vn_de, var2Data->h_vn_de, v_vn_de);
  mergePoints(normalData->h_vn_tr, var1Data->h_vn_tr, var2Data->h_vn_tr, v_vn_tr);

  mergePoints(normalData->h_vn_yCM_00to10_pp, var1Data->h_vn_yCM_00to10_pp, var2Data->h_vn_yCM_00to10_pp, v_vn_yCM_00to10_pp);
  mergePoints(normalData->h_vn_yCM_10to40_pp, var1Data->h_vn_yCM_10to40_pp, var2Data->h_vn_yCM_10to40_pp, v_vn_yCM_10to40_pp);
  mergePoints(normalData->h_vn_yCM_40to60_pp, var1Data->h_vn_yCM_40to60_pp, var2Data->h_vn_yCM_40to60_pp, v_vn_yCM_40to60_pp);

  mergePoints(normalData->h_vn_yCM_00to10_pm, var1Data->h_vn_yCM_00to10_pm, var2Data->h_vn_yCM_00to10_pm, v_vn_yCM_00to10_pm);
  mergePoints(normalData->h_vn_yCM_10to40_pm, var1Data->h_vn_yCM_10to40_pm, var2Data->h_vn_yCM_10to40_pm, v_vn_yCM_10to40_pm);
  mergePoints(normalData->h_vn_yCM_40to60_pm, var1Data->h_vn_yCM_40to60_pm, var2Data->h_vn_yCM_40to60_pm, v_vn_yCM_40to60_pm);

  mergePoints(normalData->h_vn_yCM_00to10_kp, var1Data->h_vn_yCM_00to10_kp, var2Data->h_vn_yCM_00to10_kp, v_vn_yCM_00to10_kp);
  mergePoints(normalData->h_vn_yCM_10to40_kp, var1Data->h_vn_yCM_10to40_kp, var2Data->h_vn_yCM_10to40_kp, v_vn_yCM_10to40_kp);
  mergePoints(normalData->h_vn_yCM_40to60_kp, var1Data->h_vn_yCM_40to60_kp, var2Data->h_vn_yCM_40to60_kp, v_vn_yCM_40to60_kp);

  mergePoints(normalData->h_vn_yCM_00to10_km, var1Data->h_vn_yCM_00to10_km, var2Data->h_vn_yCM_00to10_km, v_vn_yCM_00to10_km);
  mergePoints(normalData->h_vn_yCM_10to40_km, var1Data->h_vn_yCM_10to40_km, var2Data->h_vn_yCM_10to40_km, v_vn_yCM_10to40_km);
  mergePoints(normalData->h_vn_yCM_40to60_km, var1Data->h_vn_yCM_40to60_km, var2Data->h_vn_yCM_40to60_km, v_vn_yCM_40to60_km);

  mergePoints(normalData->h_vn_yCM_00to10_pr, var1Data->h_vn_yCM_00to10_pr, var2Data->h_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  mergePoints(normalData->h_vn_yCM_10to40_pr, var1Data->h_vn_yCM_10to40_pr, var2Data->h_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  mergePoints(normalData->h_vn_yCM_40to60_pr, var1Data->h_vn_yCM_40to60_pr, var2Data->h_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  mergePoints(normalData->h_vn_yCM_00to10_pr_symm, var1Data->h_vn_yCM_00to10_pr_symm, var2Data->h_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to40_pr_symm, var1Data->h_vn_yCM_10to40_pr_symm, var2Data->h_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to60_pr_symm, var1Data->h_vn_yCM_40to60_pr_symm, var2Data->h_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  mergePoints(normalData->h_vn_pT_00to10_pp, var1Data->h_vn_pT_00to10_pp, var2Data->h_vn_pT_00to10_pp, v_vn_pT_00to10_pp);
  mergePoints(normalData->h_vn_pT_10to40_pp, var1Data->h_vn_pT_10to40_pp, var2Data->h_vn_pT_10to40_pp, v_vn_pT_10to40_pp);
  mergePoints(normalData->h_vn_pT_40to60_pp, var1Data->h_vn_pT_40to60_pp, var2Data->h_vn_pT_40to60_pp, v_vn_pT_40to60_pp);

  mergePoints(normalData->h_vn_pT_00to10_pm, var1Data->h_vn_pT_00to10_pm, var2Data->h_vn_pT_00to10_pm, v_vn_pT_00to10_pm);
  mergePoints(normalData->h_vn_pT_10to40_pm, var1Data->h_vn_pT_10to40_pm, var2Data->h_vn_pT_10to40_pm, v_vn_pT_10to40_pm);
  mergePoints(normalData->h_vn_pT_40to60_pm, var1Data->h_vn_pT_40to60_pm, var2Data->h_vn_pT_40to60_pm, v_vn_pT_40to60_pm);

  mergePoints(normalData->h_vn_pT_00to10_kp, var1Data->h_vn_pT_00to10_kp, var2Data->h_vn_pT_00to10_kp, v_vn_pT_00to10_kp);
  mergePoints(normalData->h_vn_pT_10to40_kp, var1Data->h_vn_pT_10to40_kp, var2Data->h_vn_pT_10to40_kp, v_vn_pT_10to40_kp);
  mergePoints(normalData->h_vn_pT_40to60_kp, var1Data->h_vn_pT_40to60_kp, var2Data->h_vn_pT_40to60_kp, v_vn_pT_40to60_kp);

  mergePoints(normalData->h_vn_pT_00to10_km, var1Data->h_vn_pT_00to10_km, var2Data->h_vn_pT_00to10_km, v_vn_pT_00to10_km);
  mergePoints(normalData->h_vn_pT_10to40_km, var1Data->h_vn_pT_10to40_km, var2Data->h_vn_pT_10to40_km, v_vn_pT_10to40_km);
  mergePoints(normalData->h_vn_pT_40to60_km, var1Data->h_vn_pT_40to60_km, var2Data->h_vn_pT_40to60_km, v_vn_pT_40to60_km);

  mergePoints(normalData->h_vn_pT_00to10_pr, var1Data->h_vn_pT_00to10_pr, var2Data->h_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  mergePoints(normalData->h_vn_pT_10to40_pr, var1Data->h_vn_pT_10to40_pr, var2Data->h_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  mergePoints(normalData->h_vn_pT_40to60_pr, var1Data->h_vn_pT_40to60_pr, var2Data->h_vn_pT_40to60_pr, v_vn_pT_40to60_pr);
}
