#include "../include/SetupAttributes.h"


Int_t SetupAttributes::getRunIteration()
{
  return runIteration;
}

Bool_t SetupAttributes::resolutionsWereFound()
{
  return resolutionsFound;
}

Bool_t SetupAttributes::tpcEfficienciesWereFound()
{
  return tpcEfficienciesFound;
}

Bool_t SetupAttributes::tofEfficienciesWereFound()
{
  return tofEfficienciesFound;
}

Bool_t SetupAttributes::v1WeightsWereFound()
{
  return v1WeightsFound;
}

void SetupAttributes::setCorrectionFileAndRunIteration(TString fileName)
{
  if (!ifstream(fileName).good())
    {
      std::cout << "No correction file found." << std::endl
		<< "Re-centering and shifting will not be performed." << std::endl;
    }
  else
    {
      correctionFile = TFile::Open(fileName.Data(), "READ");
      
      TKey *key;
      TIter next(correctionFile->GetListOfKeys());
      TProfile *profile;

      while( (key = (TKey*)next()) )
	{
	  TClass *cl = gROOT->GetClass(key->GetClassName());
	  
	  if (cl->InheritsFrom("TProfile"))
	    {
	      profile = (TProfile*)key->ReadObj();
	      if (profile->GetEntries() == 0)
		{
		  std::cout << "TProfiles are empty!" << std::endl
			    << "Re-centering will be performed and TProfiles will be filled." << std::endl;
		  runIteration = 1;
		  break;
		}
	      else if (profile->GetEntries() != 0)
		{
		  std::cout << "Non-empty TProfiles found!" << std::endl
			    << "Re-centering and event plane shifting will be performed." << std::endl;
		  runIteration = 2;
		  break;
		}
	    }
	}
    }
}// End setCorrectionFileAndRunIteration()

void SetupAttributes::setTPCEfficiencyFile(TString fileName)
{
  if (!ifstream(fileName).good()) 
    { 
      std::cout << "One or both efficiency files missing! All efficiencies will default to 1!" << std::endl; 
    }
  else 
    { 
      pikpEfficiencyFile = TFile::Open(fileName.Data(), "READ");
      tpcEfficienciesFound = true;
      std::cout << "TPC efficiency files were found!" << std::endl; 

      h2_tracking_pp = (TH2D*)pikpEfficiencyFile->Get("h2_ratio_pp");
      h2_tracking_pm = (TH2D*)pikpEfficiencyFile->Get("h2_ratio_pm");
      h2_tracking_kp = (TH2D*)pikpEfficiencyFile->Get("h2_ratio_kp");
      h2_tracking_km = (TH2D*)pikpEfficiencyFile->Get("h2_ratio_km");	  
      h2_tracking_pr = (TH2D*)pikpEfficiencyFile->Get("h2_ratio_pr");

      if (!h2_tracking_pp ||
	  !h2_tracking_pm ||
	  !h2_tracking_kp ||
	  !h2_tracking_km ||
	  !h2_tracking_pr)
	{ 
	  std::cout << "FAILED TO RETRIEVE ALL EFFICIENCY HISTOGRAMS!" << std::endl
		    << "ALL EFFICIENCIES WILL DEFAULT TO 1!" << std::endl;

	  tpcEfficienciesFound = false;
	}
    }
}// End setTPCEfficiencyFile()


void SetupAttributes::setProtonEfficiencyFile(TString fileName)
{
  if (!ifstream(fileName).good()) 
    { 
      std::cout << "Proton efficiency file missing! All efficiencies will default to 1!" << std::endl; 
    }
  else 
    { 
      protonEfficiencyFile = TFile::Open(fileName.Data(), "READ");
      tpcEfficienciesFound = true;
      std::cout << "TPC efficiency files were found!" << std::endl; 

      h2_tracking_pp = 0;
      h2_tracking_pm = 0;
      h2_tracking_kp = 0;
      h2_tracking_km = 0;
      h2_tracking_pr = (TH2D*)protonEfficiencyFile->Get("h2_ratio_pr");

      if (!h2_tracking_pr)
	{ 
	  std::cout << "FAILED TO RETRIEVE PROTON EFFICIENCY HISTOGRAM!" << std::endl
		    << "ALL EFFICIENCIES WILL DEFAULT TO 1!" << std::endl;

	  tpcEfficienciesFound = false;
	}
    }
}// End setProtonEfficiencyFile()


void SetupAttributes::setTOFEfficiencyFile(TString fileName)
{
  if (!ifstream(fileName).good()) 
    { 
      std::cout << "TOF efficiency file missing! All TOF efficiencies will default to 1!" << std::endl; 
    }
  else 
    { 
      tofEfficiencyFile = TFile::Open(fileName.Data(), "READ");
      tofEfficienciesFound = true;
      std::cout << "TOF efficiency file was found!" << std::endl; 

      h2_ratio_tof = (TH2D*)tofEfficiencyFile->Get("h2_ratio_tof");

      if (!h2_ratio_tof)
	{ 
	  std::cout << "FAILED TO RETRIEVE TOF EFFICIENCY HISTOGRAM!" << std::endl
		    << "ALL EFFICIENCIES WILL DEFAULT TO 1!" << std::endl;

	  tofEfficienciesFound = false;
	}
    }
}// End setTOFEfficiencyFile()


void SetupAttributes::setv1WeightsFile(TString fileName)
{
  if (!ifstream(fileName).good()) 
    { 
      std::cout << "No v1 weight file was found! No v1 weights will be applied!" << std::endl; 
    }
  else 
    { 
      v1WeightsInputFile = TFile::Open(fileName.Data(), "READ"); 
      v1WeightsFound = true;
      std::cout << "v1 weight file found! Weights will be applied." << std::endl; 

      p2_TPCv1Weights = (TProfile2D*)v1WeightsInputFile->Get("p2_v1_eta_cent_TPC");
      p2_EPDv1Weights = (TProfile2D*)v1WeightsInputFile->Get("p2_v1_ring_cent_EPD");

      if (!p2_TPCv1Weights || !p2_EPDv1Weights)
	{
	  std::cout << "Failed to retrieve v1 weights from file!" << std::endl
		    << "v1 weights will not be applied!" << std::endl;
	  v1WeightsFound = false;
	}
    }
}// End setv1WeightsFile()


void SetupAttributes::setResolutionFile(TString fileName)
{
  if (!ifstream(fileName).good()) 
    std::cout << "No resolution file was found!" << std::endl; 
  else 
    { 
      resolutionFile = TFile::Open(fileName.Data(), "READ"); 
      resolutionsFound = true;
      std::cout << "Resolution file found!" << std::endl; 
    }
}// End setResolutionFile()
