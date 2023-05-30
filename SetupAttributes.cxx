#include <iostream>
#include "SetupAttributes.h"

Int_t SetupAttributes::getRunIteration(TFile* correctionFile)
{
  Int_t runIteration = 0;

  if (!correctionFile)
    {
      std::cout << "No correction file found." << std::endl
		<< "Re-centering and shifting will not be performed." << std::endl;
    }
  else
    {
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

  return runIteration;
}// End getRunIteration()


Bool_t SetupAttributes::setResolutionFile(TFile*& filePointerAddress, TString fileName)
{
  Bool_t resolutionsFound = false;
  filePointerAddress = TFile::Open(fileName, "READ"); 
  if (!filePointerAddress) 
    std::cout << "No resolution file was found!" << std::endl; 
  else 
    { 
      resolutionsFound = true;
      std::cout << "Resolution file found!" << std::endl; 
    }
  return resolutionsFound;
}// End setResolutionFile()


