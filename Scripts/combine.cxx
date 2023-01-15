void combine(/*TString fitType,*/TString jobId, Int_t fileNumber, Int_t keyJumpIndex = 1) 
{
  Int_t quarterMark = (Int_t)fileNumber/4;
  Int_t halfMark = (Int_t)fileNumber/2;
  Int_t threeQuarterMark = (Int_t)fileNumber*3/4;

  Bool_t quarterMessageGiven = false;
  Bool_t halfMessageGiven = false;
  Bool_t threeQuarterMessageGiven = false;

  // output file
  char name[1000];

  TString outputPicoFileName = jobId;
  outputPicoFileName.Append(".picoDst.result.combined.root");
  TFile *outputPicoFile = new TFile(outputPicoFileName,"recreate");

  // input file
  TFile *tmpFile[6000]; //Array of root files
  TString tmpFilename;
  TObjArray *ObjArray = new TObjArray(500); 
  Int_t Nkeys = 0;
  Int_t Nhist = 0;
  Int_t errorCount = 0;

  for(Int_t i=0; i < fileNumber; i++) 
    {
      sprintf(name,"_%d.picoDst.result.root",i);
      tmpFilename = jobId;
      tmpFilename.Append(name);
      tmpFile[i] = new TFile(tmpFilename,"read");

      if(tmpFile[i]->IsZombie()) { errorCount++; continue; }

      if(i == 0) Nkeys = tmpFile[i]->GetNkeys();
      if(Nkeys == 0) { std::cout << "No keys!" << std::endl; return; }

      TList *list = tmpFile[i]->GetListOfKeys();
      TIter next((TList*)list);
      TKey  *key;
      Int_t histCounter = 0;
      Int_t keyJumper  = 0;

      while ( (key = (TKey*)next()) ) 
	{
	  keyJumper++;
	  if(keyJumper < keyJumpIndex) continue;

	  TObject *obj = tmpFile[i]->Get(key->GetName());
	  if(!obj) { std::cout << "Could not get object out of file " << i << std::endl; return; }
	  if(obj->IsA() == TDirectory::Class()) continue;
	  if(obj->InheritsFrom(TH1::Class())) 
	    {
	      //std::cout << key->GetName() << std::endl;
	      if(i == 0) ObjArray->Add((TH1*)obj);
	      if(i >= 1) ((TH1*)ObjArray->At(histCounter))->Add((TH1*)obj); //Add like histograms together
	      histCounter++;
	    }
	}

      if(i == 0) Nhist = histCounter;

      // clean up
      if(i >= 1)
	{
	  tmpFile[i]->Close();
	  delete tmpFile[i];
	}

      //std::cout << "File " << i << " read in." << std::endl;
      if (i > quarterMark && !quarterMessageGiven) { std::cout << "25% done..." << std::endl; quarterMessageGiven = true; }
      else if (i > halfMark && !halfMessageGiven) { std::cout << "50% done..." << std::endl; halfMessageGiven = true; }
      else if (i > threeQuarterMark && !threeQuarterMessageGiven) { std::cout << "75% done..." << std::endl; threeQuarterMessageGiven = true; }
    }

  // write merged histograms
  std::cout << "Writing merged histograms." << std::endl;

  for(Int_t k=0; k <= Nhist-keyJumpIndex; k++) 
    {
      outputPicoFile->cd();
      ((TH1*)ObjArray->At(k))->Write();
    }

  ObjArray->Clear();
  tmpFile[0]->Close();
  delete tmpFile[0];

  Double_t errorPercent = (Double_t)errorCount/fileNumber * 100;

  std::cout << "Missing files: " << errorCount << " (" << errorPercent << "%)" << std::endl;
  std::cout << "Done!" << std::endl;
}
