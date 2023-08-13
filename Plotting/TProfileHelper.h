#ifndef TPROFILEHELPER_H
#define TPROFILEHELPER_H

#include <iostream>
#include <algorithm>
#include <vector>
#include "TH1D.h"
#include "TMath.h"


struct TProfileHelper
{
  std::vector<Double_t> contents;
  std::vector<Double_t> uncertainties;
  std::vector<Double_t> entries;
  Double_t lowBound;
  Double_t highBound;
  Double_t oddAverage;
  Double_t oddAverageUncertainty;
  Double_t evenAverage;
  Double_t evenAverageUncertainty;
  Int_t nBins;
    
  // Used to start, or REstart, this helper struct
  void initialize(TProfile* p1)
  {
    nBins = p1->GetNbinsX();
    lowBound = p1->GetXaxis()->GetXmin();
    highBound = p1->GetXaxis()->GetXmax();

    std::vector<Double_t>().swap(contents);
    std::vector<Double_t>().swap(uncertainties);
    std::vector<Double_t>().swap(entries);
      
    contents.resize(nBins);
    uncertainties.resize(nBins);
    entries.resize(nBins);

    for (int bin = 1; bin <= nBins; bin++)
      {
	contents.at(bin-1) = p1->GetBinContent(bin);
	uncertainties.at(bin-1) = p1->GetBinError(bin);
	entries.at(bin-1) = p1->GetBinEffectiveEntries(bin);
      }
  }

  void initialize(std::vector<Double_t> initialContents,
		  std::vector<Double_t> initialUncertainties,
		  std::vector<Double_t> initialEntries,
		  Double_t initialLowBound,
		  Double_t initialHighBound)
  {
    std::vector<Double_t>().swap(contents);
    std::vector<Double_t>().swap(uncertainties);
    std::vector<Double_t>().swap(entries);
    
    contents = initialContents;
    uncertainties = initialUncertainties;
    entries = initialEntries;

    nBins = contents.size();
    lowBound = initialLowBound;
    highBound = initialHighBound;
  }

  void reverse()
  {
    std::reverse(contents.begin(), contents.end());
    std::reverse(uncertainties.begin(), uncertainties.end());
    std::reverse(entries.begin(), entries.end());
  }

  void scale(Double_t factor)
  {
    for (int bin = 1; bin <= nBins; bin++)
      {
	contents.at(bin-1) *= factor;
	uncertainties.at(bin-1) *= TMath::Abs(factor);
	entries.at(bin-1) *= TMath::Abs(factor);
      }
  }

  TH1D* getHistogram(TString name)
  {
    TH1D* histoVersion = new TH1D(name, "", nBins, lowBound, highBound);
      
    for (int bin = 1; bin <= nBins; bin++)
      {
	histoVersion->SetBinContent(bin, contents.at(bin-1));
	histoVersion->SetBinError(bin, uncertainties.at(bin-1));
      }
      
    return histoVersion;
  }

  TGraphErrors* getGraph(TString name)
  {
    TH1D* histo = getHistogram(name);
    TGraphErrors* graphVersion = new TGraphErrors(histo);
    delete histo;
    return graphVersion;
  }

  void copy(TProfileHelper profileToCopy)
  {
    nBins = profileToCopy.nBins;

    std::vector<Double_t>().swap(contents);
    std::vector<Double_t>().swap(uncertainties);
    std::vector<Double_t>().swap(entries);
      
    contents.resize(nBins);
    uncertainties.resize(nBins);
    entries.resize(nBins);

    for (int element = 0; element < nBins; element++)
      {
	contents.at(element) = profileToCopy.contents.at(element);
	uncertainties.at(element) = profileToCopy.uncertainties.at(element);
	entries.at(element) = profileToCopy.entries.at(element);
      }
  }

  void removeEvenComponent()
  {
    // Make a copy of the values that is Reversed and Negative
    std::vector<Double_t> contentsRN = contents;
    std::vector<Double_t> uncertaintiesRN = uncertainties;
    std::vector<Double_t> entriesRN = entries;

    std::reverse(contentsRN.begin(), contentsRN.end());
    std::reverse(uncertaintiesRN.begin(), uncertaintiesRN.end());
    std::reverse(entriesRN.begin(), entriesRN.end());

    for (int element = 0; element < contentsRN.size(); element++)
	contentsRN.at(element) *= -1.0;
    ////

    // Replace contents with the ODD component only
    for (int element = 0; element < contents.size(); element++)
      {
	contents.at(element) = ( contents.at(element) + contentsRN.at(element) ) / 2.0;
	uncertainties.at(element) = 0.5 * TMath::Sqrt( TMath::Power(uncertainties.at(element), 2.0) + TMath::Power(uncertaintiesRN.at(element), 2.0) );
	entries.at(element) = ( entries.at(element) + entriesRN.at(element) ) / 2.0;
      }

    std::vector<Double_t>().swap(contentsRN);
    std::vector<Double_t>().swap(uncertaintiesRN);
    std::vector<Double_t>().swap(entriesRN);
    ////

    // Flip the sign of half and save the average
    std::vector<Double_t> contentsCopy = contents;
    std::vector<Double_t> uncertaintiesCopy = uncertainties;
    std::vector<Double_t> entriesCopy = entries;

    for (int element = 0; element < contentsCopy.size()/2; element++)
	contentsCopy.at(element) *= -1.0;

    Double_t numerator = 0.0;
    Double_t denominator = 0.0;
    Double_t sumW = 0.0;
    Double_t uncertaintySum = 0.0;
      
    for (int element = 0; element < contentsCopy.size(); element++)
      {
	if (contentsCopy.at(element) == 0.0 &&
	    uncertaintiesCopy.at(element) == 0.0 &&
	    entriesCopy.at(element) == 0.0)
	  continue;

	// weighted average
	numerator += entriesCopy.at(element) * contentsCopy.at(element);
	sumW += entriesCopy.at(element);
	uncertaintySum += TMath::Power( entriesCopy.at(element) * uncertaintiesCopy.at(element), 2.0 );

	// non-weighted average
	//numerator += contents.at(element);
	//denominator += 1.0;;
      }

    oddAverage = numerator / sumW;
    oddAverageUncertainty = (1.0/sumW) * TMath::Sqrt( uncertaintySum );

    std::vector<Double_t>().swap(contentsCopy);
    std::vector<Double_t>().swap(uncertaintiesCopy);
    std::vector<Double_t>().swap(entriesCopy);
  }


  void removeOddComponent()
  {
    // Make a copy of the values that is Reversed and Negative
    std::vector<Double_t> contentsRN = contents;
    std::vector<Double_t> uncertaintiesRN = uncertainties;
    std::vector<Double_t> entriesRN = entries;

    std::reverse(contentsRN.begin(), contentsRN.end());
    std::reverse(uncertaintiesRN.begin(), uncertaintiesRN.end());
    std::reverse(entriesRN.begin(), entriesRN.end());

    for (int element = 0; element < contentsRN.size(); element++)
	contentsRN.at(element) *= -1.0;
    ////

    // Replace contents with the EVEN component only
    for (int element = 0; element < contents.size(); element++)
      {
	contents.at(element) = ( contents.at(element) - contentsRN.at(element) ) / 2.0;
	uncertainties.at(element) = 0.5 * TMath::Sqrt( TMath::Power(uncertainties.at(element), 2.0) + TMath::Power(uncertaintiesRN.at(element), 2.0) );
	entries.at(element) = TMath::Abs( entries.at(element) - entriesRN.at(element) ) / 2.0;
      }

    std::vector<Double_t>().swap(contentsRN);
    std::vector<Double_t>().swap(uncertaintiesRN);
    std::vector<Double_t>().swap(entriesRN);
    ////

    // Save the average
    std::vector<Double_t> contentsCopy = contents;
    std::vector<Double_t> uncertaintiesCopy = uncertainties;
    std::vector<Double_t> entriesCopy = entries;

    Double_t numerator = 0.0;
    Double_t denominator = 0.0;
    Double_t sumW = 0.0;
    Double_t uncertaintySum = 0.0;
      
    for (int element = 0; element < contentsCopy.size(); element++)
      {
	if (contentsCopy.at(element) == 0.0 &&
	    uncertaintiesCopy.at(element) == 0.0 &&
	    entriesCopy.at(element) == 0.0)
	  continue;

	// weighted average
	numerator += entriesCopy.at(element) * contentsCopy.at(element);
	sumW += entriesCopy.at(element);
	uncertaintySum += TMath::Power( entriesCopy.at(element) * uncertaintiesCopy.at(element), 2.0 );

	// non-weighted average
	//numerator += contents.at(element);
	//denominator += 1.0;;
      }

    evenAverage = numerator / sumW;
    evenAverageUncertainty = (1.0/sumW) * TMath::Sqrt( uncertaintySum );

    std::vector<Double_t>().swap(contentsCopy);
    std::vector<Double_t>().swap(uncertaintiesCopy);
    std::vector<Double_t>().swap(entriesCopy);
  }

  ~TProfileHelper()
  {
    std::vector<Double_t>().swap(contents);
    std::vector<Double_t>().swap(uncertainties);
    std::vector<Double_t>().swap(entries);
  }
}; // End struct TProfileHelper



#endif
