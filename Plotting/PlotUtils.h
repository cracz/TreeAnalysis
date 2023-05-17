#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include <iostream>
#include "TH1D.h"
#include "TMath.h"

namespace PlotUtils
{
  TH1D* flipHisto(TH1D* histo)
  {
    TH1D *h_flipped = (TH1D*)histo->Clone((TString)histo->GetName() + "_flip");
    Int_t bins = histo->GetNbinsX();
    
    Int_t j = 1;
    for (int i = bins; i >= 1; i--)
      {
	h_flipped->SetBinContent(j, histo->GetBinContent(i));
	h_flipped->SetBinError(j, histo->GetBinError(i));
	j++;
      }

    return h_flipped;
  };// End flipHisto()


  // This function will remove any points with y = ey = 0
  void cleanGraph(TGraphErrors* graph)
  {
    // Location of points that need to be removed.
    // Once one is removed, the position of the rest must
    // be updated before removal.
    std::vector<int> emptyPointPositions;

    const Int_t n = graph->GetN();

    // Find the location of each empty point
    for (int i = 0; i < n; i++)
      {
	if (graph->GetPointY(i) == 0.0 && graph->GetErrorY(i) == 0.0)
	  emptyPointPositions.push_back(i);
      }

    // Loop through positions of empty points
    for (int j = 0; j < emptyPointPositions.size(); j++)
      {
	graph->RemovePoint(emptyPointPositions.at(j)); // remove empty point
	
	// Update the new position of the rest of the empty points
	for (int k = j+1; k < emptyPointPositions.size(); k++)
	  {
	    emptyPointPositions.at(k) -= 1;
	  }
      }
  };// End cleanGraph()


  TGraphErrors* flipGraph(TGraphErrors* graph)
  {
    TString newName = (TString)graph->GetName() + "_flip";
    const Int_t n = graph->GetN();

    // constant
    Double_t* x = graph->GetX();
    Double_t* ex = graph->GetEX();
    if (!ex) ex = graph->GetEXhigh();

    // flipping order
    Double_t* y = graph->GetY();
    Double_t* ey = graph->GetEY();
    if (!ey) ey = graph->GetEYhigh();

    // Final values
    Double_t xf[n];
    Double_t yf[n];
    Double_t exf[n];
    Double_t eyf[n]; 

    Int_t j = 0;
    for (int i = n-1; i >= 0; i--)
      {
	xf[j] = -1.0*x[i];
	exf[j] = ex[i];
	
	yf[j] = y[i];
	eyf[j] = ey[i];
	
	j++;
      }

    TGraphErrors* g_flipped = new TGraphErrors(n, xf, yf, exf, eyf);
    g_flipped->SetName(newName);
    return g_flipped;
  };// End flipGraph()



  
  // Cut off everything above 60% centrality
  TH1D* trimCentralityPlot(TH1D *histo)
  {
    TH1D *h_trimmed = (TH1D*)histo->Clone();
    Int_t oldBins = histo->GetNbinsX();
    Int_t newBins = (oldBins == 16) ? 12 : 6; // Usually 12, 6 for rebinned kaons.
    h_trimmed->SetBins(newBins, 0, 60);
    h_trimmed->GetXaxis()->SetTitle("Centrality (%)");

    for (int i = 1; i < newBins; i++)
      {
	h_trimmed->SetBinContent(i, histo->GetBinContent(i));
	h_trimmed->SetBinError(i, histo->GetBinError(i));
      }

    return h_trimmed;
  };// End trimCentralityPlot()



  
  // Cut off everything below y_cm = 0
  TH1D* trimRapidityPlot(TH1D* histo)
  {
    TH1D *h_trimmed = (TH1D*)histo->Clone();
    Int_t oldBins = histo->GetNbinsX();
    Int_t newBins = (oldBins == 20) ? 10 : 5; //Cut Nbins in half, 5 for rebinned kaons.
    Int_t binOffset = newBins; // old and new are offset by this amount
    h_trimmed->SetBins(newBins, 0, 1);
    h_trimmed->GetXaxis()->SetTitle("y-y_{mid}");

    for (int i = 1; i <= newBins; i++)
      {
	h_trimmed->SetBinContent(i, histo->GetBinContent(i+binOffset));
	h_trimmed->SetBinError(i, histo->GetBinError(i+binOffset));
      }

    return h_trimmed;
  };// End trimRapidityPlot()


  
  // Take a histogram that is symmetric across
  // x = 0 and return the odd component of that histogram.
  TH1D* getOddComponent(TH1D* hist)
  {
    TH1D* h_reflected = flipHisto(hist);
    h_reflected->SetName("h_reflected");

    TH1D* h_reflectedAndNegative = (TH1D*)h_reflected->Clone("h_reflectedAndNegative");
    h_reflectedAndNegative->Scale(-1.0);

    TH1D* h_oddComponent = (TH1D*)hist->Clone("h_oddComponent");
    h_oddComponent->Add(h_reflectedAndNegative, 1.0);

    for (int i = 1; i <= h_oddComponent->GetNbinsX(); i++)
      {
	h_oddComponent->SetBinContent(i, h_oddComponent->GetBinContent(i)/2.0);
	h_oddComponent->SetBinError(i, h_oddComponent->GetBinError(i)/2.0);
      }

    return h_oddComponent;
  };// End getOddComponent()


  

  // Take a histogram that is symmetric across
  // x = 0 and return the even component of that histogram.  
  TH1D* getEvenComponent(TH1D* hist)
  {
    TH1D* h_reflected = flipHisto(hist);
    h_reflected->SetName("h_reflected");

    TH1D* h_reflectedAndNegative = (TH1D*)h_reflected->Clone("h_reflectedAndNegative");
    h_reflectedAndNegative->Scale(-1.0);

    TH1D* h_evenComponent = (TH1D*)hist->Clone("h_evenComponent");
    h_evenComponent->Add(h_reflectedAndNegative, -1.0);

    for (int i = 1; i <= h_evenComponent->GetNbinsX(); i++)
      {
	h_evenComponent->SetBinContent(i, h_evenComponent->GetBinContent(i)/2.0);
	h_evenComponent->SetBinError(i, h_evenComponent->GetBinError(i)/2.0);
      }

    return h_evenComponent;
  };// End getEvenComponent()


  // Add two graphs with a multiplicative factor on the second.
  // Both graphs must be guaranteed to have the same range and same number of
  // points before using this!
  TGraphErrors* addGraphs(TGraphErrors* g1, TGraphErrors* g2, Double_t factor)
  {
    // g1
    const Int_t n1 = g1->GetN();

    Double_t* x1 = g1->GetX();
    Double_t* ex1 = g1->GetEX();
    if (!ex1) ex1 = g1->GetEXhigh();
    Double_t x1Beginning = x1[0];
    Double_t x1End = x1[n1-1];

    Double_t* y1 = g1->GetY();
    Double_t* ey1 = g1->GetEY();
    if(!ey1) ey1 = g1->GetEYhigh();

    // g2
    const Int_t n2 = g2->GetN();

    Double_t* x2 = g2->GetX();
    Double_t* ex2 = g2->GetEX();
    if (!ex2) ex2 = g2->GetEXhigh();
    Double_t x2Beginning = x2[0];
    Double_t x2End = x2[n2-1];

    Double_t* y2 = g2->GetY();
    Double_t* ey2 = g2->GetEY();
    if(!ey2) ey2 = g2->GetEYhigh();

    if (n1 != n2)
      {
	std::cout << "Added graphs do not have the same binning!" << std::endl;
	return 0;
      }

    // Final values
    Double_t yf[n1];
    Double_t eyf[n1]; 

    for (int i = 0; i < n1; i++)
      {
	yf[i] = y1[i] + (factor * y2[i]);
	eyf[i] = TMath::Sqrt(ey1[i]*ey1[i] + factor*factor*ey2[i]*ey2[i]);
      }

    TGraphErrors* g_sum = new TGraphErrors(n1, x1, yf, ex1, eyf);
    g_sum->SetName("g_sum");
    return g_sum;
  }; //End addGraphs()
  
  /*
  TGraphErrors* addGraphs(TGraphErrors* g1, TGraphErrors* g2, Double_t factor)
  {
    // g1
    const Int_t n1 = g1->GetN();

    Double_t* x1 = g1->GetX();
    Double_t* ex1 = g1->GetEXhigh();
    if (!ex1) ex1 = g1->GetEX();
    Double_t x1Beginning = x1[0];
    Double_t x1End = x1[n1-1];

    Double_t* y1 = g1->GetY();
    Double_t* ey1 = g1->GetEYhigh();
    if(!ey1) ey1 = g1->GetEY();

    // g2
    const Int_t n2 = g2->GetN();

    Double_t* x2 = g2->GetX();
    Double_t* ex2 = g2->GetEXhigh();
    if (!ex2) ex2 = g2->GetEX();
    Double_t x2Beginning = x2[0];
    Double_t x2End = x2[n2-1];

    Double_t* y2 = g2->GetY();
    Double_t* ey2 = g2->GetEYhigh();
    if(!ey2) ey2 = g2->GetEY();

    std::cout << "x1Beginning = " << x1Beginning << std::endl;
    std::cout << "x1End = " << x1End << std::endl;

    std::cout << "x2Beginning = " << x2Beginning << std::endl;
    std::cout << "x2End = " << x2End << std::endl;

    if (n1 != n2)
      {
	std::cout << "Added graphs do not have the same binning!" << std::endl;
	return 0;
      }

    // Start and end the sum where the bins match up on the x-axis,
    // ignore x values where the two graphs don't both have values at the same time.
    Int_t firstBin;
    Int_t x1FirstBin;
    Int_t x2FirstBin;
    Int_t lastBin;
    Int_t x1LastBin;
    Int_t x2LastBin;
    Int_t binsDropped;

    ////// Find beginning point
    // Here, g1 x-values extend farther down than g2.
    if (x1Beginning < x2Beginning)
      {
	std::cout << "HEREA" << std::endl;
	for (int i = 0; i < n1; i++)
	  {
	    if (x1[i] == x2[0])
	      {
		x1FirstBin = i; // first bin for x1
		x2FirstBin = 0;
		binsDropped = i;
		break;
	      }
	  }
      }
    // Here, g2 x-values extend farther down than g1.
    else if (x2Beginning < x1Beginning)
      {
	std::cout << "HEREB" << std::endl;
	for (int i = 0; i < n2; i++)
	  {
	    if (x2[i] == x1[0])
	      {
		x2FirstBin = i; // first bin for x2
		x1FirstBin = 0;
		binsDropped = i;
		break;
	      }
	  } 
      }
    else
      {
	std::cout << "HERE1" << std::endl;
	x1FirstBin = 0;
	x2FirstBin = 0;
	binsDropped = 0;
      }
    //////

    ////// Find end point
    // Here, g1 x-values extend farther up than g2.
    if (x1End > x2End)
      {
	for (int i = n1-1; i >= 0; i--)
	  {
	    if (x1[i] == x2[n2-1])
	      {
		x1LastBin = i;
		x2LastBin = n2-1;
		//binsDropped = i; redundant
		break;
	      }
	  }
      }
    // Here, g2 x-values extend farther up than g1.
    else if (x2End > x1End)
      {
	for (int i = n2-1; i >= 0; i--)
	  {
	    if (x2[i] == x1[n1-1])
	      {
		x2LastBin = i;
		x1LastBin = n1-1;
		//binsDropped = i; redundant
		break;
	      }
	  } 
      }
    else
      {
	std::cout << "HERE2" << std::endl;
	x1LastBin = n1-1;
	x2LastBin = n2-1;
	//binsDropped = 0;
      }
    //////


    // Final values
    const Int_t nf = n1-binsDropped;
    std::cout << "nf = " << nf << std::endl;
    Double_t xf[nf];
    Double_t yf[nf];
    Double_t exf[nf];
    Double_t eyf[nf]; 

    Int_t i = 0;
    // i is the index through the original points, from firstBin to lastBin.
    // j is the index of the sum points, from 0 to nf-1.
    for (int j = 0; j < nf; j++)
      {
	xf[j]  = x1[x1FirstBin+i];
	exf[j] = 0.0;     // Ignore x axis uncertainties.

	yf[j] = y1[x1FirstBin+i] + (factor * y2[x2FirstBin+i]);
	eyf[j] = TMath::Sqrt(ey1[x1FirstBin+i]*ey1[x1FirstBin+i] + factor*factor*ey2[x2FirstBin+i]*ey2[x2FirstBin+i]);
	i++;
      }

    TGraphErrors* g_sum = new TGraphErrors(nf, xf, yf, exf, eyf);
    g_sum->SetName("g_sum");
    return g_sum;
  };// End addGraphs()
  */


  void scaleGraph(TGraphErrors* graph, Double_t factor)
  {
    Bool_t normalError = false;
    Bool_t highError = false;
    Bool_t lowError = false;
    
    Double_t* errorCheck = graph->GetEY();
    normalError = true;
    
    if (!errorCheck)
      {
	normalError = false;
	errorCheck = graph->GetEYhigh();
	highError = true;
      }
    if (!errorCheck)
      {
	highError = false;
	errorCheck = graph->GetEYlow();
	lowError = true;
      }
    if (!errorCheck)
      {
	lowError = false;
      }

    if (!normalError && !highError && !lowError)
      {
	std::cout << "scaleGraph() can't find any uncertainties! Scaling aborted!" << std::endl;
	return;
      }

    
    for (int i = 0; i <= graph->GetN(); i++)
      {
	graph->GetY()[i] *= factor;

	if (normalError)
	  graph->GetEY()[i] *= TMath::Abs(factor);
	else if (highError)
	  graph->GetEYhigh()[i] *= TMath::Abs(factor);
	else if (lowError)
	  graph->GetEYlow()[i] *= TMath::Abs(factor);
      }
  };// End scaleGraph()
  


  // Take a graph with errors that is symmetric across
  // x = 0 and return the odd component of that graph.
  TGraphErrors* getOddComponentGraph(TGraphErrors* graph)
  {
    TString newName = (TString)graph->GetName() + "_odd";

    cleanGraph(graph);
    
    TGraphErrors* g_reflected = flipGraph(graph);
    g_reflected->SetName("g_reflected");

    TGraphErrors* g_reflectedAndNegative = (TGraphErrors*)g_reflected->Clone("g_reflectedAndNegative");
    scaleGraph(g_reflectedAndNegative, -1.0);

    TGraphErrors* g_oddComponent = addGraphs(graph, g_reflectedAndNegative, 1.0);
    g_oddComponent->SetName(newName);
    scaleGraph(g_oddComponent, 1.0/2.0);

    return g_oddComponent;
  };// End getOddComponentGraph()


  

  // Take a graph with errors that is symmetric across
  // x = 0 and return the even component of that graph.  
  TGraphErrors* getEvenComponentGraph(TGraphErrors* graph)
  {
    TString newName = (TString)graph->GetName() + "_even";

    cleanGraph(graph);
	
    TGraphErrors* g_reflected = flipGraph(graph);
    g_reflected->SetName("g_reflected");

    TGraphErrors* g_reflectedAndNegative = (TGraphErrors*)g_reflected->Clone("g_reflectedAndNegative");
    scaleGraph(g_reflectedAndNegative, -1.0);

    TGraphErrors* g_evenComponent = addGraphs(graph, g_reflectedAndNegative, -1.0);
    g_evenComponent->SetName(newName);
    scaleGraph(g_evenComponent, 1.0/2.0);

    return g_evenComponent;
  };// End getEvenComponentGraph()



  // Apply resolutions to standard v_n vs centrality plot (histo)
  // USE THIS BEFORE FLIPPING OR TRIMMING THE ORIGINAL HISTOGRAM!!!
  TH1D* applyResolutions(TH1D* histo, TH1D* resolutions)
  {
    TH1D *newHisto = (TH1D*)histo->Clone((TString)histo->GetName() + "_resCor");
    int bins = histo->GetNbinsX();

    for (int i = 1; i <= bins; i++)
      {
	if (resolutions->GetBinContent(i) == 0.0) continue;

	Double_t originalValue = histo->GetBinContent(i);
	Double_t originalError = histo->GetBinError(i);
	Double_t currentResolution = resolutions->GetBinContent(i);
	Double_t currentResolutionError = resolutions->GetBinError(i);

	Double_t newValue = originalValue / currentResolution;
	Double_t newError = newValue * TMath::Sqrt( TMath::Power(originalError/originalValue, 2) +
						    TMath::Power(currentResolutionError/currentResolution, 2) );

	newHisto->SetBinContent(i, newValue);
	newHisto->SetBinError(i, newError);
      }

    return newHisto;
  };// End applyResolutions()



  
  // Apply resolutions to 2D v_n plot with y vs centrality or p_T vs centrality (y-axis vs x-axis)
  // USE THIS BEFORE PROJECTING INTO WIDER CENTRALITY BINS!!!
  TH2D* applyResolutions2D(TH2D* histo, TH1D* resolutions)
  {
    TH2D *newHisto = (TH2D*)histo->Clone((TString)histo->GetName() + "_resCor");
    int nBinsX = histo->GetNbinsX();
    int nBinsY = histo->GetNbinsY();

    for (int xBin = 1; xBin <= nBinsX; xBin++)
      {
	if (resolutions->GetBinContent(xBin) == 0.0) continue;

	Double_t currentResolution = resolutions->GetBinContent(xBin);
	Double_t currentResolutionError = resolutions->GetBinError(xBin);
	
	for (int yBin = 1; yBin <= nBinsY; yBin++)
	  {
	    if (histo->GetBinContent(xBin, yBin) == 0.0) continue;

	    Double_t originalValue = histo->GetBinContent(xBin, yBin);
	    Double_t originalError = histo->GetBinError(xBin, yBin);

	    Double_t newValue = originalValue / currentResolution;
	    Double_t newError = newValue * TMath::Sqrt( TMath::Power(originalError/originalValue, 2) +
							TMath::Power(currentResolutionError/currentResolution, 2) );

	    newHisto->SetBinContent(xBin, yBin, newValue);
	    newHisto->SetBinError(xBin, yBin, newError);
	  }
      }

    return newHisto;
  };// End applyResolutions2D()
}

#endif
