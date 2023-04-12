#ifndef FLOWUTILS_H
#define FLOWUTILS_H

#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"
#include "TMath.h"


namespace FlowUtils
{
  const Int_t I_BAD_VALUE    = -999;
  const Double_t D_BAD_VALUE = -999.0;

  // Custom type to hold info for each good particle
  struct Particle
  {
    Double_t phi;
    Double_t eta;
    Double_t pT;
    Double_t KT;
    Double_t weight;
    Double_t rapidity;
    Double_t EPDTnMIP;

    Int_t EPDring;

    Bool_t ppTag;
    Bool_t pmTag;
    Bool_t kpTag;
    Bool_t kmTag;
    Bool_t prTag;
    Bool_t deTag;
    Bool_t trTag;

    Bool_t isInTpcA;
    Bool_t isInTpcB;
    Bool_t isInEpdA;
    Bool_t isInEpdB;

    void reset()
    {
      phi = D_BAD_VALUE;
      eta = D_BAD_VALUE;
      pT  = D_BAD_VALUE;
      KT  = D_BAD_VALUE;
      weight = 0;
      rapidity = D_BAD_VALUE;
      EPDTnMIP = D_BAD_VALUE;

      EPDring = I_BAD_VALUE;

      ppTag = false;
      pmTag = false;
      kpTag = false;
      kmTag = false;
      prTag = false;
      deTag = false;
      trTag = false;

      isInTpcA = false;
      isInTpcB = false;
      isInEpdA = false;
      isInEpdB = false;
    }
  }; // End struct Particle




  // Custom type to hold important info for every good event
  struct Event
  {
    std::vector<Particle> tpcParticles;
    std::vector<Particle> epdParticles;

    Double_t XnTpc;
    Double_t YnTpc;
    Double_t psiTpc; 
    Double_t XnTpcA;
    Double_t YnTpcA;
    Double_t psiTpcA;       // Overall EP angle without removing autocorrelations
    Double_t XnTpcB;
    Double_t YnTpcB;
    Double_t psiTpcB;
    Double_t XnEpd;
    Double_t YnEpd;
    Double_t psiEpd;
    Double_t XnEpdA;
    Double_t YnEpdA;
    Double_t psiEpdA;
    Double_t XnEpdB;
    Double_t YnEpdB;
    Double_t psiEpdB;

    Int_t nTracksTpc;
    Int_t nTracksTpcA;      // Number of GOOD tracks in the sub-event
    Int_t nTracksTpcB;
    Int_t nHitsEpd;
    Int_t nHitsEpdA;
    Int_t nHitsEpdB;

    Int_t centID;
    Int_t primTracks;   // Number of primary tracks before track cuts (used for centrality)
    bool badEvent;      // Flag for marking events to ignore

    void incrementQvectorTPC(Bool_t ODD_PLANE, Double_t ORDER_M, Double_t Y_MID, Double_t eta, Double_t phi, Double_t weight)
    {
      if (ODD_PLANE)
	{
	  if (eta > Y_MID)        // Account for Q vector sign change past mid-rapidity.
	    {
	      XnTpc += weight * TMath::Cos(ORDER_M * phi);
	      YnTpc += weight * TMath::Sin(ORDER_M * phi);
	    }
	  else if (eta < Y_MID)
	    {
	      XnTpc -= weight * TMath::Cos(ORDER_M * phi);
	      YnTpc -= weight * TMath::Sin(ORDER_M * phi);
	    }
	}
      else
	{
	  XnTpc += weight * TMath::Cos(ORDER_M * phi);
	  YnTpc += weight * TMath::Sin(ORDER_M * phi);
	}
    }// End incrementQvectorTPC()


    void incrementQvectorTPCA(Bool_t ODD_PLANE, Double_t ORDER_M, Double_t Y_MID, Double_t eta, Double_t phi, Double_t weight)
    {
      if (ODD_PLANE)
	{
	  if (eta > Y_MID)        // Account for Q vector sign change past mid-rapidity.
	    {
	      XnTpcA += weight * TMath::Cos(ORDER_M * phi);
	      YnTpcA += weight * TMath::Sin(ORDER_M * phi);
	    }
	  else if (eta < Y_MID)
	    {
	      XnTpcA -= weight * TMath::Cos(ORDER_M * phi);
	      YnTpcA -= weight * TMath::Sin(ORDER_M * phi);
	    }
	}
      else
	{
	  XnTpcA += weight * TMath::Cos(ORDER_M * phi);
	  YnTpcA += weight * TMath::Sin(ORDER_M * phi);
	}
    }// End incrementQvectorTPCA()

    void incrementQvectorTPCB(Bool_t ODD_PLANE, Double_t ORDER_M, Double_t Y_MID, Double_t eta, Double_t phi, Double_t weight)
    {
      if (ODD_PLANE)
	{
	  if (eta > Y_MID)        // Account for Q vector sign change past mid-rapidity.
	    {
	      XnTpcB += weight * TMath::Cos(ORDER_M * phi);
	      YnTpcB += weight * TMath::Sin(ORDER_M * phi);
	    }
	  else if (eta < Y_MID)
	    {
	      XnTpcB -= weight * TMath::Cos(ORDER_M * phi);
	      YnTpcB -= weight * TMath::Sin(ORDER_M * phi);
	    }
	}
      else
	{
	  XnTpcB += weight * TMath::Cos(ORDER_M * phi);
	  YnTpcB += weight * TMath::Sin(ORDER_M * phi);
	}
    }// End incrementQvectorTPCB()

    void incrementQvectorEPD(Bool_t ODD_PLANE, Double_t ORDER_M, Double_t Y_MID, Double_t eta, Double_t phi, Double_t weight)
    {
      if (ODD_PLANE)
	{
	  if (eta > Y_MID)        // Account for Q vector sign change past mid-rapidity.
	    {
	      XnEpd += weight * TMath::Cos(ORDER_M * phi);
	      YnEpd += weight * TMath::Sin(ORDER_M * phi);
	    }
	  else if (eta < Y_MID)
	    {
	      XnEpd -= weight * TMath::Cos(ORDER_M * phi);
	      YnEpd -= weight * TMath::Sin(ORDER_M * phi);
	    }
	}
      else
	{
	  XnEpd += weight * TMath::Cos(ORDER_M * phi);
	  YnEpd += weight * TMath::Sin(ORDER_M * phi);
	}
    }// End incrementQvectorEPD()

    void incrementQvectorEPDA(Bool_t ODD_PLANE, Double_t ORDER_M, Double_t Y_MID, Double_t eta, Double_t phi, Double_t weight)
    {
      if (ODD_PLANE)
	{
	  if (eta > Y_MID)        // Account for Q vector sign change past mid-rapidity.
	    {
	      XnEpdA += weight * TMath::Cos(ORDER_M * phi);
	      YnEpdA += weight * TMath::Sin(ORDER_M * phi);
	    }
	  else if (eta < Y_MID)
	    {
	      XnEpdA -= weight * TMath::Cos(ORDER_M * phi);
	      YnEpdA -= weight * TMath::Sin(ORDER_M * phi);
	    }
	}
      else
	{
	  XnEpdA += weight * TMath::Cos(ORDER_M * phi);
	  YnEpdA += weight * TMath::Sin(ORDER_M * phi);
	}
    }// End incrementQvectorEPDA()

    void incrementQvectorEPDB(Bool_t ODD_PLANE, Double_t ORDER_M, Double_t Y_MID, Double_t eta, Double_t phi, Double_t weight)
    {
      if (ODD_PLANE)
	{
	  if (eta > Y_MID)        // Account for Q vector sign change past mid-rapidity.
	    {
	      XnEpdB += weight * TMath::Cos(ORDER_M * phi);
	      YnEpdB += weight * TMath::Sin(ORDER_M * phi);
	    }
	  else if (eta < Y_MID)
	    {
	      XnEpdB -= weight * TMath::Cos(ORDER_M * phi);
	      YnEpdB -= weight * TMath::Sin(ORDER_M * phi);
	    }
	}
      else
	{
	  XnEpdB += weight * TMath::Cos(ORDER_M * phi);
	  YnEpdB += weight * TMath::Sin(ORDER_M * phi);
	}
    }// End incrementQvectorEPDB()

    void reset()//Reset all values in the struct to reuse
    {
      std::vector<Particle>().swap(tpcParticles);
      std::vector<Particle>().swap(epdParticles);

      XnTpc = 0;
      YnTpc = 0;
      psiTpc = D_BAD_VALUE;        //Just some number to use that is out of bounds
      XnTpcA = 0;
      YnTpcA = 0;
      psiTpcA = D_BAD_VALUE;
      XnTpcB = 0;
      YnTpcB = 0;
      psiTpcB = D_BAD_VALUE;
      XnEpd = 0;
      YnEpd = 0;
      psiEpd = D_BAD_VALUE;
      XnEpdA = 0;
      YnEpdA = 0;
      psiEpdA = D_BAD_VALUE;
      XnEpdB = 0;
      YnEpdB = 0;
      psiEpdB = D_BAD_VALUE;

      nTracksTpc = 0;
      nTracksTpcA = 0;
      nTracksTpcB = 0;
      nHitsEpd = 0;
      nHitsEpdA = 0;
      nHitsEpdB = 0;

      badEvent  = false;
      primTracks = 0;
      centID = I_BAD_VALUE;
    }
  }; // End struct Event


  struct NewNSigmaProton3p22GeV
  {
    //**************************************************************************************************
    //**************** The following code can be used to implement the new nσ_p ************************
    //********************************************** Code **********************************************
  
    const static int numRapidityWindows = 17;
    const static int numRapidityWindowsToUsePointwiseCorrection = 3;
    TF1* bichselFunction_mean[numRapidityWindows];
    TF1* bichselFunction_1sig[numRapidityWindows];
  
    TGraph* mean_graph[3];
    TGraph* sig1_graph[3];
    
    // Pointwise-correction parameters to be used at low momenta and low rapidities 
    double lowestMomWhichHasAFit[3];
    double lowestMomUsedForBichselFit[3];
    
    double pointwise_momentum_0[20] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195};
    double pointwise_mean_0[20] = {4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 
				   4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.0297, 3.98871, 3.93897, 3.88467, 3.82839, 3.77154};
    double pointwise_1sig_0[20] = {4.01181, 4.01181, 4.01181, 4.01181, 4.01181, 4.01181, 4.01181, 4.01181, 
				   4.01181, 4.01181, 4.01181, 4.01181, 4.01181, 4.01181, 3.96494, 3.92333, 3.87653, 3.82684, 3.77564, 3.72078};
    
    double pointwise_momentum_1[24] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 
				       0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235};
    double pointwise_mean_1[24] = {3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 
				   3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.84422, 3.79324, 3.74542, 3.69679, 3.64574, 3.59203, 3.537};
    double pointwise_1sig_1[24] = {3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 
				   3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.77462, 3.72506, 3.68064, 3.63513, 3.58583, 3.53271, 3.47861};
    
    double pointwise_momentum_2[32] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 
				       0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235, 0.245, 0.255, 0.265, 0.275, 0.285, 0.295, 0.305, 0.315};
    double pointwise_mean_2[32] = {3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 
				   3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 
				   3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.40009, 3.35262, 3.30464, 3.25334, 3.20198, 3.1514, 3.10295};
    double pointwise_1sig_2[32] = {3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 
				   3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 
				   3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.33009, 3.28262, 3.23464, 3.18334, 3.13342, 3.08537, 3.03957};
    

    // Parameters from Bichsel curve fits (H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195) 
    double meanBichselParameters[numRapidityWindows][6] = {
      {-5.141330,-4.605085,1.986481,0.951782,1.908647,0.938000},
      {-0.174776,-1.929899,1.681636,1.133727,2.534694,0.938000},
      {-5.128737,-4.605543,2.012423,0.921413,2.070305,0.938000},
      {-5.094534,-4.660380,1.951221,0.933111,2.088334,0.938000},
      {-5.069162,-4.721723,1.897411,0.920773,2.138482,0.938000},
      {-5.072757,-4.779960,1.874566,0.727895,2.358116,0.938000},
      {-5.047109,-4.832496,1.840607,0.403700,2.646944,0.938000},
      {-0.144720,-2.946088,0.989209,0.101381,3.264879,0.938000},
      {2.141951,-2.171655,0.497909,-0.299899,3.761835,0.938000},
      {-5.033564,-5.124145,1.688028,0.383665,2.885763,0.938000},
      {-4.927111,-5.316412,1.571083,1.245790,2.441705,0.938000},
      {-4.880357,-5.864842,1.377256,1.406169,2.670642,0.938000},
      {-4.185822,-6.984285,0.997941,-0.176737,4.723153,0.938000},
      {-3.833604,-7.124846,0.946295,1.724259,2.819881,0.938000},
      {-14.571006,-12.731894,1.444600,1.155655,2.393369,0.938000},
      {-27.749556,-18.025495,2.010033,0.663528,2.479596,0.938000},
      {-6.083697,-5.821953,1.712394,0.299414,3.116444,0.938000}
    };
    double sig1BichselParameters[numRapidityWindows][6] = {
      {-5.253659,-3.860611,2.722785,0.938063,1.831117,1.121600},
      {-0.337862,-3.795008,0.814606,1.042412,2.534717,0.643936},
      {3.657587,-2.147362,0.233230,1.517358,2.983071,0.523100},
      {-5.150241,-4.649925,1.944546,0.845950,2.229480,0.906125},
      {-2.614199,-4.245613,1.250001,0.923550,2.404442,0.784636},
      {-5.149159,-4.653147,1.939991,0.630485,2.466005,0.934273},
      {-5.147741,-4.637226,1.955622,0.121486,2.933515,0.945826},
      {-3.298422,-3.040149,2.581071,-0.053511,3.403947,1.222049},
      {-3.721689,-3.095312,2.857154,-0.654774,4.103774,1.309822},
      {-14.781475,-6.670506,5.016989,0.427622,2.424976,1.353251},
      {-5.133429,-4.657034,1.945658,0.492435,2.762268,0.992960},
      {-5.119586,-4.668256,1.951927,1.488254,2.388292,1.053244},
      {-5.128534,-4.660921,1.946557,-0.338899,6.280077,1.238813},
      {-10.012040,-8.059375,1.800773,1.406721,2.520363,1.163635},
      {-9.990270,-8.015409,1.803408,1.325243,2.521330,1.218257},
      {-5.125269,-4.662303,1.949754,0.908734,2.720948,1.186670},
      {-5.063713,-4.698326,1.975382,1.417511,3.683013,1.544924}
    };

    TF1* getLnBichselFunction(std::string a_functionName,double a_mass, bool IsNegativeMomentum=false){
      //H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195
      TF1* funct;
      if(IsNegativeMomentum) funct = new TF1(a_functionName.c_str(), "log(pow(sqrt([5]*[5] + x*x)/(-1.0*x),[3]) * ([0] - [1]*log([2] + pow([5]/(-1.0*x),[4]))) - [1])",-10.0,-0.5);
      else funct = new TF1(a_functionName.c_str(), "log(pow(sqrt([5]*[5] + x*x)/x,[3]) * ([0] - [1]*log([2] + pow([5]/x,[4]))) - [1])",0.0,10.0);
      //   ln(   (beta)^-D * (A - B ln(C + betagamma^-E)) - B   )
      /*funct->SetParNames("A","B","C","D","E","m");
	funct->SetParameter(0,-5.18614e+00);// A
	funct->SetParameter(1,-4.51905e+00);// B
	funct->SetParameter(2, 2.25999e+00); //C
	funct->SetParameter(3, 7.86756e-01); //D
	funct->SetParameter(4, 2.32399e+00);// E
	funct->SetParameter(5,a_mass); // m*/
      funct->SetParNames("A","B","C","D","E","m");
      funct->SetParameter(0,-5.14284e+00);// A
      funct->SetParameter(1,-4.65154e+00);// B
      funct->SetParameter(2, 1.94239e+00); //C
      funct->SetParameter(3,-9.90463e-02); //D
      funct->SetParameter(4, 3.19269e+00);// E
      funct->SetParameter(5,a_mass); // m


      funct->SetNpx(1000);
      return funct;
    }

    void initialize()
    {
      lowestMomWhichHasAFit[0]=0.130000;
      lowestMomUsedForBichselFit[0]=0.200000;

      mean_graph[0] = new TGraph(20,pointwise_momentum_0,pointwise_mean_0);
      sig1_graph[0] = new TGraph(20,pointwise_momentum_0,pointwise_1sig_0);
      lowestMomWhichHasAFit[1]=0.170000;
      lowestMomUsedForBichselFit[1]=0.240000;
    
      mean_graph[1] = new TGraph(24,pointwise_momentum_1,pointwise_mean_1);
      sig1_graph[1] = new TGraph(24,pointwise_momentum_1,pointwise_1sig_1);
      lowestMomWhichHasAFit[2]=0.250000;
      lowestMomUsedForBichselFit[2]=0.320000;
    
      mean_graph[2] = new TGraph(32,pointwise_momentum_2,pointwise_mean_2);
      sig1_graph[2] = new TGraph(32,pointwise_momentum_2,pointwise_1sig_2);
        
      // Loop over rapidity steps 
      for(int i=0; i<numRapidityWindows; i++)
	{
	  bichselFunction_mean[i] = getLnBichselFunction(Form("bichselFunction_mean_%d",i),0.938);
	  bichselFunction_1sig[i] = getLnBichselFunction(Form("bichselFunction_1sig_%d",i),0.938);
	  for(int j=0; j<6; j++){
	    bichselFunction_mean[i]->SetParameter(j,meanBichselParameters[i][j]);
	    bichselFunction_1sig[i]->SetParameter(j,sig1BichselParameters[i][j]);
	  }
	}
    } // End initialize()


    double getNewNSigmaProton(double rapidity, double momentum, double dedx)
    {
      // Later when the rapidity, momentum, and dedx of a track has been grabbed, calculate the nσ_p value  
      double nSigmaProton = 999; 
      double lndedx = log(dedx); 
      rapidity = fabs(rapidity); // take absolute value of rapidity
      //int rapIndex = TMath::FloorNint(10.0*rapidity);
      int rapIndex_lo = TMath::FloorNint(10.0*(-0.05+TMath::Abs(rapidity)));
      int rapIndex_hi = rapIndex_lo + 1;
      if(rapIndex_lo==-1) rapIndex_lo=0;
      if(rapIndex_hi==-1) rapIndex_hi=0;
      if(rapIndex_hi>=17) rapIndex_hi=16;
      if(rapIndex_lo>=17) rapIndex_lo=16;
      //double rapPointwiseCutoff = 10.0*numRapidityWindowsToUsePointwiseCorrection;
      double shiftValue_lo = 0.0;
      double shiftValue_hi = 0.0;
      double shiftValue = 0.0;
      double stretchValue_lo = 1.0;
      double stretchValue_hi = 1.0;
      double stretchValue = 1.0;
      if(rapIndex_lo < numRapidityWindowsToUsePointwiseCorrection && momentum < lowestMomUsedForBichselFit[rapIndex_lo])
	{
	  shiftValue_lo = mean_graph[rapIndex_lo]->Eval(momentum); 
	  stretchValue_lo = sig1_graph[rapIndex_lo]->Eval(momentum); 
	  if(momentum < lowestMomWhichHasAFit[rapIndex_lo])
	    {
	      shiftValue_lo = mean_graph[rapIndex_lo]->GetY()[0]; 
	      stretchValue_lo = sig1_graph[rapIndex_lo]->GetY()[0];
	    } 
	}
      else 
	{ 
	  shiftValue_lo = bichselFunction_mean[rapIndex_lo]->Eval(momentum); 
	  stretchValue_lo = bichselFunction_1sig[rapIndex_lo]->Eval(momentum); 
	}
      if(rapIndex_hi < numRapidityWindowsToUsePointwiseCorrection && momentum < lowestMomUsedForBichselFit[rapIndex_hi])
	{
	  shiftValue_hi = mean_graph[rapIndex_hi]->Eval(momentum); 
	  stretchValue_hi = sig1_graph[rapIndex_hi]->Eval(momentum);
	  if(momentum < lowestMomWhichHasAFit[rapIndex_hi])
	    {
	      shiftValue_hi = mean_graph[rapIndex_hi]->GetY()[0]; 
	      stretchValue_hi = sig1_graph[rapIndex_hi]->GetY()[0];
	    } 
	}
      else 
	{ shiftValue_hi = bichselFunction_mean[rapIndex_hi]->Eval(momentum); 
	  stretchValue_hi = bichselFunction_1sig[rapIndex_hi]->Eval(momentum); 
	}
     
      double rapidity_lo = 0.05+0.1*(rapIndex_lo);
      double rapidity_hi = 0.05+0.1*(rapIndex_hi);
      double relativeWeight;
      if(rapidity_lo == rapidity_hi) relativeWeight = 0.5;

      if(rapidity < rapidity_lo) relativeWeight = 0.0;
      else if(rapidity > rapidity_hi) relativeWeight = 1.0;
      else relativeWeight = (rapidity-rapidity_lo)/(rapidity_hi-rapidity_lo);

      shiftValue = relativeWeight*shiftValue_hi + (1.0-relativeWeight)*shiftValue_lo;
      stretchValue = relativeWeight*stretchValue_hi + (1.0-relativeWeight)*stretchValue_lo;
     
      // Now shift and stretch ln(de/dx) to get nSigmaProton 
      double shiftedLnDedx = lndedx-shiftValue;
      nSigmaProton = shiftedLnDedx/(shiftValue-stretchValue);
     
      // Reject tracks with the following traits 
      //if(nHitsDedx<20) nSigmaProton = 999; 
      if(momentum<0.0 || momentum > 10.0) nSigmaProton = 999; 

      return nSigmaProton;
    } // End getNewNSigmaProton()
  }; // End struct NewNSigmaProton3p22GeV


  
  ////////
  //   Calculates the event plane angle in every subevent.
  ////////
  void getAllPsi(Event &eventInfo, Double_t order_m)
  {
    eventInfo.psiTpc  = TMath::ATan2(eventInfo.YnTpc,  eventInfo.XnTpc)  / order_m;
    eventInfo.psiTpcA = TMath::ATan2(eventInfo.YnTpcA, eventInfo.XnTpcA) / order_m;
    eventInfo.psiTpcB = TMath::ATan2(eventInfo.YnTpcB, eventInfo.XnTpcB) / order_m;
    eventInfo.psiEpd  = TMath::ATan2(eventInfo.YnEpd,  eventInfo.XnEpd)  / order_m;
    eventInfo.psiEpdA = TMath::ATan2(eventInfo.YnEpdA, eventInfo.XnEpdA) / order_m;
    eventInfo.psiEpdB = TMath::ATan2(eventInfo.YnEpdB, eventInfo.XnEpdB) / order_m;
  }

  ////////
  //   Moves event plane angles back into the -pi to pi period.
  ////////
  Double_t angleShift(Double_t angle, Int_t order)
  {
    if (angle < -TMath::Pi()/(Double_t)order) { angle += TMath::TwoPi()/(Double_t)order; }
    else if (angle >  TMath::Pi()/(Double_t)order) { angle -= TMath::TwoPi()/(Double_t)order; }
    return angle;
  }

  ////////
  //   Moves all event plane angles into the -pi to pi period.
  ////////
  void setAllPeriods(Event &eventInfo, Double_t order_m)
  {
    eventInfo.psiTpc  = angleShift(eventInfo.psiTpc,  order_m);
    eventInfo.psiTpcA = angleShift(eventInfo.psiTpcA, order_m);
    eventInfo.psiTpcB = angleShift(eventInfo.psiTpcB, order_m);
    eventInfo.psiEpd  = angleShift(eventInfo.psiEpd,  order_m);
    eventInfo.psiEpdA = angleShift(eventInfo.psiEpdA, order_m);
    eventInfo.psiEpdB = angleShift(eventInfo.psiEpdB, order_m);
  }

  ////////
  //   Checks if an event has any flow vectors equal to zero. Updates the event's member variable "badEvent".
  ////////
  void checkZeroQ(Event &event)
  {
    //if (event.XnTpc == 0 && event.YnTpc == 0) { event.badEvent = true; }
    if (event.XnTpcA == 0 && event.YnTpcA == 0) { event.badEvent = true; }
    if (event.XnTpcB == 0 && event.YnTpcB == 0) { event.badEvent = true; }
    //else if (event.XnEpd  == 0 && event.YnEpd  == 0) { event.badEvent = true; }
    if (event.XnEpdA == 0 && event.YnEpdA == 0) { event.badEvent = true; }
    if (event.XnEpdB == 0 && event.YnEpdB == 0) { event.badEvent = true; }
  }


  ////////
  //   Using px, py, pz, and rest mass, return rapidity
  ////////
  Double_t rapidity(Double_t px, Double_t py, Double_t pz, Double_t mass)
  {
    Double_t rapidity, energy, momentum = 0;
    momentum = TMath::Sqrt(px*px + py*py + pz*pz);
    energy   = TMath::Sqrt(momentum*momentum + mass*mass);
    rapidity = TMath::ATanH(pz/energy);
    return rapidity;
  }

  Double_t pseudorapidity(Double_t px, Double_t py, Double_t pz)
  {
    Double_t momentum = TMath::Sqrt(px*px + py*py + pz*pz);
    return 0.5 * TMath::Log( (momentum + pz)/(momentum - pz) );
  }

  ////////
  //   Using px, py, pz, and rest mass, return transverse mass
  ////////
  Double_t transMass(Double_t px, Double_t py, Double_t mass) 
  { return TMath::Sqrt(mass*mass + px*px + py*py); }


  Double_t transMomentum(Double_t px, Double_t py)
  { return TMath::Sqrt(px*px + py*py); }

 
  Double_t totalMomentum(Double_t px, Double_t py, Double_t pz)
  {
    return TMath::Sqrt(px*px + py*py + pz*pz);
  }


  Double_t phi(Double_t px, Double_t py)
  { return TMath::ATan2(py, px); }

  Int_t epdSector(Short_t tileID)
  { return TMath::Abs( tileID/100 ); }


  Int_t epdRow(Short_t tileID)
  {
    Short_t absID = TMath::Abs(tileID);
    Short_t tileNum = absID % 100;
    Int_t rowNum = 0;
    
    if (tileNum == 1) rowNum = 1;
    else if (tileNum == 2 || tileNum == 3) rowNum = 2;
    else if (tileNum == 4 || tileNum == 5) rowNum = 3;
    else if (tileNum == 6 || tileNum == 7) rowNum = 4;
    else if (tileNum == 8 || tileNum == 9) rowNum = 5;
    else if (tileNum == 10 || tileNum == 11) rowNum = 6;
    else if (tileNum == 12 || tileNum == 13) rowNum = 7;
    else if (tileNum == 14 || tileNum == 15) rowNum = 8;
    else if (tileNum == 16 || tileNum == 17) rowNum = 9;
    else if (tileNum == 18 || tileNum == 19) rowNum = 10;
    else if (tileNum == 20 || tileNum == 21) rowNum = 11;
    else if (tileNum == 22 || tileNum == 23) rowNum = 12;
    else if (tileNum == 24 || tileNum == 25) rowNum = 13;
    else if (tileNum == 26 || tileNum == 27) rowNum = 14;
    else if (tileNum == 28 || tileNum == 29) rowNum = 15;
    else if (tileNum == 30 || tileNum == 31) rowNum = 16;

    return rowNum;
  }

  ////////
  //   Using px, py, pz, and rest mass, fill histograms raw of dN/dy,
  // and dN/dmT (shifted left by m0), and a 2D histogram of mT-m0 vs y.
  ////////
  void fillRawSpect(Double_t px, Double_t py, Double_t pz, Double_t mass, TH1D *dndy, TH1D *dndm, TH2D *MvsY)
  {
    Double_t y  = rapidity(px, py, pz, mass);
    Double_t mT = transMass(px, py, mass);
    Double_t M  = mT - mass;
    dndy->Fill(y);
    dndm->Fill(M);
    MvsY->Fill(y,M, 1/(TMath::TwoPi() * mT));
  }


  Double_t getTpcEff(Double_t yCM, Double_t pT, TH2D *h2_ratio)
  {
    Int_t xBin = h2_ratio->GetXaxis()->FindBin(yCM);
    Int_t yBin = h2_ratio->GetYaxis()->FindBin(pT);
    Double_t efficiency = h2_ratio->GetBinContent(xBin, yBin);
    return (efficiency == 0.0 || efficiency > 1.3) ? -1.0 : efficiency;
  }


  Double_t getTofEff(Double_t eta, Double_t pT, TH2D *h2_ratio)
  {
    Int_t xBin = h2_ratio->GetXaxis()->FindBin(eta);
    Int_t yBin = h2_ratio->GetYaxis()->FindBin(pT);
    Double_t efficiency = h2_ratio->GetBinContent(xBin, yBin);
    return (efficiency == 0.0) ? -1.0 : efficiency;
  }

  Double_t getTPCv1Weight(Double_t centID, Double_t eta, TProfile2D* p2_TPCv1Weights)
  {
    Int_t xBin = p2_TPCv1Weights->GetXaxis()->FindBin(centID);
    Int_t yBin = p2_TPCv1Weights->GetYaxis()->FindBin(eta);
    return p2_TPCv1Weights->GetBinContent(xBin, yBin);
  }

  Double_t getEPDv1Weight(Double_t centID, Double_t ring, TProfile2D* p2_EPDv1Weights)
  {
    Int_t xBin = p2_EPDv1Weights->GetXaxis()->FindBin(centID);
    Int_t yBin = p2_EPDv1Weights->GetYaxis()->FindBin(ring);
    return p2_EPDv1Weights->GetBinContent(xBin, yBin);
  }


  ////////
  //   Recenters the flow vectors of every subevent region in an event using the averages 
  //  found over all events in the full dataset and then recalculates event plane angles.
  ////////
  void recenterQ(Event &eventInfo, TFile *correctionInputFile, Double_t order_m)
  {
    TH1D *h_XnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_XnTpc");
    TH1D *h_XnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcA");
    TH1D *h_XnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcB");
    TH1D *h_XnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_XnEpd");
    TH1D *h_XnEpdA_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdA");
    TH1D *h_XnEpdB_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdB");

    TH1D *h_YnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_YnTpc");
    TH1D *h_YnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcA");
    TH1D *h_YnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcB");
    TH1D *h_YnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_YnEpd");
    TH1D *h_YnEpdA_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdA");
    TH1D *h_YnEpdB_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdB");

    Double_t d_XnTpc_Avg  = h_XnTpc_INPUT->GetMean();
    Double_t d_XnTpcA_Avg = h_XnTpcA_INPUT->GetMean();
    Double_t d_XnTpcB_Avg = h_XnTpcB_INPUT->GetMean();
    Double_t d_XnEpd_Avg  = h_XnEpd_INPUT->GetMean();
    Double_t d_XnEpdA_Avg = h_XnEpdA_INPUT->GetMean();
    Double_t d_XnEpdB_Avg = h_XnEpdB_INPUT->GetMean();

    Double_t d_YnTpc_Avg  = h_YnTpc_INPUT->GetMean();
    Double_t d_YnTpcA_Avg = h_YnTpcA_INPUT->GetMean();
    Double_t d_YnTpcB_Avg = h_YnTpcB_INPUT->GetMean();
    Double_t d_YnEpd_Avg  = h_YnEpd_INPUT->GetMean();
    Double_t d_YnEpdA_Avg = h_YnEpdA_INPUT->GetMean();
    Double_t d_YnEpdB_Avg = h_YnEpdB_INPUT->GetMean();


    eventInfo.XnTpc  -= d_XnTpc_Avg;
    eventInfo.XnTpcA -= d_XnTpcA_Avg;
    eventInfo.XnTpcB -= d_XnTpcB_Avg;
    eventInfo.XnEpd  -= d_XnEpd_Avg;
    eventInfo.XnEpdA -= d_XnEpdA_Avg;
    eventInfo.XnEpdB -= d_XnEpdB_Avg;

    eventInfo.YnTpc  -= d_YnTpc_Avg;
    eventInfo.YnTpcA -= d_YnTpcA_Avg;
    eventInfo.YnTpcB -= d_YnTpcB_Avg;
    eventInfo.YnEpd  -= d_YnEpd_Avg;
    eventInfo.YnEpdA -= d_YnEpdA_Avg;
    eventInfo.YnEpdB -= d_YnEpdB_Avg;

    checkZeroQ(eventInfo);

    getAllPsi(eventInfo, order_m);
    setAllPeriods(eventInfo, order_m);
  }// End recenterQ()


  ////////
  //   Performs the event-by-event shifting described in the Poskanzer paper to flatten the event 
  //  plane angle distributions of each subevent.
  ////////
  void shiftPsi(Event &eventInfo, TFile *correctionInputFile, Double_t order_m, Int_t shiftTerms)
  {
    TProfile *p_sinAvgsTpc_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsTpc");
    TProfile *p_cosAvgsTpc_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsTpc");
    TProfile *p_sinAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcA");
    TProfile *p_cosAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcA");
    TProfile *p_sinAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcB");
    TProfile *p_cosAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcB");
    TProfile *p_sinAvgsEpd_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsEpd");
    TProfile *p_cosAvgsEpd_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsEpd");
    TProfile *p_sinAvgsEpdA_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdA");
    TProfile *p_cosAvgsEpdA_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdA");
    TProfile *p_sinAvgsEpdB_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdB");
    TProfile *p_cosAvgsEpdB_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdB");


    // Get corrected event plane angles //

    Double_t psiTpc_delta  = 0;
    Double_t psiTpcA_delta = 0;
    Double_t psiTpcB_delta = 0;
    Double_t psiEpd_delta  = 0;
    Double_t psiEpdA_delta = 0;
    Double_t psiEpdB_delta = 0;

    Double_t jthSinAvg_Tpc  = 0;
    Double_t jthCosAvg_Tpc  = 0;
    Double_t jthSinAvg_TpcA = 0;
    Double_t jthCosAvg_TpcA = 0;
    Double_t jthSinAvg_TpcB = 0;
    Double_t jthCosAvg_TpcB = 0;
    Double_t jthSinAvg_Epd  = 0;
    Double_t jthCosAvg_Epd  = 0;
    Double_t jthSinAvg_EpdA = 0;
    Double_t jthCosAvg_EpdA = 0;
    Double_t jthSinAvg_EpdB = 0;
    Double_t jthCosAvg_EpdB = 0;


    for (Int_t j = 1; j <= shiftTerms; j++)    // Build the correction sums
      {
	jthSinAvg_Tpc  = p_sinAvgsTpc_INPUT->GetBinContent(j);
	jthCosAvg_Tpc  = p_cosAvgsTpc_INPUT->GetBinContent(j);
	jthSinAvg_TpcA = p_sinAvgsTpcA_INPUT->GetBinContent(j);
	jthCosAvg_TpcA = p_cosAvgsTpcA_INPUT->GetBinContent(j);
	jthSinAvg_TpcB = p_sinAvgsTpcB_INPUT->GetBinContent(j);
	jthCosAvg_TpcB = p_cosAvgsTpcB_INPUT->GetBinContent(j);
	jthSinAvg_Epd  = p_sinAvgsEpd_INPUT->GetBinContent(j);
	jthCosAvg_Epd  = p_cosAvgsEpd_INPUT->GetBinContent(j);
	jthSinAvg_EpdA = p_sinAvgsEpdA_INPUT->GetBinContent(j);
	jthCosAvg_EpdA = p_cosAvgsEpdA_INPUT->GetBinContent(j);
	jthSinAvg_EpdB = p_sinAvgsEpdB_INPUT->GetBinContent(j);
	jthCosAvg_EpdB = p_cosAvgsEpdB_INPUT->GetBinContent(j);

	psiTpc_delta  += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_Tpc * TMath::Cos((Double_t)j * order_m * eventInfo.psiTpc) 
							+jthCosAvg_Tpc * TMath::Sin((Double_t)j * order_m * eventInfo.psiTpc));
	psiTpcA_delta += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_TpcA * TMath::Cos((Double_t)j * order_m * eventInfo.psiTpcA) 
							+jthCosAvg_TpcA * TMath::Sin((Double_t)j * order_m * eventInfo.psiTpcA));
	psiTpcB_delta += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_TpcB * TMath::Cos((Double_t)j * order_m * eventInfo.psiTpcB) 
							+jthCosAvg_TpcB * TMath::Sin((Double_t)j * order_m * eventInfo.psiTpcB));
	psiEpd_delta  += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_Epd * TMath::Cos((Double_t)j * order_m * eventInfo.psiEpd)
							+jthCosAvg_Epd * TMath::Sin((Double_t)j * order_m * eventInfo.psiEpd));
	psiEpdA_delta += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_EpdA * TMath::Cos((Double_t)j * order_m * eventInfo.psiEpdA)
							+jthCosAvg_EpdA * TMath::Sin((Double_t)j * order_m * eventInfo.psiEpdA));
	psiEpdB_delta += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_EpdB * TMath::Cos((Double_t)j * order_m * eventInfo.psiEpdB)
							+jthCosAvg_EpdB * TMath::Sin((Double_t)j * order_m * eventInfo.psiEpdB));
      }

    // Shift event plane angles
    eventInfo.psiTpc  += psiTpc_delta;
    eventInfo.psiTpcA += psiTpcA_delta;
    eventInfo.psiTpcB += psiTpcB_delta;
    eventInfo.psiEpd  += psiEpd_delta;
    eventInfo.psiEpdA += psiEpdA_delta;
    eventInfo.psiEpdB += psiEpdB_delta;

    // Keep angles in the correct period
    setAllPeriods(eventInfo, order_m);
  }// End shiftPsi()


  Bool_t momDepDeuteronID(Double_t sqrt_s_NN, Double_t momentum, Double_t zDeuteron, 
			  Bool_t tofTrack, Double_t m2, 
			  Double_t zLowBoundWithTof, Double_t zHighBoundWithTof, Double_t m2LowBoundWithTof, Double_t m2HighBoundWithTof)
  {
    Bool_t deuteron = false;

    if (sqrt_s_NN == 3.0)
      {
	if (momentum >= 0.4 && momentum < 3.0)
	  {
	    if      (momentum >= 0.4 && momentum < 0.5 && zDeuteron > -0.476112 && zDeuteron < 0.248539) deuteron = true;
	    else if (momentum >= 0.5 && momentum < 0.6 && zDeuteron > -0.445644 && zDeuteron < 0.311067) deuteron = true;
	    else if (momentum >= 0.6 && momentum < 0.7 && zDeuteron > -0.43008  && zDeuteron < 0.331624) deuteron = true;
	    else if (momentum >= 0.7 && momentum < 0.8 && zDeuteron > -0.416061 && zDeuteron < 0.341399) deuteron = true;
	    else if (momentum >= 0.8 && momentum < 0.9 && zDeuteron > -0.404842 && zDeuteron < 0.338091) deuteron = true;
	    else if (momentum >= 0.9 && momentum < 1.0 && zDeuteron > -0.37419  && zDeuteron < 0.337724) deuteron = true;
	    else if (momentum >= 1.0 && momentum < 1.1 && zDeuteron > -0.32986  && zDeuteron < 0.332241) deuteron = true;
	    else if (momentum >= 1.1 && momentum < 1.2 && zDeuteron > -0.332995 && zDeuteron < 0.325582) deuteron = true;
	    else if (momentum >= 1.2 && momentum < 1.3 && zDeuteron > -0.306145 && zDeuteron < 0.319532) deuteron = true;
	    else if (momentum >= 1.3 && momentum < 1.4 && zDeuteron > -0.275987 && zDeuteron < 0.313227) deuteron = true;
	    else if (momentum >= 1.4 && momentum < 1.5 && zDeuteron > -0.250464 && zDeuteron < 0.301911) deuteron = true;
	    else if (momentum >= 1.5 && momentum < 1.6 && zDeuteron > -0.215135 && zDeuteron < 0.302149) deuteron = true;
	    else if (momentum >= 1.6 && momentum < 1.7 && zDeuteron > -0.176733 && zDeuteron < 0.308644) deuteron = true;
	    else if (momentum >= 1.7 && momentum < 1.8 && zDeuteron > -0.160866 && zDeuteron < 0.29673) deuteron = true;
	    else if (momentum >= 1.8 && momentum < 1.9 && zDeuteron > -0.149249 && zDeuteron < 0.281362) deuteron = true;
	    else if (momentum >= 1.9 && momentum < 2.0 && zDeuteron > -0.0830817 && zDeuteron < 0.273483) deuteron = true;
	    else if (momentum >= 2.0 && momentum < 2.1 && zDeuteron > -0.065219 && zDeuteron < 0.269654) deuteron = true;
	    else if (momentum >= 2.1 && momentum < 2.2 && zDeuteron > -0.04952  && zDeuteron < 0.265074) deuteron = true;
	    else if (momentum >= 2.2 && momentum < 2.3 && zDeuteron > -0.0358834 && zDeuteron < 0.258749) deuteron = true;
	    else if (momentum >= 2.3 && momentum < 2.4 && zDeuteron > -0.0218641 && zDeuteron < 0.25294) deuteron = true;
	    else if (momentum >= 2.4 && momentum < 2.5 && zDeuteron > -0.0114193 && zDeuteron < 0.244108) deuteron = true;
	    else if (momentum >= 2.5 && momentum < 2.6 && zDeuteron > -0.000659632 && zDeuteron < 0.205416) deuteron = true;
	    else if (momentum >= 2.6 && momentum < 2.7 && zDeuteron > 0.010662  && zDeuteron < 0.198006) deuteron = true;
	    else if (momentum >= 2.7 && momentum < 2.8 && zDeuteron > 0.0203815 && zDeuteron < 0.189092) deuteron = true;
	    else if (momentum >= 2.8 && momentum < 2.9 && zDeuteron > 0.0313737 && zDeuteron < 0.181285) deuteron = true;
	    else if (momentum >= 2.9 && momentum < 3.0 && zDeuteron > 0.0446902 && zDeuteron < 0.174561) deuteron = true;
	  }
	else if (tofTrack)
	  {
	    if (zDeuteron > zLowBoundWithTof && zDeuteron < zHighBoundWithTof &&
		m2 > m2LowBoundWithTof && m2 < m2HighBoundWithTof)
	      deuteron = true;
	  }
      }// End if (3.0 GeV)

    return deuteron;
  }// End momDepDeuteronID()



  Bool_t momDepDeuteronID_lowSystematics(Double_t sqrt_s_NN, Double_t momentum, Double_t zDeuteron, 
					 Bool_t tofTrack, Double_t m2, 
					 Double_t zLowBoundWithTof, Double_t zHighBoundWithTof, Double_t m2LowBoundWithTof, Double_t m2HighBoundWithTof)
  {
    Bool_t deuteron = false;

    if (sqrt_s_NN == 3.0)
      {
	if (momentum >= 0.4 && momentum < 3.0)
	  {
	    if      (momentum >= 0.4 && momentum < 0.5 && zDeuteron > -0.4036469 && zDeuteron < 0.1760739) deuteron = true;
	    else if (momentum >= 0.5 && momentum < 0.6 && zDeuteron > -0.3699729 && zDeuteron < 0.2353959) deuteron = true;
	    else if (momentum >= 0.6 && momentum < 0.7 && zDeuteron > -0.3539096 && zDeuteron < 0.2554536) deuteron = true;
	    else if (momentum >= 0.7 && momentum < 0.8 && zDeuteron > -0.340315 && zDeuteron < 0.265653) deuteron = true;
	    else if (momentum >= 0.8 && momentum < 0.9 && zDeuteron > -0.3305487 && zDeuteron < 0.2637977) deuteron = true;
	    else if (momentum >= 0.9 && momentum < 1.0 && zDeuteron > -0.3029986  && zDeuteron < 0.2665326) deuteron = true;
	    else if (momentum >= 1.0 && momentum < 1.1 && zDeuteron > -0.2636499  && zDeuteron < 0.2660309) deuteron = true;
	    else if (momentum >= 1.1 && momentum < 1.2 && zDeuteron > -0.2671373 && zDeuteron < 0.2597243) deuteron = true;
	    else if (momentum >= 1.2 && momentum < 1.3 && zDeuteron > -0.2435773 && zDeuteron < 0.2569643) deuteron = true;
	    else if (momentum >= 1.3 && momentum < 1.4 && zDeuteron > -0.2170656 && zDeuteron < 0.2543056) deuteron = true;
	    else if (momentum >= 1.4 && momentum < 1.5 && zDeuteron > -0.1952265 && zDeuteron < 0.2466735) deuteron = true;
	    else if (momentum >= 1.5 && momentum < 1.6 && zDeuteron > -0.1634066 && zDeuteron < 0.2504206) deuteron = true;
	    else if (momentum >= 1.6 && momentum < 1.7 && zDeuteron > -0.1281953 && zDeuteron < 0.2601063) deuteron = true;
	    else if (momentum >= 1.7 && momentum < 1.8 && zDeuteron > -0.1151064 && zDeuteron < 0.2509704) deuteron = true;
	    else if (momentum >= 1.8 && momentum < 1.9 && zDeuteron > -0.1061879 && zDeuteron < 0.2383009) deuteron = true;
	    else if (momentum >= 1.9 && momentum < 2.0 && zDeuteron > -0.04742523 && zDeuteron < 0.23782653) deuteron = true;
	    else if (momentum >= 2.0 && momentum < 2.1 && zDeuteron > -0.0317317 && zDeuteron < 0.2361667) deuteron = true;
	    else if (momentum >= 2.1 && momentum < 2.2 && zDeuteron > -0.0180606  && zDeuteron < 0.2336146) deuteron = true;
	    else if (momentum >= 2.2 && momentum < 2.3 && zDeuteron > -0.00642016 && zDeuteron < 0.22928576) deuteron = true;
	    else if (momentum >= 2.3 && momentum < 2.4 && zDeuteron > 0.00561631 && zDeuteron < 0.22545959) deuteron = true;
	    else if (momentum >= 2.4 && momentum < 2.5 && zDeuteron > 0.01413343 && zDeuteron < 0.21855527) deuteron = true;
	    else if (momentum >= 2.5 && momentum < 2.6 && zDeuteron > 0.019947933 && zDeuteron < 0.184808437) deuteron = true;
	    else if (momentum >= 2.6 && momentum < 2.7 && zDeuteron > 0.0293964  && zDeuteron < 0.1792716) deuteron = true;
	    else if (momentum >= 2.7 && momentum < 2.8 && zDeuteron > 0.03725255 && zDeuteron < 0.17222095) deuteron = true;
	    else if (momentum >= 2.8 && momentum < 2.9 && zDeuteron > 0.04636483 && zDeuteron < 0.16629387) deuteron = true;
	    else if (momentum >= 2.9 && momentum < 3.0 && zDeuteron > 0.05767728 && zDeuteron < 0.16157392) deuteron = true;
	  }
	else if (tofTrack)
	  {
	    if (zDeuteron > zLowBoundWithTof && zDeuteron < zHighBoundWithTof &&
		m2 > m2LowBoundWithTof && m2 < m2HighBoundWithTof)
	      deuteron = true;
	  }
      }// End if (3.0 GeV)

    return deuteron;
  }// End momDepDeuteronID_lowSystematics()



  Bool_t momDepDeuteronID_highSystematics(Double_t sqrt_s_NN, Double_t momentum, Double_t zDeuteron, 
					  Bool_t tofTrack, Double_t m2, 
					  Double_t zLowBoundWithTof, Double_t zHighBoundWithTof, Double_t m2LowBoundWithTof, Double_t m2HighBoundWithTof)
  {
    Bool_t deuteron = false;

    if (sqrt_s_NN == 3.0)
      {
	if (momentum >= 0.4 && momentum < 3.0)
	  {
	    if      (momentum >= 0.4 && momentum < 0.5 && zDeuteron > -0.5485771 && zDeuteron < 0.3210041) deuteron = true;
	    else if (momentum >= 0.5 && momentum < 0.6 && zDeuteron > -0.5213151 && zDeuteron < 0.3867381) deuteron = true;
	    else if (momentum >= 0.6 && momentum < 0.7 && zDeuteron > -0.5062504 && zDeuteron < 0.4077944) deuteron = true;
	    else if (momentum >= 0.7 && momentum < 0.8 && zDeuteron > -0.491807 && zDeuteron < 0.417145) deuteron = true;
	    else if (momentum >= 0.8 && momentum < 0.9 && zDeuteron > -0.4791353 && zDeuteron < 0.4123843) deuteron = true;
	    else if (momentum >= 0.9 && momentum < 1.0 && zDeuteron > -0.4453814  && zDeuteron < 0.4089154) deuteron = true;
	    else if (momentum >= 1.0 && momentum < 1.1 && zDeuteron > -0.3960701  && zDeuteron < 0.3984511) deuteron = true;
	    else if (momentum >= 1.1 && momentum < 1.2 && zDeuteron > -0.3988527 && zDeuteron < 0.3914397) deuteron = true;
	    else if (momentum >= 1.2 && momentum < 1.3 && zDeuteron > -0.3687127 && zDeuteron < 0.3820997) deuteron = true;
	    else if (momentum >= 1.3 && momentum < 1.4 && zDeuteron > -0.3349084 && zDeuteron < 0.3721484) deuteron = true;
	    else if (momentum >= 1.4 && momentum < 1.5 && zDeuteron > -0.3057015 && zDeuteron < 0.3571485) deuteron = true;
	    else if (momentum >= 1.5 && momentum < 1.6 && zDeuteron > -0.2668634 && zDeuteron < 0.3538774) deuteron = true;
	    else if (momentum >= 1.6 && momentum < 1.7 && zDeuteron > -0.2252707 && zDeuteron < 0.3571817) deuteron = true;
	    else if (momentum >= 1.7 && momentum < 1.8 && zDeuteron > -0.2066256 && zDeuteron < 0.3424896) deuteron = true;
	    else if (momentum >= 1.8 && momentum < 1.9 && zDeuteron > -0.1923101 && zDeuteron < 0.3244231) deuteron = true;
	    else if (momentum >= 1.9 && momentum < 2.0 && zDeuteron > -0.11873817 && zDeuteron < 0.30913947) deuteron = true;
	    else if (momentum >= 2.0 && momentum < 2.1 && zDeuteron > -0.0987063 && zDeuteron < 0.3031413) deuteron = true;
	    else if (momentum >= 2.1 && momentum < 2.2 && zDeuteron > -0.0809794  && zDeuteron < 0.2965334) deuteron = true;
	    else if (momentum >= 2.2 && momentum < 2.3 && zDeuteron > -0.06534664 && zDeuteron < 0.28821224) deuteron = true;
	    else if (momentum >= 2.3 && momentum < 2.4 && zDeuteron > -0.04934451 && zDeuteron < 0.28042041) deuteron = true;
	    else if (momentum >= 2.4 && momentum < 2.5 && zDeuteron > -0.03697203 && zDeuteron < 0.26966073) deuteron = true;
	    else if (momentum >= 2.5 && momentum < 2.6 && zDeuteron > -0.021267193 && zDeuteron < 0.226023563) deuteron = true;
	    else if (momentum >= 2.6 && momentum < 2.7 && zDeuteron > -0.0080724  && zDeuteron < 0.2167404) deuteron = true;
	    else if (momentum >= 2.7 && momentum < 2.8 && zDeuteron > 0.00351045 && zDeuteron < 0.20596305) deuteron = true;
	    else if (momentum >= 2.8 && momentum < 2.9 && zDeuteron > 0.01638257 && zDeuteron < 0.19627613) deuteron = true;
	    else if (momentum >= 2.9 && momentum < 3.0 && zDeuteron > 0.03170312 && zDeuteron < 0.18754808) deuteron = true;
	  }
	else if (tofTrack)
	  {
	    if (zDeuteron > zLowBoundWithTof && zDeuteron < zHighBoundWithTof &&
		m2 > m2LowBoundWithTof && m2 < m2HighBoundWithTof)
	      deuteron = true;
	  }
      }// End if (3.0 GeV)

    return deuteron;
  }// End momDepDeuteronID_highSystematics()



  Bool_t momDepTritonID(Double_t sqrt_s_NN, Double_t momentum, Double_t zTriton, 
			Bool_t tofTrack, Double_t m2, 
			Double_t zLowBoundWithTof, Double_t zHighBoundWithTof, Double_t m2LowBoundWithTof, Double_t m2HighBoundWithTof)
  {
    Bool_t triton = false;

    if (sqrt_s_NN == 3.0)
      {
	if (momentum >= 1.0 && momentum < 4.0)
	  {
	    if      (momentum >= 1.0 && momentum < 1.1 && zTriton > -0.332011 && zTriton < 0.251103) triton = true;
	    else if (momentum >= 1.1 && momentum < 1.2 && zTriton > -0.310412 && zTriton < 0.296090) triton = true;
	    else if (momentum >= 1.2 && momentum < 1.3 && zTriton > -0.293322 && zTriton < 0.334467) triton = true;
	    else if (momentum >= 1.3 && momentum < 1.4 && zTriton > -0.270550 && zTriton < 0.373857) triton = true;
	    else if (momentum >= 1.4 && momentum < 1.5 && zTriton > -0.248412 && zTriton < 0.406237) triton = true;
	    else if (momentum >= 1.5 && momentum < 1.6 && zTriton > -0.228044 && zTriton < 0.333261) triton = true;
	    else if (momentum >= 1.6 && momentum < 1.7 && zTriton > -0.210093 && zTriton < 0.343588) triton = true;
	    else if (momentum >= 1.7 && momentum < 1.8 && zTriton > -0.190900 && zTriton < 0.332586) triton = true;
	    else if (momentum >= 1.8 && momentum < 1.9 && zTriton > -0.183153 && zTriton < 0.334197) triton = true;
	    else if (momentum >= 1.9 && momentum < 2.0 && zTriton > -0.166020 && zTriton < 0.323303) triton = true;
	    else if (momentum >= 2.0 && momentum < 2.1 && zTriton > -0.102334 && zTriton < 0.307724) triton = true;
	    else if (momentum >= 2.1 && momentum < 2.2 && zTriton > -0.091053 && zTriton < 0.294345) triton = true;
	    else if (momentum >= 2.2 && momentum < 2.3 && zTriton > -0.076457 && zTriton < 0.285978) triton = true;
	    else if (momentum >= 2.3 && momentum < 2.4 && zTriton > -0.055669 && zTriton < 0.253769) triton = true;
	    else if (momentum >= 2.4 && momentum < 2.5 && zTriton > -0.035848 && zTriton < 0.254487) triton = true;
	    else if (momentum >= 2.5 && momentum < 2.6 && zTriton > -0.027266 && zTriton < 0.249350) triton = true;
	    else if (momentum >= 2.6 && momentum < 2.7 && zTriton > -0.028152 && zTriton < 0.236713) triton = true;
	    else if (momentum >= 2.7 && momentum < 2.8 && zTriton > -0.027867 && zTriton < 0.227672) triton = true;
	    else if (momentum >= 2.8 && momentum < 2.9 && zTriton > -0.024675 && zTriton < 0.222215) triton = true;
	    else if (momentum >= 2.9 && momentum < 3.0 && zTriton > -0.019179 && zTriton < 0.227362) triton = true;
	    else if (momentum >= 3.0 && momentum < 3.1 && zTriton > -0.013267 && zTriton < 0.236052) triton = true;
	    else if (momentum >= 3.1 && momentum < 3.2 && zTriton > -0.007851 && zTriton < 0.246071) triton = true;
	    else if (momentum >= 3.2 && momentum < 3.3 && zTriton > -0.006311 && zTriton < 0.254907) triton = true;
	    else if (momentum >= 3.3 && momentum < 3.4 && zTriton > 0.019834  && zTriton < 0.244291) triton = true;
	    else if (momentum >= 3.4 && momentum < 3.5 && zTriton > 0.031221  && zTriton < 0.273652) triton = true;
	    else if (momentum >= 3.5 && momentum < 3.6 && zTriton > 0.068248  && zTriton < 0.257484) triton = true;
	    else if (momentum >= 3.6 && momentum < 3.7 && zTriton > 0.088804  && zTriton < 0.260799) triton = true;
	    else if (momentum >= 3.7 && momentum < 3.8 && zTriton > 0.091490  && zTriton < 0.271776) triton = true;
	    else if (momentum >= 3.8 && momentum < 3.9 && zTriton > 0.106161  && zTriton < 0.285652) triton = true;
	    else if (momentum >= 3.9 && momentum < 4.0 && zTriton > 0.103653  && zTriton < 0.299234) triton = true;
	  }
	else if (tofTrack)
	  {
	    if (zTriton > zLowBoundWithTof && zTriton < zHighBoundWithTof &&
		m2 > m2LowBoundWithTof && m2 < m2HighBoundWithTof)
	      triton = true;
	  }
      }// End if (3.0 GeV)

    return triton;
  }// End momDepTritonID()


  Bool_t momDepTritonID_lowSystematics(Double_t sqrt_s_NN, Double_t momentum, Double_t zTriton, 
				       Bool_t tofTrack, Double_t m2, 
				       Double_t zLowBoundWithTof, Double_t zHighBoundWithTof, Double_t m2LowBoundWithTof, Double_t m2HighBoundWithTof)
  {
    Bool_t triton = false;

    if (sqrt_s_NN == 3.0)
      {
	if (momentum >= 1.0 && momentum < 4.0)
	  {
	    if      (momentum >= 1.0 && momentum < 1.1 && zTriton > -0.2736996 && zTriton < 0.1927916) triton = true;
	    else if (momentum >= 1.1 && momentum < 1.2 && zTriton > -0.2497618 && zTriton < 0.2354398) triton = true;
	    else if (momentum >= 1.2 && momentum < 1.3 && zTriton > -0.2305431 && zTriton < 0.2716881) triton = true;
	    else if (momentum >= 1.3 && momentum < 1.4 && zTriton > -0.2061093 && zTriton < 0.3094163) triton = true;
	    else if (momentum >= 1.4 && momentum < 1.5 && zTriton > -0.1829471 && zTriton < 0.3407721) triton = true;
	    else if (momentum >= 1.5 && momentum < 1.6 && zTriton > -0.1719135 && zTriton < 0.2771305) triton = true;
	    else if (momentum >= 1.6 && momentum < 1.7 && zTriton > -0.1547249 && zTriton < 0.2882199) triton = true;
	    else if (momentum >= 1.7 && momentum < 1.8 && zTriton > -0.1385514 && zTriton < 0.2802374) triton = true;
	    else if (momentum >= 1.8 && momentum < 1.9 && zTriton > -0.131418 && zTriton < 0.282462) triton = true;
	    else if (momentum >= 1.9 && momentum < 2.0 && zTriton > -0.1170877 && zTriton < 0.2743707) triton = true;
	    else if (momentum >= 2.0 && momentum < 2.1 && zTriton > -0.0613282 && zTriton < 0.2667182) triton = true;
	    else if (momentum >= 2.1 && momentum < 2.2 && zTriton > -0.0525132 && zTriton < 0.2558052) triton = true;
	    else if (momentum >= 2.2 && momentum < 2.3 && zTriton > -0.0402135 && zTriton < 0.2497345) triton = true;
	    else if (momentum >= 2.3 && momentum < 2.4 && zTriton > -0.0247252 && zTriton < 0.2228252) triton = true;
	    else if (momentum >= 2.4 && momentum < 2.5 && zTriton > -0.0068145 && zTriton < 0.2254535) triton = true;
	    else if (momentum >= 2.5 && momentum < 2.6 && zTriton > 0.0003956 && zTriton < 0.2216884) triton = true;
	    else if (momentum >= 2.6 && momentum < 2.7 && zTriton > -0.0016655 && zTriton < 0.2102265) triton = true;
	    else if (momentum >= 2.7 && momentum < 2.8 && zTriton > -0.0023131 && zTriton < 0.2021181) triton = true;
	    else if (momentum >= 2.8 && momentum < 2.9 && zTriton > 0.000014 && zTriton < 0.197526) triton = true;
	    else if (momentum >= 2.9 && momentum < 3.0 && zTriton > 0.0054751 && zTriton < 0.2027079) triton = true;
	    else if (momentum >= 3.0 && momentum < 3.1 && zTriton > 0.0116649 && zTriton < 0.2111201) triton = true;
	    else if (momentum >= 3.1 && momentum < 3.2 && zTriton > 0.0175412 && zTriton < 0.2206788) triton = true;
	    else if (momentum >= 3.2 && momentum < 3.3 && zTriton > 0.0198108 && zTriton < 0.2287852) triton = true;
	    else if (momentum >= 3.3 && momentum < 3.4 && zTriton > 0.0422797  && zTriton < 0.2218453) triton = true;
	    else if (momentum >= 3.4 && momentum < 3.5 && zTriton > 0.0554641  && zTriton < 0.2494089) triton = true;
	    else if (momentum >= 3.5 && momentum < 3.6 && zTriton > 0.0871716  && zTriton < 0.2385604) triton = true;
	    else if (momentum >= 3.6 && momentum < 3.7 && zTriton > 0.1060035  && zTriton < 0.2435995) triton = true;
	    else if (momentum >= 3.7 && momentum < 3.8 && zTriton > 0.1095186  && zTriton < 0.2537474) triton = true;
	    else if (momentum >= 3.8 && momentum < 3.9 && zTriton > 0.1241101  && zTriton < 0.2677029) triton = true;
	    else if (momentum >= 3.9 && momentum < 4.0 && zTriton > 0.1232111  && zTriton < 0.2796759) triton = true;
	  } 
	else if (tofTrack)
	  {
	    if (zTriton > zLowBoundWithTof && zTriton < zHighBoundWithTof &&
		m2 > m2LowBoundWithTof && m2 < m2HighBoundWithTof)
	      triton = true;
	  }
      }// End if (3.0 GeV)

    return triton;
  }// End momDepTritonID_lowSystematics()



  Bool_t momDepTritonID_highSystematics(Double_t sqrt_s_NN, Double_t momentum, Double_t zTriton, 
				       Bool_t tofTrack, Double_t m2, 
				       Double_t zLowBoundWithTof, Double_t zHighBoundWithTof, Double_t m2LowBoundWithTof, Double_t m2HighBoundWithTof)
  {
    Bool_t triton = false;

    if (sqrt_s_NN == 3.0)
      {
	if (momentum >= 1.0 && momentum < 4.0)
	  {
	    if      (momentum >= 1.0 && momentum < 1.1 && zTriton > -0.3903224 && zTriton < 0.3094144) triton = true;
	    else if (momentum >= 1.1 && momentum < 1.2 && zTriton > -0.3710622 && zTriton < 0.3567402) triton = true;
	    else if (momentum >= 1.2 && momentum < 1.3 && zTriton > -0.3561009 && zTriton < 0.3972459) triton = true;
	    else if (momentum >= 1.3 && momentum < 1.4 && zTriton > -0.3349907 && zTriton < 0.4382977) triton = true;
	    else if (momentum >= 1.4 && momentum < 1.5 && zTriton > -0.3138769 && zTriton < 0.4717019) triton = true;
	    else if (momentum >= 1.5 && momentum < 1.6 && zTriton > -0.2841745 && zTriton < 0.3893915) triton = true;
	    else if (momentum >= 1.6 && momentum < 1.7 && zTriton > -0.2654611 && zTriton < 0.3989561) triton = true;
	    else if (momentum >= 1.7 && momentum < 1.8 && zTriton > -0.2432486 && zTriton < 0.3849346) triton = true;
	    else if (momentum >= 1.8 && momentum < 1.9 && zTriton > -0.234888 && zTriton < 0.385932) triton = true;
	    else if (momentum >= 1.9 && momentum < 2.0 && zTriton > -0.2149523 && zTriton < 0.3722353) triton = true;
	    else if (momentum >= 2.0 && momentum < 2.1 && zTriton > -0.1433398 && zTriton < 0.3487298) triton = true;
	    else if (momentum >= 2.1 && momentum < 2.2 && zTriton > -0.1295928 && zTriton < 0.3328848) triton = true;
	    else if (momentum >= 2.2 && momentum < 2.3 && zTriton > -0.1127005 && zTriton < 0.3222215) triton = true;
	    else if (momentum >= 2.3 && momentum < 2.4 && zTriton > -0.0866128 && zTriton < 0.2847128) triton = true;
	    else if (momentum >= 2.4 && momentum < 2.5 && zTriton > -0.0648815 && zTriton < 0.2835205) triton = true;
	    else if (momentum >= 2.5 && momentum < 2.6 && zTriton > -0.0549276 && zTriton < 0.2770116) triton = true;
	    else if (momentum >= 2.6 && momentum < 2.7 && zTriton > -0.0546385 && zTriton < 0.2631995) triton = true;
	    else if (momentum >= 2.7 && momentum < 2.8 && zTriton > -0.0534209 && zTriton < 0.2532259) triton = true;
	    else if (momentum >= 2.8 && momentum < 2.9 && zTriton > -0.049364 && zTriton < 0.246904) triton = true;
	    else if (momentum >= 2.9 && momentum < 3.0 && zTriton > -0.0438331 && zTriton < 0.2520161) triton = true;
	    else if (momentum >= 3.0 && momentum < 3.1 && zTriton > -0.0381989 && zTriton < 0.2609839) triton = true;
	    else if (momentum >= 3.1 && momentum < 3.2 && zTriton > -0.0332432 && zTriton < 0.2714632) triton = true;
	    else if (momentum >= 3.2 && momentum < 3.3 && zTriton > -0.0324328 && zTriton < 0.2810288) triton = true;
	    else if (momentum >= 3.3 && momentum < 3.4 && zTriton > -0.0026117  && zTriton < 0.2667367) triton = true;
	    else if (momentum >= 3.4 && momentum < 3.5 && zTriton > 0.0069779  && zTriton < 0.2978951) triton = true;
	    else if (momentum >= 3.5 && momentum < 3.6 && zTriton > 0.0493244  && zTriton < 0.2764076) triton = true;
	    else if (momentum >= 3.6 && momentum < 3.7 && zTriton > 0.0716045  && zTriton < 0.2779985) triton = true;
	    else if (momentum >= 3.7 && momentum < 3.8 && zTriton > 0.0734614  && zTriton < 0.2898046) triton = true;
	    else if (momentum >= 3.8 && momentum < 3.9 && zTriton > 0.0882119  && zTriton < 0.3036011) triton = true;
	    else if (momentum >= 3.9 && momentum < 4.0 && zTriton > 0.0840949  && zTriton < 0.3187921) triton = true;
	  }
	else if (tofTrack)
	  {
	    if (zTriton > zLowBoundWithTof && zTriton < zHighBoundWithTof &&
		m2 > m2LowBoundWithTof && m2 < m2HighBoundWithTof)
	      triton = true;
	  }
      }// End if (3.0 GeV)

    return triton;
  }// End momDepTritonID_highSystematics()


}// End namespace FlowUtils

#endif
