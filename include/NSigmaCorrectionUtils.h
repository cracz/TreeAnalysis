#ifndef NSIGMACORRECTIONUTILS
#define NSIGMACORRECTIONUTILS

namespace NSigmaCorrectionUtils
{
  struct NewNSigmaProton3p2GeV
  {
    const static  int numRapidityWindows = 17;
    const static int numRapidityWindowsToUsePointwiseCorrection = 3;
    TF1* bichselFunction_mean[numRapidityWindows];
    TF1* bichselFunction_1sig[numRapidityWindows];
  
    TGraph* mean_graph[3];
    TGraph* sig1_graph[3];  

    // Pointwise-correction parameters to be used at low momenta and low rapidities 
    double lowestMomWhichHasAFit[3];
    double lowestMomUsedForBichselFit[3];
    
    double pointwise_momentum_0[20] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195};
    double pointwise_mean_0[20] = {4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.03425, 3.99311, 3.94308, 3.88901, 3.83291, 3.77598};
    double pointwise_1sig_0[20] = {4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 4.01109, 3.97051, 3.92808, 3.88058, 3.83105, 3.77989, 3.72508};

    double pointwise_momentum_1[24] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235};
    double pointwise_mean_1[24] = {3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.84933, 3.7977, 3.74982, 3.70112, 3.65017, 3.59674, 3.54167};
    double pointwise_1sig_1[24] = {3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.77977, 3.72948, 3.68502, 3.63939, 3.59034, 3.53756, 3.48334};
    
    double pointwise_momentum_2[32] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235, 0.245, 0.255, 0.265, 0.275, 0.285, 0.295, 0.305, 0.315};
    double pointwise_mean_2[32] = {3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.40402, 3.35661, 3.30871, 3.2577, 3.20611, 3.15557, 3.107};
    double pointwise_1sig_2[32] = {3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.33402, 3.28661, 3.23871, 3.1877, 3.13764, 3.08949, 3.0435};
    
    // Parameters from Bichsel curve fits (H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195) 
    double meanBichselParameters[numRapidityWindows][6] = {
      {-0.227966,-1.947452,1.698877,1.116508,2.491739,0.938000},
      {-0.906363,-2.348223,1.743461,1.108302,2.359716,0.938000},
      {-0.614651,-2.241335,1.624326,1.144715,2.380232,0.938000},
      {-1.186437,-2.620453,1.622527,1.118374,2.333279,0.938000},
      {-5.089407,-4.698291,1.938589,0.831378,2.262067,0.938000},
      {-5.073031,-4.778549,1.877451,0.762094,2.338385,0.938000},
      {-5.046706,-4.836419,1.837368,0.512003,2.557916,0.938000},
      {-4.965796,-4.879391,1.818065,0.347630,3.087001,0.938000},
      {-13.505080,-8.529639,2.498983,0.507081,2.637780,0.938000},
      {-5.035899,-5.161085,1.645866,-0.991880,4.372225,0.938000},
      {-7.988699,-6.848821,1.785826,0.857071,2.614882,0.938000},
      {-4.097592,-6.227094,1.083402,4.447512,1.566225,0.938000},
      {-3.368052,-6.739509,0.907986,5.251134,1.669462,0.938000},
      {-3.701538,-6.701420,0.979808,1.051356,3.165626,0.938000},
      {-4.439844,-6.816681,1.067793,1.580762,2.548197,0.938000},
      {-4.895664,-5.324965,1.554491,-0.287449,3.132595,0.938000},
      {-5.978530,-5.690012,1.740206,0.189245,3.177314,0.938000}
    };
    double sig1BichselParameters[numRapidityWindows][6] = {
      {4.942490,-1.732692,0.103650,1.584292,3.443852,0.458243},
      {-2.048825,-4.388962,1.075924,0.960740,2.413125,0.688201},
      {-5.140414,-4.654629,1.951665,0.819828,2.264229,0.887815},
      {-5.141626,-4.653979,1.947455,0.774334,2.331612,0.901746},
      {-3.995552,-4.761255,1.473468,0.714432,2.502485,0.813790},
      {-5.142394,-4.655081,1.943287,0.508679,2.604774,0.931181},
      {-3.718692,-3.581791,2.149504,0.305906,2.891777,1.063736},
      {-71.403803,-18.809103,19.233168,0.185852,1.072852,5.512266},
      {-15.970397,-9.042774,2.903327,0.563681,2.554410,0.987476},
      {1.518479,-0.603197,2.093214,0.349728,4.765742,1.611703},
      {-29.532566,-19.070781,1.988991,0.513077,2.446239,0.666395},
      {-5.144622,-4.653428,1.939297,-0.821817,6.644838,1.134653},
      {-5.121391,-4.665033,1.949753,-0.272327,6.340915,1.243531},
      {-5.094469,-4.682484,1.962738,1.485384,3.028247,1.248443},
      {-9.973230,-7.819392,1.867982,1.023731,2.730247,1.173559},
      {-5.136450,-4.659384,1.948466,0.589825,2.644892,1.066994},
      {-5.074660,-4.691923,1.970780,1.243181,3.600409,1.463297}
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
  }; // End struct NewNSigmaProton3p2GeV





  
  struct NewNSigmaProton3p5GeV
  {
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
    double pointwise_mean_0[20] = {4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.07434, 4.04196, 3.99361, 3.94046, 3.88689, 3.83252, 3.77603};
    double pointwise_1sig_0[20] = {4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 4.01008, 3.97848, 3.92785, 3.8764, 3.8277, 3.77849, 3.72457};
    
    double pointwise_momentum_1[24] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235};
    double pointwise_mean_1[24] = {3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.85984, 3.80459, 3.75502, 3.70479, 3.65265, 3.59745, 3.54021};
    double pointwise_1sig_1[24] = {3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.78984, 3.73459, 3.68829, 3.64167, 3.59196, 3.53768, 3.48136};
    
    double pointwise_momentum_2[32] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235, 0.245, 0.255, 0.265, 0.275, 0.285, 0.295, 0.305, 0.315};
    double pointwise_mean_2[32] = {3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.41102, 3.36142, 3.31105, 3.25866, 3.20617, 3.15464, 3.10548};
    double pointwise_1sig_2[32] = {3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.34102, 3.29142, 3.24105, 3.18881, 3.13825, 3.08904, 3.04236};
    
    // Parameters from Bichsel curve fits (H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195) 
    double meanBichselParameters[numRapidityWindows][6] = {
      {-5.092943,-4.640572,1.908917,0.978064,1.837915,0.938000},
      {-0.659163,-2.209373,1.678841,1.154059,2.290688,0.938000},
      {-0.313533,-2.035206,1.584761,1.230861,2.285030,0.938000},
      {-0.893688,-2.405065,1.603568,1.236016,2.212285,0.938000},
      {-5.067733,-4.718756,1.890024,0.972786,2.090088,0.938000},
      {-5.060938,-4.835593,1.801288,1.015586,2.057320,0.938000},
      {-5.045866,-4.880064,1.782700,0.842702,2.226889,0.938000},
      {-5.033661,-4.890768,1.799017,0.157471,2.869778,0.938000},
      {-0.956699,-3.468721,1.071352,-0.170151,3.351006,0.938000},
      {-4.982008,-5.098390,1.678426,0.386340,2.865584,0.938000},
      {-4.930145,-5.464084,1.505620,1.832087,2.178729,0.938000},
      {-4.667087,-5.612072,1.401343,1.278266,2.723832,0.938000},
      {-6.853560,-8.485027,1.141610,-0.227246,5.167554,0.938000},
      {-2.368033,-5.187790,0.995863,3.689650,2.167992,0.938000},
      {-4.316409,-6.463155,1.103758,0.478797,3.220366,0.938000},
      {-13.437760,-10.140298,1.823939,0.662773,2.510908,0.938000},
      {-5.051530,-4.695740,1.979253,0.674798,2.785134,0.938000}
    };
    double sig1BichselParameters[numRapidityWindows][6] = {
      {-5.249671,-4.645090,1.904914,0.960111,1.854916,0.921248},
      {-5.194184,-4.641954,1.937924,0.930961,2.007013,0.901315},
      {-5.181766,-4.637340,1.930889,0.939909,2.021784,0.913086},
      {0.579304,-3.268654,0.672862,1.317425,2.430856,0.642275},
      {-5.165629,-4.645111,1.933057,0.893594,2.180889,0.928042},
      {-5.174665,-4.644717,1.924935,0.908074,2.154252,0.951100},
      {-5.152355,-4.651824,1.937027,0.345109,2.695883,0.945865},
      {-4.467375,-3.971009,2.153613,-1.133300,4.556280,1.072056},
      {-3.516175,-3.189363,2.495112,-1.299702,4.860644,1.292207},
      {-1.752874,-2.696696,1.860238,-0.513426,3.830402,1.253298},
      {-5.270542,-4.714190,1.958827,0.991582,2.427328,1.005969},
      {-5.123778,-4.664506,1.950183,1.220159,2.506098,1.031556},
      {-5.122731,-4.664056,1.949713,-0.056186,4.798412,1.158343},
      {-5.100615,-4.678819,1.960143,1.521345,2.944816,1.227547},
      {-5.110781,-4.672200,1.955830,1.486902,2.845572,1.288074},
      {-5.143146,-4.653472,1.944525,0.412436,2.622885,0.992467},
      {-5.080330,-4.688606,1.968640,1.380224,3.451587,1.434801}
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
  }; // End struct NewNSigmaProton3p5GeV


  struct NewNSigmaProton3p9GeV
  {
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
    double pointwise_mean_0[20] = {4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.07005, 4.03022, 3.98377, 3.93277, 3.88062, 3.82665, 3.77048};
    double pointwise_1sig_0[20] = {4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 4.0092, 3.96611, 3.91731, 3.86877, 3.82141, 3.77262, 3.71876};
    
    double pointwise_momentum_1[24] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235};
    double pointwise_mean_1[24] = {3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.85093, 3.79815, 3.74984, 3.70108, 3.64975, 3.59552, 3.53934};
    double pointwise_1sig_1[24] = {3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.78093, 3.72845, 3.68411, 3.63898, 3.58984, 3.53651, 3.48123};
    
double pointwise_momentum_2[32] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235, 0.245, 0.255, 0.265, 0.275, 0.285, 0.295, 0.305, 0.315};
    double pointwise_mean_2[32] = {3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.41134, 3.3626, 3.31301, 3.26118, 3.20877, 3.1577, 3.10853};
    double pointwise_1sig_2[32] = {3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.34134, 3.2926, 3.24301, 3.19157, 3.14087, 3.09196, 3.04508};
    
    // Parameters from Bichsel curve fits (H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195) 
    double meanBichselParameters[numRapidityWindows][6] = {
      {-5.177325,-4.565114,2.043163,0.926697,1.969583,0.938000},
      {-0.625984,-2.174822,1.711623,1.135181,2.374554,0.938000},
      {-0.467159,-2.130139,1.618613,1.174733,2.376978,0.938000},
      {-0.831247,-2.391516,1.597115,1.167913,2.353091,0.938000},
      {-5.093345,-4.689092,1.943219,0.832131,2.276891,0.938000},
      {-5.074470,-4.764597,1.891964,0.663462,2.460011,0.938000},
      {-5.038901,-4.829667,1.845495,0.133349,2.944921,0.938000},
      {-0.043268,-2.919098,0.963756,0.088509,3.300096,0.938000},
      {0.006582,-3.898368,0.733902,-2.145955,5.128702,0.938000},
      {-14.107090,-9.562289,2.176158,0.763134,2.597167,0.938000},
      {-12.014173,-9.346169,1.819754,1.009344,2.611237,0.938000},
      {-4.640004,-6.112600,1.269502,1.735888,2.751780,0.938000},
      {-5.806170,-7.747801,1.139622,1.930779,2.599439,0.938000},
      {-2.268550,-7.013657,0.778746,2.618637,2.811688,0.938000},
      {-2.280379,-7.061412,0.771718,3.320349,2.448393,0.938000},
      {-5.293035,-6.518209,1.284588,0.567845,3.170296,0.938000},
      {-4.786069,-5.157588,1.620347,1.073383,2.594362,0.938000}
    };
    double sig1BichselParameters[numRapidityWindows][6] = {
      {-2.199537,-4.551449,1.048161,0.971736,2.322050,0.670206},
      {-2.110362,-4.574579,1.037392,0.964539,2.397126,0.669572},
      {-5.157718,-4.646264,1.943121,0.872093,2.172467,0.893694},
      {0.499905,-3.367670,0.682757,1.219520,2.519937,0.645850},
      {-3.552719,-4.598840,1.400043,0.791708,2.475654,0.803871},
      {-5.145075,-4.653687,1.942072,0.529104,2.597866,0.930610},
      {-3.551871,-3.207324,2.514574,0.150560,3.145590,1.161430},
      {-3.790708,-3.559189,2.217708,-0.791251,4.182170,1.134552},
      {-3.573812,-3.217454,2.481291,-1.893540,5.855862,1.318876},
      {-15.401816,-12.081929,1.649789,0.844405,2.561858,0.779057},
      {-9.819051,-8.753436,1.546401,1.105451,2.607804,0.869509},
      {-7.657118,-7.020629,1.625261,1.379417,2.609855,0.974502},
      {-6.966518,-6.414624,1.685056,1.513236,2.722290,1.075895},
      {-7.293939,-6.579775,1.719777,1.834657,2.665309,1.231057},
      {-7.212817,-6.609524,1.686214,1.698213,2.833678,1.318012},
      {-6.809715,-5.983380,1.829462,1.612812,2.721876,1.364471},
      {-5.521483,-5.039877,1.904576,1.867710,2.590536,1.324373}
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
  }; // End struct NewNSigmaProton3p9GeV



  struct NewNSigmaProton4p5GeV
  {
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
    double pointwise_mean_0[20] = {4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.09365, 4.05405, 4.00152, 3.94765, 3.89523, 3.841, 3.78386};
    double pointwise_1sig_0[20] = {4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 4.03159, 3.99234, 3.93627, 3.88369, 3.83597, 3.78675, 3.73152};
    
    double pointwise_momentum_1[24] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235};
    double pointwise_mean_1[24] = {3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.87467, 3.81896, 3.76826, 3.71725, 3.66345, 3.60687, 3.54826};
    double pointwise_1sig_1[24] = {3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.80467, 3.74896, 3.70126, 3.65398, 3.60227, 3.54655, 3.48891};
    
    double pointwise_momentum_2[32] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.205, 0.215, 0.225, 0.235, 0.245, 0.255, 0.265, 0.275, 0.285, 0.295, 0.305, 0.315};
    double pointwise_mean_2[32] = {3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.42323, 3.37141, 3.32037, 3.26727, 3.21402, 3.16228, 3.11243};
    double pointwise_1sig_2[32] = {3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.35323, 3.30141, 3.25037, 3.19792, 3.14654, 3.09698, 3.04921};
    
    // Parameters from Bichsel curve fits (H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195) 
    double meanBichselParameters[numRapidityWindows][6] = {
      {-0.878993,-2.317102,1.722582,1.131647,2.220628,0.938000},
      {-0.189078,-1.908437,1.721055,1.147400,2.525231,0.938000},
      {-0.201507,-1.955316,1.625150,1.205107,2.420558,0.938000},
      {-5.105297,-4.642927,1.968073,0.958584,2.077747,0.938000},
      {-5.090398,-4.677399,1.952695,0.863164,2.250091,0.938000},
      {-5.078298,-4.744772,1.912748,0.657149,2.486766,0.938000},
      {-5.044129,-4.816336,1.857193,0.199179,2.895074,0.938000},
      {-2.055400,-3.892668,1.267918,-0.338285,3.492077,0.938000},
      {0.115918,-3.330273,0.816915,-0.749739,3.901644,0.938000},
      {-9.467907,-7.285930,1.984900,0.705668,2.633194,0.938000},
      {-9.670624,-8.318109,1.667156,1.196107,2.541994,0.938000},
      {-7.983004,-7.772080,1.494031,1.363233,2.594165,0.938000},
      {-4.618411,-6.698504,1.134318,1.952621,2.673747,0.938000},
      {-2.918459,-7.065833,0.844667,1.953744,2.994180,0.938000},
      {-3.030857,-6.939436,0.866608,1.491180,3.116373,0.938000},
      {-4.487478,-5.637798,1.355464,1.133786,2.853095,0.938000},
      {-4.525409,-5.446432,1.435297,1.149635,2.880920,0.938000}
    };
    double sig1BichselParameters[numRapidityWindows][6] = {
      {0.010680,-3.345721,0.791844,1.141033,2.392233,0.653260},
      {2.632076,-2.518065,0.374227,1.368559,2.766669,0.555622},
      {-5.155654,-4.647320,1.945173,0.891506,2.147227,0.893156},
      {-5.152697,-4.648901,1.942672,0.870726,2.212334,0.904123},
      {-5.206472,-4.723595,1.919716,0.764264,2.366571,0.906424},
      {-5.133057,-4.544795,2.021496,0.523465,2.623447,0.946695},
      {-3.614053,-3.271727,2.472167,0.115895,3.188194,1.144807},
      {-3.893105,-3.493578,2.356486,-1.101787,4.693655,1.172724},
      {-3.835657,-3.098863,2.942318,-1.382106,5.278179,1.360753},
      {-12.259466,-8.502958,2.128565,0.662530,2.589974,0.908177},
      {-7.948264,-7.194245,1.619977,1.135863,2.607592,0.915820},
      {-8.205761,-7.608636,1.554044,1.315144,2.595075,0.930365},
      {-10.661223,-8.815795,1.693172,1.352099,2.532983,1.016327},
      {-7.006314,-6.497088,1.672735,1.666056,2.822535,1.189091},
      {-6.573408,-6.101615,1.720192,1.776161,2.888317,1.336846},
      {-6.233619,-5.769010,1.758411,1.894112,2.652388,1.370290},
      {-5.086869,-4.686561,1.965677,1.962093,2.724107,1.392304}
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
  }; // End struct NewNSigmaProton4p5GeV



}// End namespace NSigmaCorrectionUtils

#endif
