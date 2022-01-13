#include "TMath.h"
#include "IEventPlane.h"

ClassImp(IEventPlane)

//______________
IEventPlane::IEventPlane(){
	mPhi = 0;
	SetWeight(0);
	Construct();
}

//______________
IEventPlane::IEventPlane(double in_phi, double in_weight){
	mPhi = in_phi;
	SetWeight(in_weight);
	Construct();
}

void IEventPlane::Construct(){
	mEta = 0.0;
	mRingNumber = 0;
	mpT = 0.0;
	mPID = 0;
	mTPCCharge = 0;
	mToF = 0.0;
}

//______________
double IEventPlane::QxTerm(int harmonic){
	double Cosine = cos(mPhi*(double)harmonic);
	return mWeight * Cosine;
}

//______________
double IEventPlane::QyTerm(int harmonic){
	double Sine = sin(mPhi*(double)harmonic);
	return mWeight * Sine;
}
