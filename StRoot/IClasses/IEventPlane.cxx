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
IEventPlane::IEventPlane(float in_phi, float in_weight){
	mPhi = in_phi;
	SetWeight(in_weight);
	Construct();
}

void IEventPlane::Construct(){
	mEta = 0.0;
	mMomentum.SetXYZ(0,0,0);
	mPID = 0;
	mCharge = 0;
	mToFBeta = 0.0;
	nHitsFit = 0;
	nHitsPoss = 0;
	mtileID = 0;
	nMIP = 0.0;
}

Float_t IEventPlane::GetRapidity(){ //Must have PID and momentum
	
	Float_t mass = 0;
	
	switch (abs(mPID)){
		case 2212 : //Proton
			mass = 0.938272081;
			break;
		case 211 : //Pi+-
			mass = 0.13957061;
			break;
		case 321 : //K +-
			mass = 0.493677;
			break;
		default:
			return 0;
	}


	Float_t energy = TMath::Sqrt(mMomentum.Mag2() + mass*mass);
	
	return 0.5*TMath::Log( (energy + mMomentum.Z()) / (energy - mMomentum.Z()) );
}

//______________
Float_t IEventPlane::QxTerm(int harmonic){
	Float_t Cosine = cos(mPhi*(float)harmonic);
	return mWeight * Cosine;
}

//______________
Float_t IEventPlane::QyTerm(int harmonic){
	Float_t Sine = sin(mPhi*(float)harmonic);
	return mWeight * Sine;
}

//0 = east, 1 = west
Int_t IEventPlane::GetEPDew(){
	if (mtileID > 0){return 1;}
	return 0;
}