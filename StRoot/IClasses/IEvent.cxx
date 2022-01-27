#include "IEvent.h"
#include "stdio.h"
#include <iostream>

#include "ITrack.h"
#include "IParticle.h"
#include "IEventPlane.h"


ClassImp(IEvent)

//_________________
IEvent::IEvent(){
	mRunNumber = -99;
	mChronologicalRunId = -99;
	mEventID = -99;
	mRefMult = -99;
	mCentralityID9 = -99;
	mCentralityID16 = -99;
	mPrimaryVertex.SetXYZ(-99.0,-99.0,-99.0);
	mBfield = -99.;
	mVpdVz = -99.;
	mqCenter[0] = 0.0;
	mqCenter[1] = 0.0;
}

//__________________
IEvent::~IEvent(){
  ClearEvent();
  
}

//_________________
void IEvent::Init(){
}
//_________________
void IEvent::ClearEvent(){

  mRunNumber = -99;
  mChronologicalRunId = -99;
  mEventID = -99;
  mRefMult = -99;
  mCentralityID9 = -99;
  mCentralityID16 = -99;
  mPrimaryVertex.SetXYZ(-99.0,-99.0,-99.0);
  mBfield = -99.;
  mVpdVz = -99.;
  mTriggers.clear();
  mEPParticles.clear();
  mqCenter[0] = 0.0;
  mqCenter[1] = 0.0;
}

void IEvent::CalcQVector(int harmonic, std::vector<double> &nqx, std::vector<double> &nqy){
	std::vector<double> Qx;
	std::vector<double> Qy;
	
	
	Qx.push_back(-1*mqCenter[0]);
	Qy.push_back(-1*mqCenter[1]);
	
	for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
		double sign = 1;
	
		//For odd harmonics, flip the sign in the negative rapidity region
		if (harmonic % 2 == 1 && mEPParticles[parts].GetEta() < 0){sign = -1;}
		
		double xterm = sign*mEPParticles[parts].QxTerm(harmonic);
		double yterm = sign*mEPParticles[parts].QyTerm(harmonic);
		
		Qx[0] += xterm;
		Qy[0] += yterm;
		
		Qx.push_back(-1*xterm);
		Qy.push_back(-1*yterm);
	}
	
	for (unsigned int parts = 1; parts <= mEPParticles.size(); parts++){
		Qx[parts] += Qx[0];
		Qy[parts] += Qy[0];
	}
	
	nqx = Qx;
	nqy = Qy;
}

double IEvent::GetQx(int harmonic){
	std::vector<double> Qx;
	std::vector<double> Qy;
	
	CalcQVector(harmonic, Qx, Qy);
	
	return Qx[0];
}

double IEvent::GetQy(int harmonic){
	std::vector<double> Qx;
	std::vector<double> Qy;
	
	CalcQVector(harmonic, Qx, Qy);
	
	return Qy[0];
}

//Returns a vector for the psis for the event. Index 0 is for the whole event and all others
//are with autocorrelations removed.
std::vector<double>  IEvent::CalcPsi(int harmonic){
	std::vector<double> Qx;
	std::vector<double> Qy;
	
	std::vector<double> Psi;
	
	CalcQVector(harmonic, Qx, Qy);
	
	for (unsigned int parts = 0; parts <= mEPParticles.size(); parts++){
		if(Qx[parts] || Qy[parts] ){
			Psi.push_back(TMath::ATan2(Qy[parts],Qx[parts])/((Double_t)harmonic));
		}
		else{
			Psi.push_back(-99); //This shouldn't happen, but bad psi for this particle
			std::cout << "Warning! Error particle [" << parts;
			std::cout << " / " << mEPParticles.size() << "] in this event!" << std::endl;
		}
		
	}
	
	return Psi;
}

std::vector<double>  IEvent::GetPsi(int harmonic){	
	return CalcPsi(harmonic);
}

std::vector<double>  IEvent::GetPsi(int harmonic, char option, double param1, double param2 = 0.0){
	
	IEvent subEvent = GetSubEvent(option, param1, param2);

	return subEvent.CalcPsi(harmonic);
}

IEvent IEvent::GetSubEvent(char option, double param1, double param2 = 0.0){
	IEvent subEvent = *this; //Use itself as a base; this probably black magic
	
	std::vector<IEventPlane> mEPParticlesCopy;
	
	
	switch (option){
		case 'e' : //eta
			for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
				if (mEPParticles[parts].GetEta() >= param1 && mEPParticles[parts].GetEta() <= param2){
					mEPParticlesCopy.push_back(mEPParticles[parts]);
				}
			}
			break;
		case 'r' : //ring number
			for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
				if (mEPParticles[parts].GetRingNumber() >= (int)param1 && mEPParticles[parts].GetRingNumber() <= (int)param2){
					mEPParticlesCopy.push_back(mEPParticles[parts]);
				}
			}				
			break;
		case 'p' : //p_T
			for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
				if (mEPParticles[parts].GetTransverseMomentum() >= param1 && mEPParticles[parts].GetTransverseMomentum() <= param2){
					mEPParticlesCopy.push_back(mEPParticles[parts]);
				}
			}				
			break;
		case 'i' : //ID
			for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
				if (mEPParticles[parts].GetParticleID() == (int)param1){
					mEPParticlesCopy.push_back(mEPParticles[parts]);
				}
			}				
			break;
		case 'c' : //charge
			for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
				if (mEPParticles[parts].GetTPCCharge() >= param1 && mEPParticles[parts].GetTPCCharge() <= param2){
					mEPParticlesCopy.push_back(mEPParticles[parts]);
				}
			}				
			break;
		case 't' : //ToF	
			for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
				if (mEPParticles[parts].GetTimeOfFlight() >= param1 && mEPParticles[parts].GetTimeOfFlight() <= param2){
					mEPParticlesCopy.push_back(mEPParticles[parts]);
				}
			}		
			break;
		default: 
			//std::cerr << "Error! Not a valid option! Ignoring parameters!" << std::endl;
			mEPParticlesCopy = mEPParticles;
	}
	
	subEvent.SetEPParticles(mEPParticlesCopy);
	
	return subEvent;
}


void IEvent::SetQCenter(double nqx, double nqy){
	mqCenter[0] = nqx;
	mqCenter[1] = nqy;
}
