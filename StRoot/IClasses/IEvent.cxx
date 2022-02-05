#include "IEvent.h"
#include "stdio.h"
#include <iostream>

#include "IEventPlane.h"
#include "TMath.h"


//EPD geometry header file
#include "../StEpdUtil/StEpdGeom.h"

ClassImp(IEvent)

//_________________
IEvent::IEvent(){
	ClearEvent();
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

	/*
	mRunNumber = -99;
	mChronologicalRunId = -99;
	mEventID = -99;
	mRefMult = -99;
	mCentrality = -99;
	mPrimaryVertex.SetXYZ(-99.0,-99.0,-99.0);
	mBfield = -99.;
	mVpdVz = -99.;
	nEPDhits = -99;
	mTriggers.clear();
	EPDid.fill(0);
	EPDnMip.fill(0);
	*/
	
	mqCenter[0] = 0.0;
	mqCenter[1] = 0.0;
	
	mEPParticles.clear();
}


void IEvent::CalcQVector(int harmonic, std::vector<float> &nqx, std::vector<float> &nqy){
	std::vector<float> Qx;
	std::vector<float> Qy;
	
	
	Qx.push_back(-1*mqCenter[0]);
	Qy.push_back(-1*mqCenter[1]);
	
	for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
		float sign = 1;
	
		//For odd harmonics, flip the sign in the negative rapidity region
		if (harmonic % 2 == 1 && mEPParticles[parts].GetEta() < 0){sign = -1;}
		
		float xterm = sign*mEPParticles[parts].QxTerm(harmonic);
		float yterm = sign*mEPParticles[parts].QyTerm(harmonic);
		
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

Float_t IEvent::GetQx(int harmonic){
	std::vector<float> Qx;
	std::vector<float> Qy;
	
	CalcQVector(harmonic, Qx, Qy);
	
	return Qx[0];
}

Float_t IEvent::GetQy(int harmonic){
	std::vector<float> Qx;
	std::vector<float> Qy;
	
	CalcQVector(harmonic, Qx, Qy);
	
	return Qy[0];
}

//Returns a vector for the psis for the event. Index 0 is for the whole event and all others
//are with autocorrelations removed.
std::vector<float>  IEvent::CalcPsi(int harmonic){
	std::vector<float> Qx;
	std::vector<float> Qy;
	
	std::vector<float> Psi;
	
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

std::vector<float>  IEvent::GetPsi(int harmonic){	
	return CalcPsi(harmonic);
}

std::vector<float>  IEvent::GetPsi(int harmonic, char option, float param1, float param2 = 0.0){
	
	IEvent subEvent = GetSubEvent(option, param1, param2);

	return subEvent.CalcPsi(harmonic);
}

//Any sort of subevent can be constructed in this manner
//But these are the provided types
IEvent IEvent::GetSubEvent(char option, float param1, float param2 = 0.0){
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
		case 'p' : //p_T
			for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
				if (mEPParticles[parts].GetMomentum().Pt() >= param1 && mEPParticles[parts].GetMomentum().Pt() <= param2){
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
				if (mEPParticles[parts].GetCharge() >= param1 && mEPParticles[parts].GetCharge() <= param2){
					mEPParticlesCopy.push_back(mEPParticles[parts]);
				}
			}				
			break;
		case 't' : //ToF	
			for (unsigned int parts = 0; parts < mEPParticles.size(); parts++){
				if (mEPParticles[parts].GetToFBeta() >= param1 && mEPParticles[parts].GetToFBeta() <= param2){
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


void IEvent::SetQCenter(float nqx, float nqy){
	mqCenter[0] = nqx;
	mqCenter[1] = nqy;
}

//Convert EPD nMIP information into an IEventPlane vector
std::vector<IEventPlane> IEvent::EPDVector(TVector3 primaryVertex, float COMrapidity = 0.0){
	std::vector<IEventPlane> mEPParticlesCopy;
	
	Double_t mThresh = 0.3; // EPD EP by hand
	Double_t mMax = 3.0; // EPD EP by hand
	StEpdGeom *mEpdGeom = new StEpdGeom();
	
	int num_sides = 1; //1 for east only, 2 for east and west
	
	for (int ew = 0; ew < num_sides; ew++){
		for (int pp = 1; pp <=12; pp++){
			for (int tt = 1; tt <=31; tt++){
				
				//Get hit
				if (EPDnMip[ew][pp - 1][tt - 1]){

					//Construct tileID using ew pp tt
					int tileID = (pp * 100 + tt);
					if (ew == 0){tileID *= -1;}
					
					Float_t nMip = EPDnMip[ew][pp - 1][tt - 1];


					if (nMip<mThresh) continue;
					float TileWeight = (nMip<mMax)?nMip:mMax;
					
					TVector3 StraightLine = mEpdGeom->TileCenter(tileID) - primaryVertex;
					float phi = StraightLine.Phi();
					float eta = StraightLine.Eta() - COMrapidity;
					
					IEventPlane eventPlane(phi, TileWeight);
					eventPlane.SetEta(eta);
					eventPlane.SetTileID(tileID);		
					eventPlane.SetnMIP(nMip);				
				
					mEPParticlesCopy.push_back(eventPlane);		
				}
			}
		}
	}
	
	delete mEpdGeom;
	
	return mEPParticlesCopy;
}

void IEvent::setEPDnMip(int in_tileID, float nMip){
	int ew;
	if (in_tileID > 0){ew = 1;}
	else {ew = 0;}
	
	int tt = TMath::Abs(in_tileID) % 100;
	int pp = (TMath::Abs(in_tileID) - tt) / 100;
	
	setEPDnMip(ew, pp, tt, nMip);
}

//Call this this once per event at maximum!
void IEvent::AddEPDtoTracks(TVector3 primaryVertex, float COMrapidity = 0.0){
	std::vector<IEventPlane> epdtracks = EPDVector(primaryVertex, COMrapidity);
	  
	for (unsigned int hit = 0; hit < epdtracks.size(); hit++){
		AddEPParticle(epdtracks[hit]);
	}
	
}