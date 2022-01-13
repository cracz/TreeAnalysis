// Tiny DST class for calculating event planes
#ifndef IEVENTPLANE_H
#define IEVENTPLANE_H

#include "TObject.h"

class IEventPlane : public TObject {

	private:
		double mPhi;
		double mWeight;
		
		double mEta;
		int mRingNumber;
		double mpT;
		int mPID;
		int mTPCCharge;
		double mToF;
		
		void Construct();


	public:
		IEventPlane();
		IEventPlane(double, double);
		
		~IEventPlane(){} //This is needed apparently

		void SetWeight(double in_weight){mWeight = in_weight;}
		
		void SetEta(double in_eta){mEta = in_eta;}
		void SetRingNumber(int in_ring){mRingNumber = in_ring;}
		void SetTransverseMomentum(double in_pT){mpT = in_pT;}
		void SetParticleID(int in_PID){mPID = in_PID;}
		void SetTPCCharge(int in_TPCCharge){mTPCCharge = in_TPCCharge;}
		void SetTimeOfFlight(double in_ToF){mWeight = in_ToF;}		
		
		double GetPhi(){return mPhi;}
		double GetWeight(){return mWeight;}
		
		double GetEta(){return mEta;}
		int GetRingNumber(){return mRingNumber;}
		double GetTransverseMomentum(){return mpT;}
		int GetParticleID(){return mPID;}
		int GetTPCCharge(){return mTPCCharge;}
		double GetTimeOfFlight(){return mToF;}			
		
		double QxTerm(int);
		double QyTerm(int);

		ClassDef(IEventPlane,1)  // my EventPlane
};
    
#endif
