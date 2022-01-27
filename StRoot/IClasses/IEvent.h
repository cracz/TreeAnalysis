// author: Mike Lisa 14 feb 2013 // edited for just lambdas Dec 5 by Isaac Upsal

#ifndef IEVENT_H
#define IEVENT_H

#include "TVector3.h"
#include "TClonesArray.h"
#include "TObject.h"
#include <vector>
#include "IEventPlane.h"


/*
  Store simple event summary information
*/

class IEvent : public TObject {

 protected:

  Int_t mRunNumber;
  Int_t mChronologicalRunId;
  Int_t mEventID;  // within run
  Int_t mRefMult;
  Int_t mCentralityID9;		//centrality ID number from StRefMultCorr RefMult9
  Int_t mCentralityID16;	//centrality ID number from StRefMultCorr RefMult16
  TVector3 mPrimaryVertex;
  Float_t mBfield;
  Float_t mVpdVz;
  double mqCenter[2];
  std::vector<unsigned int> mTriggers;
  
  std::vector<IEventPlane> mEPParticles;
  
  void CalcQVector(int, std::vector<double>&, std::vector<double>&);
  std::vector<double> CalcPsi(int);

 public:

  IEvent();
  ~IEvent();

  void ClearEvent();
  void Init();

  
  void AddEPParticle(IEventPlane n){mEPParticles.push_back(n);}
  void AddEPParticle(std::vector<IEventPlane> n){mEPParticles = n;}
  
  std::vector<double> GetPsi(int);
  std::vector<double> GetPsi(int, char, double, double);
  double GetQx(int);
  double GetQy(int);

  TVector3 PrimaryVertex() {return mPrimaryVertex;}
  void SetPrimaryVertex(TVector3 pv){mPrimaryVertex = pv;}

  Int_t GetRunNumber(){return mRunNumber;}
  Int_t GetChronologicalRunId(){return mChronologicalRunId;}
  Int_t GetEventID(){return mEventID;}
  Int_t GetRefMult(){return mRefMult;}
  Int_t GetCentralityID9(){return mCentralityID9;}
  Int_t GetCentralityID16(){return mCentralityID16;}
  Float_t Bfield(){return mBfield;}
  Float_t GetVpdVz(){return mVpdVz;}
  std::vector<unsigned int> GetTriggers(){return mTriggers;}
  std::vector<IEventPlane> GetEPParticles(){return mEPParticles;}


  void SetRunNumber(Int_t rn){mRunNumber=rn;}
  void SetChronologicalRunId(Int_t rn){mChronologicalRunId=rn;}
  void SetEventID(Int_t en){mEventID=en;}
  void SetRefMult(Int_t rm){mRefMult = rm;}
  void SetCentralityID9(Int_t cent){mCentralityID9=cent;}
  void SetCentralityID16(Int_t cent){mCentralityID16=cent;}
  void SetBfield(Float_t f){mBfield=f;}
  void SetVpdVz(Float_t n){mVpdVz=n;}
  void SetTriggers( const std::vector<unsigned int> & n){mTriggers = n;}
  void SetEPParticles( const std::vector<IEventPlane> & n){mEPParticles = n;}
  void SetQCenter(double, double);

  IEvent GetSubEvent(char, double, double);

  ClassDef(IEvent,1)  // my event
};

#endif


