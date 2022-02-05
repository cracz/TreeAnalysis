// author: Erik Loyd
// Based off of code by Mike Lisa

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

/*
  Int_t mRunNumber;
  Int_t mChronologicalRunId;
  Int_t mEventID;  // within run
  Float_t mBfield;
  TVector3 mPrimaryVertex;
  Int_t mRefMult;
  Int_t mCentrality;
  Float_t mVpdVz;
  std::vector<unsigned int> mTriggers;
  Short_t             EPDid[744];
  UShort_t            nEPDhits;

  UShort_t            tree_tracknumber; //What is this?
 
*/


  Float_t EPDnMip[2][12][31]; //0 = E, 1 = W
 
  Float_t mqCenter[2];  // This is the Q recentering offset for all events in this data set
  std::vector<IEventPlane> mEPParticles;
  
  void CalcQVector(int, std::vector<float>&, std::vector<float>&);
  std::vector<float> CalcPsi(int);
 
 public:

  IEvent();
  ~IEvent();

  void ClearEvent();
  void Init();

  
  void AddEPParticle(IEventPlane n){mEPParticles.push_back(n);}
  
  void SetEPParticle(std::vector<IEventPlane> n){mEPParticles = n;}
  
  std::vector<float> GetPsi(int);
  std::vector<float> GetPsi(int, char, float, float);
  Float_t GetQx(int);
  Float_t GetQy(int);

/*
  TVector3 PrimaryVertex() {return mPrimaryVertex;}
  void SetPrimaryVertex(TVector3 pv){mPrimaryVertex = pv;}
  Int_t GetRunNumber(){return mRunNumber;}
  Int_t GetChronologicalRunId(){return mChronologicalRunId;}
  Int_t GetEventID(){return mEventID;}
  Int_t GetRefMult(){return mRefMult;}
  Int_t GetCentralityID9(){return mCentrality}
  Float_t Bfield(){return mBfield;}
  Float_t GetVpdVz(){return mVpdVz;}
  
  */
  //std::vector<unsigned int> GetTriggers(){return mTriggers;}
  std::vector<IEventPlane> GetEPParticles(){return mEPParticles;}

/*
  void SetRunNumber(Int_t rn){mRunNumber=rn;}
  void SetChronologicalRunId(Int_t rn){mChronologicalRunId=rn;}
  void SetEventID(Int_t en){mEventID=en;}
  void SetRefMult(Int_t rm){mRefMult = rm;}
  void SetCentrality(Int_t cent){mCentrality=cent;}
  void SetBfield(Float_t f){mBfield=f;}
  void SetVpdVz(Float_t n){mVpdVz=n;}
  void SetTriggers( const std::vector<unsigned int> & n){mTriggers = n;}
*/ 
 
  void SetEPParticles( const std::vector<IEventPlane> & n){mEPParticles = n;}
  void SetQCenter(float, float);
  
  std::vector<IEventPlane> EPDVector(TVector3, float);
  
  void setEPDnMip(int ew, int pp, int tt, float nMip){EPDnMip[ew][pp - 1][tt - 1] = nMip;}
  void setEPDnMip(int tileID, float nMip);
  
  Float_t GetEPDnMip(int ew, int pp, int tt){return EPDnMip[ew][pp - 1][tt - 1];}

  IEvent GetSubEvent(char, float, float);
  
  void AddEPDtoTracks(TVector3, float);;

  ClassDef(IEvent,1)  // my event
};

#endif


