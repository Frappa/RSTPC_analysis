#ifndef RSTPC_HITS_HH
#define RSTPC_HITS_HH

#include "RSTPC_RunProcessor.hh"

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"



class RSTPC_Pulse: public TObject
{
private:
	static UInt_t fgPulses; //! transient value
	
public:
	WireType fWireType;
	
	UInt_t fWireNum; //Corresponds also to the channel number after the channels map is applied
	
	UInt_t fPulseID; //This is specific to the hit number on this wire (0 is the earliest in time)
	Double_t fMax, fMin;
	
	Int_t fMaxPos, fMinPos;
	Int_t fLedge, fRedge;
	
	
	RSTPC_Pulse(); //This needed by the rootsystem in order to save this class in a tree
	RSTPC_Pulse(WireType _type);
	~RSTPC_Pulse();
	
	
	RSTPC_Pulse &operator=(const RSTPC_Pulse &orig);
	const Bool_t operator==(const RSTPC_Pulse& right) const;
	const Bool_t operator<(const RSTPC_Pulse& right) const;
	
	static UInt_t GetNpulses();
	static void ResetCounter();
	
	ClassDef(RSTPC_Pulse,1)
};


class RSTPC_Hit: public TObject
{
private:
	//Counter for getting unique IDs. It is increased only by the constructor.
	static UInt_t fgNhits; //! transient value
	
public:
	
	UInt_t fHitID; //This is specific and unique to the hit number on this wire (0 means that the hit is not part of the collection of hits)
	Int_t fColWireNum, fIndWireNum;
	UInt_t fColPulseID, fIndPulseID;
	Double_t fX, fY;
	
	//The following variables are directly set from outside the class
	Double_t fZ;
	Double_t fMeanHeight; //This is related to the collection charge in the specific time slice (will be used as a weight)
	Double_t fCentreTime;
	Double_t fMeanTime; //This is weighted with the charge but should be very close to the "centre time"
	Double_t fLedge, fRedge; //Boundaries of the time slice (in samples units)
	
	
	//This constructor doesn't increase the counter 
	RSTPC_Hit();
	
	//This constructor increases the static counter. (For debugging purpouses)
	RSTPC_Hit(const RSTPC_Pulse* ColPulse, const RSTPC_Pulse* IndPulse);
	
	~RSTPC_Hit();
	
	RSTPC_Hit& operator=(const RSTPC_Hit &orig);
	
	const Bool_t operator<(const RSTPC_Hit& right) const;
	const Bool_t operator==(const RSTPC_Hit& right) const;
	
	
	static UInt_t GetNhits();
	static void ResetCounter();
	
	void SetCentreTime(Double_t centre); //In samples units. It should correspond to the center of a time slice.
	
	ClassDef(RSTPC_Hit,1)

};

/*
#ifndef RSTPC_HITS_CC
#define RSTPC_HITS_STATICS
UInt_t RSTPC_Pulse::fgPulses = 0;
UInt_t RSTPC_Hit::fgNhits = 0;
#endif
*/


#endif /* RSTPC_HIT_HH */