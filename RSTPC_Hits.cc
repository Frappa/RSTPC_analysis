#ifndef RSTPC_HITS_CC
#define RSTPC_HITS_CC


#include "RSTPC_Hits.hh"

#include "Rtypes.h"


ClassImp(RSTPC_Pulse)
ClassImp(RSTPC_Hit)


//#ifndef RSTPC_HITS_STATICS
UInt_t RSTPC_Pulse::fgPulses = 0;
//#endif



UInt_t RSTPC_Pulse::GetNpulses()
{
	return fgPulses;
}

void RSTPC_Pulse::ResetCounter()
{
	fgPulses = 0;
}

RSTPC_Pulse::RSTPC_Pulse():fWireType(kUndef)
{
	fPulseID = 0;
	return;
}

RSTPC_Pulse::RSTPC_Pulse(WireType _type)
{
	fWireType = _type;
	fPulseID = fgPulses+1;
	fgPulses++;
}

RSTPC_Pulse::~RSTPC_Pulse()
{
	TObject::Clear();
}

RSTPC_Pulse& RSTPC_Pulse::operator=(const RSTPC_Pulse &orig)
{
	fWireType = orig.fWireType;
	
	fWireNum = orig.fWireNum; //Corresponds also to the channel number after the channels map is applied
	
	fPulseID = orig.fPulseID; //This is specific to the hit number on this wire (0 is the earliest in time)
	fMax = orig.fMax;
	fMin = orig.fMin;
	
	fMaxPos = orig.fMaxPos;
	fMinPos = orig.fMinPos;
	fLedge = orig.fLedge;
	fRedge = orig.fRedge;
	
	return *this;
}


const Bool_t RSTPC_Pulse::operator<(const RSTPC_Pulse& right) const
{
	/*
	if( (fgPeakingTime>0.) && (abs(left.fCentre-right.fCentre)<=RSTPC_RunProcessor::GetPeakingTime()) )
	{
		return left.fWireNum<right.fWireNum;
	}
	return left.fCentre<right.fMassCentre;
	*/
	return (*this).fPulseID < right.fPulseID;
}

const Bool_t RSTPC_Pulse::operator==(const RSTPC_Pulse& right) const
{
	return (*this).fPulseID == right.fPulseID;
}



//#ifndef RSTPC_HITS_STATICS
UInt_t RSTPC_Hit::fgNhits = 0;
//#endif


UInt_t RSTPC_Hit::GetNhits()
{
	return fgNhits;
}

void RSTPC_Hit::ResetCounter()
{
	fgNhits=0;
}

RSTPC_Hit::RSTPC_Hit(const RSTPC_Pulse* ColPulse, const RSTPC_Pulse* IndPulse)
{
	if( !(ColPulse && IndPulse) ) return;
	
	fHitID = fgNhits+1;
	fgNhits++;
	
	fColPulseID = ColPulse->fPulseID;
	fIndPulseID = IndPulse->fPulseID;
	
	fColWireNum = ColPulse->fWireNum;
	fIndWireNum = IndPulse->fWireNum;
	
	fX = fIndWireNum*RSTPC_RunProcessor::GetPitchSize();
	fY = fColWireNum*RSTPC_RunProcessor::GetPitchSize();
}


RSTPC_Hit::RSTPC_Hit()
{
	fHitID = 0;
}


RSTPC_Hit::~RSTPC_Hit()
{
	TObject::Clear();
}


RSTPC_Hit& RSTPC_Hit::operator=(const RSTPC_Hit &orig)
{
	TObject::operator=(orig);
	
	fHitID = orig.fHitID;
	
	fColPulseID = orig.fColPulseID;
	fIndPulseID = orig.fIndPulseID;
	
	fColWireNum = orig.fColWireNum;
	fIndWireNum = orig.fIndWireNum;
	
	fX = orig.fX;
	fY = orig.fY;
	fZ = orig.fZ;
	
	fMeanHeight = orig.fMeanHeight;
	fCentreTime = orig.fCentreTime;
	fMeanTime = orig.fMeanTime;
	
	fLedge = orig.fLedge;
	fRedge = orig.fRedge;
	
	return *this;
}


const Bool_t RSTPC_Hit::operator<(const RSTPC_Hit& right) const
{
	return (*this).fHitID < right.fHitID;
}


const Bool_t RSTPC_Hit::operator==(const RSTPC_Hit& right) const
{
	return (*this).fHitID == right.fHitID;
}


void RSTPC_Hit::SetCentreTime(Double_t centre)
{
	fCentreTime = centre; //This is units of samples
	fZ = centre * RSTPC_RunProcessor::GetDriftVel() / RSTPC_RunProcessor::GetSamplingFreq();
}


#endif /* RSTPC_HITS_CC */