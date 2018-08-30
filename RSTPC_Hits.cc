#ifndef RSTPC_HITS_CC
#define RSTPC_HITS_CC


#include "RSTPC_Hits.hh"

#include "Rtypes.h"


ClassImp(RSTPC_Pulse)

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
	
	fFWHM = -1;
	fFWTM = -1;
	fMeanTime = -1;
	fSigma = -1;
}

RSTPC_Pulse::RSTPC_Pulse(const RSTPC_Pulse &orig)
{
	*this = orig;
}

RSTPC_Pulse::RSTPC_Pulse(WireType _type)
{
	fWireType = _type;
	fPulseID = fgPulses+1;
	fgPulses++;
	
	fFWHM = -1;
	fFWTM = -1;
	fMeanTime = -1;
	fSigma = -1;
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
	
	fFWHM = orig.fFWHM;
	fFWTM = orig.fFWTM;
	fMeanTime = orig.fMeanTime;
	fSigma = orig.fSigma;
	
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


Double_t RSTPC_Pulse::SetFWHM(vector<Double_t>* flWf)
{
	if(!flWf) return -1;
	
	if(fWireType!=kCol)
	{//It makes sense only for unipolar collection pulses
		fFWHM = -1;
		return -1;
	}
	
	double xlow, xup;
	
	Int_t nSamps = flWf->size();
	
	//Find the left bin where the amplitude is half of the maximum
	Int_t jSamp = fMaxPos;
	do{
		jSamp--;
		if( (jSamp<0) || (jSamp==fLedge) ) break;
	}while( (flWf->at(jSamp) > 0.5*fMax) && (jSamp>=fLedge) );
	
	if( jSamp<=0 ){
		xlow = jSamp;
	}else if( flWf->at(jSamp) == 0.5*fMax ){
		xlow = jSamp;
	}else{
		xlow = jSamp + 0.5; //Take the middle point with the previous bin
	}
	
	//Find the right bin where the amplitude is half of the maximum
	jSamp = fMaxPos;
	do{
		jSamp++;
		if( (jSamp>=nSamps) || (jSamp==fRedge) ) break;
	}while( (flWf->at(jSamp) > 0.5*fMax) && (jSamp<=fRedge) );
	
	if( jSamp>=nSamps ){
		xup = jSamp;
	}else if( flWf->at(jSamp) == 0.5*fMax ){
		xup =jSamp;
	}else{
		xup = jSamp - 0.5; //Take the middle point with the previous bin
	}
	
	fFWHM = (xup-xlow);
	return fFWHM;
	
}


Double_t RSTPC_Pulse::SetFWTM(vector<Double_t>* flWf)
{
	if(!flWf) return -1;
	
	if(fWireType!=kCol)
	{//It makes sense only for unipolar collection pulses
		fFWTM = -1;
		return -1;
	}
	
	double xlow, xup;
	
	Int_t nSamps = flWf->size();
	
	//Find the left bin where the amplitude is half of the maximum
	Int_t jSamp = fMaxPos;
	do{
		jSamp--;
		if( (jSamp<0) || (jSamp==fLedge) ) break;
	}while( (flWf->at(jSamp) > 0.1*fMax) && (jSamp>=fLedge) );
	
	if( jSamp<=0 ){
		xlow = jSamp;
	}else if( flWf->at(jSamp) == 0.1*fMax ){
		xlow = jSamp;
	}else{
		xlow = jSamp + 0.1; //Take the middle point with the previous bin
	}
	
	//Find the right bin where the amplitude is half of the maximum
	jSamp = fMaxPos;
	do{
		jSamp++;
		if( (jSamp>=nSamps) || (jSamp==fRedge) ) break;
	}while( (flWf->at(jSamp) > 0.1*fMax) && (jSamp<=fRedge) );
	
	if( jSamp>=nSamps ){
		xup = jSamp;
	}else if( flWf->at(jSamp) == 0.1*fMax ){
		xup =jSamp;
	}else{
		xup = jSamp - 0.1; //Take the middle point with the previous bin
	}
	
	fFWTM = (xup-xlow);
	return fFWTM;
}


Double_t RSTPC_Pulse::SetMeanTime(vector<Double_t>* flWf)
{
	if(!flWf) return -1;
	
	Double_t val, sum2=0;
	
	fMeanTime = 0;
	
	for(Int_t iSamp=fLedge; iSamp<=fRedge; iSamp++)
	{
		val = pow(flWf->at(iSamp), 2);
		sum2 += val;
		fMeanTime += iSamp*val;
	}
	
	fMeanTime = fMeanTime/sum2;
	
	if(fMeanTime<=0) fMeanTime = ((Double_t)(fLedge+fRedge))/2.;
	
	return fMeanTime;
}


Double_t RSTPC_Pulse::SetSigma(vector<Double_t>* flWf)
{
	if(!flWf) return -1;
	
	if(fMeanTime<=0) SetMeanTime(flWf);
	
	Double_t val, sum2=0;
	
	fSigma = 0;
	
	for(Int_t iSamp=fLedge; iSamp<=fRedge; iSamp++)
	{
		val = pow(flWf->at(iSamp), 2);
		sum2 += val;
		fSigma += val*pow((iSamp-fMeanTime), 2);
	}
	
	fSigma = fSigma/sum2;
	
	//Sigma of a flat distribution
	if(fSigma<=0) fSigma = pow((Double_t)(fRedge-fLedge), 2)/12.;
	
	return fSigma;
}




ClassImp(RSTPC_Hit)

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


RSTPC_Hit::RSTPC_Hit(const RSTPC_Hit &orig)
{
	*this = orig;
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