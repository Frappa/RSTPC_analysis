#ifndef RSTPC_RUNPROCESSOR_HH
#define RSTPC_RUNPROCESSOR_HH

#include "RSTPC_Globals.hh"

#include "RSTPC_Analyser.hh"
#include "MppcTreeWrapper.hh"
#include "RSTPC_T1wrapper.hh"
//#include "RSTPC_Hits.hh"

#include "HistoManipulators.hh"
#include "DigitalFilters.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"

#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>


typedef enum
{
	kUndef=0,
	kCol=1,
	kInd=2
} WireType;


class RSTPC_Pulse;
class RSTPC_Hit;
class EventData;

class RSTPC_RunProcessor
{
public:
	vector<Double_t> *fTimeDiffVec;
	vector<Double_t> *fTpcTimes;
	vector<Double_t> *fFebTimes;
	vector<Int_t> *fFebEvs;
	
	
	RSTPC_RunProcessor();
	RSTPC_RunProcessor(Int_t RunNumber);
	virtual ~RSTPC_RunProcessor();
	
	
	void LoadColMap(string filename)
	{
		fTpcMan->LoadCollMap(filename);
	};
	void LoadIndMap(string filename)
	{
		fTpcMan->LoadIndcMap(filename);
	};
	
	Bool_t InitT1proc(Int_t RunNumber=-1);
	virtual void DescribeT1();
	virtual void T1process();
	
	//Int_t CheckTimesConsistency();
	
	Bool_t InitT2proc(Int_t RunNumber){return false;}; //Not implemented yet
	Bool_t InitT2proc();
	virtual void DescribeT2();
	virtual void T2Process();
	Bool_t IsProcT2init(){return fProcT2;};
	
	map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >* CombinePulses(vector<RSTPC_Pulse*>* ColPulses,  vector<RSTPC_Pulse*>* IndPulses, Bool_t debug=false);
	
	
	
	Bool_t IsRunOpen() const {return fRunOpenFlag;};
	
	//Getters to debug interactively the class
	const string GetDataDir() const {return fDataDir;};
	const string GetOutDir() const {return fOutDir;};
	
	TFile* GetOutFile() const {return fOutFile;};
	RSTPC_T1wrapper* GetT1wrapper() const {return fT1wr;};
	
	EventData* GetEventData() const {return gEventData;};
	
	
	
	//Static functions. They must be defined in .cc file otherwise ROOT does not compile (actually link) the classes including this header!!!
	static void SetDebug(Bool_t flag=true, Int_t nEvs = 100);
	
	static void SetSigmaThr(Double_t _val);
	static Double_t GetSigmaThr();
	
	static void SetPeakingTime(Double_t _val);
	static Double_t GetPeakingTime();
	
	static void SetSamplingFreq(Double_t _val);
	static Double_t GetSamplingFreq();
	
	static void SetPitchSize(Double_t _val);
	static Double_t GetPitchSize();
	
	static void SetDriftLenght(Double_t _val);
	static Double_t GetDriftLenght();
	
	static void SetDriftVel(Double_t _val);
	static Double_t GetDriftVel();
	
	
protected:
	
	//This methods should only be used by the T2Process function and not by the user.
	void FindPulses(TH2D* h, WireType type, Bool_t debug=false);
	
	Int_t HitsFinder(map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >* pulseMap, Bool_t debug=false);
	
	
	RSTPC_Analyser *fTpcMan;
	MppcTreeWrapper *fFebMan;
	
	string fDataDir, fOutDir;
	
	TFile *fOutFile;
	TTree *fOutT1, *fOutT2;
	
	//This is used to read the T1 tree, is allocated by the AddressT1() method and in general used for the T2 processing
	RSTPC_T1wrapper *fT1wr;
	
	Bool_t fRunOpenFlag;
	
	Int_t fRunNumber;
	
	Bool_t fProcT1, fProcT2; //It is true when the T1 tree can be processed
	
	//This actually behaves as a global variable
	EventData *gEventData; //It should only be allocated and deallocated by the T1/2process methods
	
	//This is actually a global variable
	vector<vector<Double_t> > gPulseWfs;
	
	//Here the waveforms of each channel are stored
	vector<vector<Double_t> > *gColWfsVect, *gIndWfsVect;
	
	//Vectors with the collection of all the pulses
	vector<RSTPC_Pulse*> *gColPulses, *gIndPulses;
	
private:
	
	static Bool_t fgDebug;
	static Int_t fgNevPrint;//Frequency of debug print outs
	
	static Double_t fgSigmaThr;
	
	static Double_t fgPeakingTime; //In micro-secs
	
	static Double_t fgSamplingFreq; //In 1./micro-secs
	
	static Double_t fgPitchSize; //In mm
	static Double_t fgDriftLenght; //In mm
	static Double_t fgDriftVel; //In mm/usec
};


class EventData
{//this is actually a structure to set branch addresses of the output tree
public:
	Int_t TpcEv, FebEv;
	Double_t TpcTime, FebTime, RmsColWires[32], RmsIndWires[32];
	UShort_t FebTopAmp[3], FebBotAmp[3];
	Double_t FebTopTotAmp, FebBotTotAmp, FebTotAmp;
	
	//Objects handled (owned) from ouside
	TH2D *hCol;
	TH2D *hInd;
	
	//This stuff is for the T2 processing and this class is the only owner
	TObjArray *ColPulses, *IndPulses;
	TObjArray *Hits;
	
	Bool_t GoodEvent;
	
	EventData();
	virtual ~EventData();
	
	void Reset();
};

/*
#ifndef RSTPC_RUNPROCESSOR_CC
//Initialise the static variables here
#define RSTPC_RUNPROC_STATICS
Double_t RSTPC_RunProcessor::fgSigmaThr = 3.0;

Double_t RSTPC_RunProcessor::fgPeakingTime = 0.0;
Double_t RSTPC_RunProcessor::fgSamplingFreq = 1./50e-3;
Double_t RSTPC_RunProcessor::fgPitchSize = 1.875;
Double_t RSTPC_RunProcessor::fgDriftLenght = 150;
Double_t RSTPC_RunProcessor::fgDriftVel = 0.0;
#endif
*/

#endif /* RSTPC_RUNPROCESSOR_HH */