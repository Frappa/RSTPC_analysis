//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 29 15:57:19 2018 by ROOT version 6.10/08
// from TTree mppc/mppc
// found on file: Hangover_run4.root
//////////////////////////////////////////////////////////

#ifndef MPPCTREEWRAPPER_HH
#define MPPCTREEWRAPPER_HH

#include "RSTPC_Globals.hh"

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

// Header file for the classes stored in the TTree if any.

class MppcTreeWrapper
{
public:
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain
	
	string fDataDir;
	
	// Fixed size dimensions of array or collections stored in the TTree if any.
	
	// Declaration of leaf types
	UChar_t         mac5;
	UShort_t        chg[32];
	UInt_t          flags;
	UInt_t          lostcpu;
	UInt_t          lostfpga;
	UInt_t          ts0;
	UInt_t          ts1;
	UInt_t          ts0_ref;
	UInt_t          ts1_ref;
	UInt_t          sec;
	UInt_t          msec;
	
	// List of branches
	TBranch        *b_mac5;   //!
	TBranch        *b_chg;   //!
	TBranch        *b_flags;   //!
	TBranch        *b_lostcpu;   //!
	TBranch        *b_lostfpga;   //!
	TBranch        *b_ts0;   //!
	TBranch        *b_ts1;   //!
	TBranch        *b_ts0_ref;   //!
	TBranch        *b_ts1_ref;   //!
	TBranch        *b_sec;   //!
	TBranch        *b_msec;   //!
	
	MppcTreeWrapper();
	MppcTreeWrapper(TTree *tree);
	MppcTreeWrapper(string filename);
	MppcTreeWrapper(TFile* f);
	virtual ~MppcTreeWrapper();
	
	virtual Bool_t   Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual Bool_t   Open(string filename);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
	Double_t         GetEventTime(); //Returns the event in seconds
	Bool_t           IsInit(){return fClassInit;}

private:
	virtual Bool_t   Init(TTree *tree);
	Bool_t fClassInit;
};


#endif // #ifndef MPPCTREEWRAPPER_HH
