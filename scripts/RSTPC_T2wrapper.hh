//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul 29 19:59:36 2018 by ROOT version 6.10/08
// from TTree T2/Tear 2 data of the merged tree for the Resistive Shell TPC data
// found on file: RSTPC_Run000002032_Merged.root
//////////////////////////////////////////////////////////

#ifndef RSTPC_T2WRAPPER_HH
#define RSTPC_T2WRAPPER_HH

#include "RSTPC_T1wrapper.hh"
#include "RSTPC_Hits.hh"

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

// Header file for the classes stored in the TTree if any.
#include "TObjArray.h"
#include "TClassTable.h"

class RSTPC_T2wrapper{
public :
	TTree           *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t            fCurrent; //!current Tree number in a TChain
	RSTPC_T1wrapper *fT1wr;
	TFile           *fInFile;
// Fixed size dimensions of array or collections stored in the TTree if any.
	
	// Declaration of leaf types
	Bool_t           GoodEvent;
	TObjArray       *ColPulses;
	TObjArray       *IndPulses;
	TObjArray       *Hits;

	// List of branches
	TBranch        *b_GoodEvent;   //!
	TBranch        *b_ColPulses;   //!
	TBranch        *b_IndPulses;   //!
	TBranch        *b_Hits;   //!
	
	
	RSTPC_T2wrapper();
	RSTPC_T2wrapper(string filename, Bool_t bMakeT1w=false);
	RSTPC_T2wrapper(TFile* f, Bool_t bMakeT1w=false);
	RSTPC_T2wrapper(TTree *tree=0, Bool_t bMakeT1w=false);
	RSTPC_T2wrapper(RSTPC_T1wrapper* t1w);


	virtual ~RSTPC_T2wrapper();
	
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
	
	
	void SetFileOwner(Bool_t flag=true);
	
	Bool_t           IsT1wrOwner(){return fT1wrOwner;};
	Bool_t           IsInFileOwner(){return fInFileOwner;};
	Bool_t           IsInit(){return fClassInit;};

private:
	Bool_t fInFileOwner;
	Bool_t fT1wrOwner;
	Bool_t fClassInit;
};

#endif /* RSTPC_T2WRAPPER_HH */
