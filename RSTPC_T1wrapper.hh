//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 23 10:01:42 2018 by ROOT version 6.10/08
// from TTree T1/Tear 1 data of the merged tree for the Resistive Shell TPC data
// found on file: RSTPC_Run000002032_Merged.root
//////////////////////////////////////////////////////////

#ifndef RSTPC_T1WRAPPER_HH
#define RSTPC_T1WRAPPER_HH

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

// Header file for the classes stored in the TTree if any.
#include "TH2.h"

class RSTPC_T1wrapper
{//This class is used to perform the analysis on the tear 1 data and eventually to produce tear 2 data tree
public :
	TFile *fInfile;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   string fMergedDataDir;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           TpcEvent;
   Double_t        TpcTime;
   Double_t        RmsColWires[32];
   Double_t        RmsIndWires[32];
   TH2D            *ColHist;
   TH2D            *IndHist;
   Int_t           FebEvent;
   Double_t        FebTime;
   UShort_t        FebTopAmp[3];
   Double_t        FebTopTotAmp;
   UShort_t        FebBotAmp[3];
   Double_t        FebBotTotAmp;
   Double_t        FebTotAmp;

   // List of branches
   TBranch        *b_TpcEvent;   //!
   TBranch        *b_TpcTime;
   TBranch        *b_RmsColWires;   //!
   TBranch        *b_RmsIndWires;   //!
   TBranch        *b_ColHist;   //!
   TBranch        *b_IndHist;   //!
   TBranch        *b_FebEvent;   //!
   TBranch        *b_FebTime;   //!
   TBranch        *b_FebTopAmp;   //!
   TBranch        *b_FebTopTotAmp;   //!
   TBranch        *b_FebBotAmp;   //!
   TBranch        *b_FebBotTotAmp;   //!
   TBranch        *b_FebTotAmp;   //!


	RSTPC_T1wrapper();
	RSTPC_T1wrapper(TTree *tree);
	RSTPC_T1wrapper(string filename);
	RSTPC_T1wrapper(TFile* f);
	virtual ~RSTPC_T1wrapper();
	
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual Bool_t   Init(TTree *tree);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
	virtual Bool_t   Open(string filename);
	
	Bool_t           IsInit(){return fClassInit;};
	void SetFileOwner(Bool_t flag=true){fInFileOwner = flag;};
	Bool_t IsFileOwner(){return fInFileOwner;};

private:
	Bool_t fClassInit;
	Bool_t fInFileOwner;
};

#endif /* RSTPC_T1WRAPPER_HH */
