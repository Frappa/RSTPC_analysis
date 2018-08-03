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
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
	RSTPC_T1wrapper * fT1wr;
// Fixed size dimensions of array or collections stored in the TTree if any.
	
   // Declaration of leaf types
   Bool_t          GoodEvent;
   TObjArray       *ColPulses;
   TObjArray       *IndPulses;
   TObjArray       *Hits;

   // List of branches
   TBranch        *b_GoodEvent;   //!
   TBranch        *b_ColPulses;   //!
   TBranch        *b_IndPulses;   //!
   TBranch        *b_Hits;   //!

   RSTPC_T2wrapper(TTree *tree=0);
   virtual ~RSTPC_T2wrapper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   void DrawPulsesLenghts();
   void DrawPulsesWfs(Int_t iEv);
   void TestPulsesOverlap(Int_t event);
   void TestPulsesOverlap(Int_t event, UInt_t colwire, UInt_t indwire);
   void PlotColPulses();
};

#endif /* RSTPC_T2WRAPPER_HH */

#ifdef RSTPC_T2WRAPPER_CC
RSTPC_T2wrapper::RSTPC_T2wrapper(TTree *tree) : fChain(0), fT1wr(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	
	if(!TClassTable::GetDict("RSTPC_Hits")) {
		gSystem->Load("RSTPC_Hits");
	}
	
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root");
      if (!f || !f->IsOpen()) {
         f = TFile::Open("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root");
      }
      f->GetObject("T2",tree);

   }
   Init(tree);
   
   fT1wr = new RSTPC_T1wrapper((TTree*)tree->GetCurrentFile()->Get("T1"));
   fT1wr->SetFileOwner(false);
}

RSTPC_T2wrapper::~RSTPC_T2wrapper()
{
	if(fT1wr) delete fT1wr;
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t RSTPC_T2wrapper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RSTPC_T2wrapper::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}


void RSTPC_T2wrapper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ColPulses = 0;
   IndPulses = 0;
   Hits = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("GoodEvent", &GoodEvent, &b_GoodEvent);
   fChain->SetBranchAddress("ColPulses", &ColPulses, &b_ColPulses);
   fChain->SetBranchAddress("IndPulses", &IndPulses, &b_IndPulses);
   fChain->SetBranchAddress("Hits", &Hits, &b_Hits);
   Notify();
}


Bool_t RSTPC_T2wrapper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RSTPC_T2wrapper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RSTPC_T2wrapper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RSTPC_T2WRAPPER_CC
