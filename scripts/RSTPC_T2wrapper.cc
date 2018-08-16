#define RSTPC_T2WRAPPER_CC
#include "RSTPC_T2wrapper.hh"
#include "TH2.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <sstream>
#include <vector>


RSTPC_T2wrapper::RSTPC_T2wrapper() : fChain(0), fT1wr(0), fInFile(0), fInFileOwner(false), fClassInit(false), fT1wrOwner(false)
{
	return;
}


RSTPC_T2wrapper::RSTPC_T2wrapper(string filename, Bool_t bMakeT1w) : fChain(0), fT1wr(0), fInFile(0), fInFileOwner(false), fClassInit(false), fT1wrOwner(false)
{
	TFile* f = (TFile*)gROOT->FindObject( filename.c_str() );
	if(f && f->IsOpen())
	{//The owner is someone else
		fInFile = f;
	}
	else
	{
		fInFile = TFile::Open( filename.c_str(), "read" );
		if(!fInFile || !fInFile->IsOpen()) return;
		fInFileOwner = true;
	}
	
	TTree* tree = (TTree*)fInFile->Get("T2");
	if(!tree)
	{
		if(fInFileOwner)
		{
			delete fInFile;
			fInFile = NULL;
			fInFileOwner = false;
		}
		return;
	}
	
	Init(tree);
	
	if(fClassInit)
	{
		fChain = tree;
		
		if(bMakeT1w)
		{
			fT1wr = new RSTPC_T1wrapper(fChain);
			fT1wr->SetFileOwner(false);
			if(!fT1wr->IsInit())
			{
				delete fT1wr;
				fT1wr = NULL;
			}
			fT1wrOwner = true;
		}
	}
}


RSTPC_T2wrapper::RSTPC_T2wrapper(TFile* f, Bool_t bMakeT1w) : fChain(0), fT1wr(0), fInFile(0), fInFileOwner(false), fClassInit(false), fT1wrOwner(false)
{
	if( (!f) || (!f->IsOpen()) ) return;
	
	fInFile = f;
	
	
	TTree* tree = (TTree*)fInFile->Get("T2");
	if(!tree)
	{
		if(fInFileOwner)
		{
			delete fInFile;
			fInFile = NULL;
			fInFileOwner = false;
		}
		return;
	}
	
	Init(tree);
	
	if(fClassInit)
	{
		fChain = tree;
		
		if(bMakeT1w)
		{
			fT1wr = new RSTPC_T1wrapper(fChain);
			fT1wr->SetFileOwner(false);
			if(!fT1wr->IsInit())
			{
				delete fT1wr;
				fT1wr = NULL;
			}
			fT1wrOwner = true;
		}
	}
}


RSTPC_T2wrapper::RSTPC_T2wrapper(RSTPC_T1wrapper* t1w) : fChain(0), fT1wr(0), fInFile(0), fInFileOwner(false), fClassInit(false), fT1wrOwner(false)
{
	if(!TClassTable::GetDict("RSTPC_Hits")) {
		gSystem->Load("RSTPC_Hits");
	}
	
	if(!t1w) return;
	
	if(!t1w->IsInit()) return;
	
	if(!t1w->fInfile) return;
	
	if(!t1w->fInfile->IsOpen()) return;
	
	fInFile = t1w->fInfile;
	
	fChain = (TTree*)fInFile->Get("T2");
	
	if(!fChain) return;
	
	Init(fChain);
	
	if(!fClassInit) return;
	
	fT1wr = t1w;
	fT1wrOwner = false;
	fInFileOwner = false;
}


RSTPC_T2wrapper::RSTPC_T2wrapper(TTree *tree, Bool_t bMakeT1w) : fChain(0), fT1wr(0), fInFile(0), fClassInit(false), fT1wrOwner(false)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	
	if (!tree) return;
	
	if(!TClassTable::GetDict("RSTPC_Hits")) {
		gSystem->Load("RSTPC_Hits");
	}
	
	
	Init(tree);
	
	fInFile = fChain->GetCurrentFile();
	fInFileOwner = false;
	
	if(fClassInit && bMakeT1w)
	{
		fT1wr = new RSTPC_T1wrapper((TTree*)tree->GetCurrentFile()->Get("T1"));
		fT1wr->SetFileOwner(false);
		if(!fT1wr->IsInit())
		{
			delete fT1wr;
			fT1wr = NULL;
		}
		fT1wrOwner = true;
	}
}


RSTPC_T2wrapper::~RSTPC_T2wrapper()
{
	if(fInFileOwner)
	{
		if (!fChain) return;
		delete fChain->GetCurrentFile();
	}
	
	if(fT1wr && fT1wrOwner)
	{
		delete fT1wr;
	}
}


void RSTPC_T2wrapper::SetFileOwner(Bool_t flag)
{
	if(fT1wr)
	{
		if(fT1wrOwner)
		{
			fT1wr->SetFileOwner(false);
		}
		else
		{
			if(!fT1wr->IsFileOwner())
			{
				fInFileOwner = flag;
			}
		}
	}
	fInFileOwner = flag;
}


Int_t RSTPC_T2wrapper::GetEntry(Long64_t entry)
{
// Read contents of entry.
	if (!fChain) return 0;
	Int_t readBytes = fChain->GetEntry(entry);
	
	if(fT1wr && fT1wr->IsInit() && fT1wrOwner)
	{
		readBytes += fT1wr->GetEntry(entry);
	}
	
	return readBytes;
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
	
	fClassInit = false;
	
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
	
	fClassInit = true;
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


void RSTPC_T2wrapper::Loop()
{
//   In a ROOT session, you can do:
//      root> .L T2.C
//      root> T2 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
