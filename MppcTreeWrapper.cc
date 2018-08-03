#ifndef MPPCTREEWRAPPER_CC
#define MPPCTREEWRAPPER_CC

#include "MppcTreeWrapper.hh"

#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"


using namespace std;

MppcTreeWrapper::MppcTreeWrapper() : fChain(0)
{
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fDataDir = RSTPC_Options::GetInstance()->GetDataDir();
	}
	
	fClassInit = false;
}

MppcTreeWrapper::MppcTreeWrapper(TTree *tree) : fChain(0)
{
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fDataDir = RSTPC_Options::GetInstance()->GetDataDir();
	}
	
	if(!tree)
	{
		fClassInit = false;
		return;
	}
	
	fClassInit = Init(tree);
}

MppcTreeWrapper::MppcTreeWrapper(string filename) : fChain(0) 
{
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fDataDir = RSTPC_Options::GetInstance()->GetDataDir();
	}
	
	fClassInit = Open(filename);
}

MppcTreeWrapper::MppcTreeWrapper(TFile* f) : fChain(0) 
{
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fDataDir = RSTPC_Options::GetInstance()->GetDataDir();
	}
	
	if(!f)
	{
		fClassInit = false;
		return;
	}
	
	fChain = (TTree*)f->Get("mppc");
	if(!fChain)
	{
		fClassInit = false;
		return;
	}
	
	fClassInit = Init(fChain);
}

MppcTreeWrapper::~MppcTreeWrapper()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Bool_t MppcTreeWrapper::Open(string filename)
{
	TFile *f = TFile::Open( filename.c_str(), "read" );
	if(!f)
	{
		fClassInit = false;
		return false;
	}
	
	fChain = (TTree*)f->Get("mppc");
	if(!fChain)
	{
		fClassInit = false;
		return false;
	}
	
	fClassInit = Init(fChain);
	return fClassInit;
}

Bool_t MppcTreeWrapper::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).
	
	// Set branch addresses and branch pointers
	if (!tree) return false;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);
	
	fChain->SetBranchAddress("mac5", &mac5, &b_mac5);
	fChain->SetBranchAddress("chg", chg, &b_chg);
	fChain->SetBranchAddress("flags", &flags, &b_flags);
	fChain->SetBranchAddress("lostcpu", &lostcpu, &b_lostcpu);
	fChain->SetBranchAddress("lostfpga", &lostfpga, &b_lostfpga);
	fChain->SetBranchAddress("ts0", &ts0, &b_ts0);
	fChain->SetBranchAddress("ts1", &ts1, &b_ts1);
	fChain->SetBranchAddress("ts0_ref", &ts0_ref, &b_ts0_ref);
	fChain->SetBranchAddress("ts1_ref", &ts1_ref, &b_ts1_ref);
	fChain->SetBranchAddress("sec", &sec, &b_sec);
	fChain->SetBranchAddress("msec", &msec, &b_msec);
	Notify();
	
	return true;
}

Int_t MppcTreeWrapper::GetEntry(Long64_t entry)
{
// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}

Long64_t MppcTreeWrapper::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent)
	{
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

Bool_t MppcTreeWrapper::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.
	
	return kTRUE;
}

void MppcTreeWrapper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}

Bool_t MppcTreeWrapper::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return true;
}

Double_t MppcTreeWrapper::GetEventTime()
{
	return ((Double_t)sec) + (((Double_t)msec)/1000.);
}


void MppcTreeWrapper::Loop()
{
//   In a ROOT session, you can do:
//      root> .L mppctree.C
//      root> mppctree t
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
	if (!fChain) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	/*
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
	}
	*/
}



#endif /* MPPCTREEWRAPPER_CC */

