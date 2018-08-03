#ifndef RSTPC_T1WRAPPER_CC
#define RSTPC_T1WRAPPER_CC


#include "RSTPC_T1wrapper.hh"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <sstream>


using namespace std;


RSTPC_T1wrapper::RSTPC_T1wrapper() : fChain(0), fInfile(NULL), fInFileOwner(true)
{
	fClassInit = false;
	
	fInfile = NULL;
	
	/*
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fMergedDataDir = RSTPC_Options::GetInstance()->GetMergedDataDir();
	}
	*/
	
	fMergedDataDir = "/home/francescop/data/ResistiveShell/merged/";
	
}


RSTPC_T1wrapper::RSTPC_T1wrapper(TTree *tree) : fChain(0), fInfile(NULL), fInFileOwner(true)
{
	//fMergedDataDir = RSTPC_Options::GetInstance()->GetMergedDataDir();
	fMergedDataDir = "/home/francescop/data/ResistiveShell/merged/";
	
	if(!tree)
	{
		fClassInit = false;
		return;
	}
	
	fInfile = tree->GetCurrentFile();
	
	fClassInit = Init(tree) && (fInfile!=NULL);
}


RSTPC_T1wrapper::RSTPC_T1wrapper(string filename) : fChain(0) , fInfile(NULL), fInFileOwner(true)
{
	/*
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fMergedDataDir = RSTPC_Options::GetInstance()->GetDataDir();
	}
	*/
	
	fMergedDataDir = "/home/francescop/data/ResistiveShell/merged/";
	
	fClassInit = Open(filename);
}

RSTPC_T1wrapper::RSTPC_T1wrapper(TFile* f) : fChain(0) , fInfile(NULL), fInFileOwner(true)
{
	fClassInit = false;
	
	if( (!f) || (!f->IsOpen()) ) return;
	
	fChain = (TTree*)f->Get("T1");
	if(fChain) fClassInit = Init(fChain);
}


RSTPC_T1wrapper::~RSTPC_T1wrapper()
{
	if (!fChain) return;
	if(fInFileOwner)
	{
		if(fChain->GetCurrentFile()) delete fChain->GetCurrentFile();
	}
}


Int_t RSTPC_T1wrapper::GetEntry(Long64_t entry)
{
// Read contents of entry.
	if (!fChain) return 0;
	
	stringstream ColHistName, IndHistName;
	
	ColHistName.str(""); ColHistName << "Col_" << entry;
	IndHistName.str(""); IndHistName << "Ind_" << entry;
	
	ColHist = (TH2D*)fInfile->Get(ColHistName.str().c_str());
	IndHist = (TH2D*)fInfile->Get(IndHistName.str().c_str());
	
	if( !(ColHist && IndHist) ) return -100;
	
	return fChain->GetEntry(entry);
}


Long64_t RSTPC_T1wrapper::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
	if(!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if(centry < 0) return centry;
	if(fChain->GetTreeNumber() != fCurrent)
	{
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}


Bool_t RSTPC_T1wrapper::Open(string filename)
{
	TFile *f = TFile::Open( filename.c_str(), "read" );
	if(!f)
	{
		fClassInit = false;
		return false;
	}
	
	fChain = (TTree*)f->Get("T1");
	if(!fChain)
	{
		fClassInit = false;
		return false;
	}
	
	fInfile = f;
	
	fClassInit = Init(fChain);
	return fClassInit;
}


Bool_t RSTPC_T1wrapper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
	
	if (!tree) return false;
	
   // Set object pointer
   ColHist = NULL;
   IndHist = NULL;
   // Set branch addresses and branch pointers
   
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
	
	//This in the case the tree was just written and now needs to re-addressed to be read
	fChain->ResetBranchAddresses();
	
   fChain->SetBranchAddress("TpcEvent", &TpcEvent, &b_TpcEvent);
   fChain->SetBranchAddress("TpcTime", &TpcTime, &b_TpcTime);
   fChain->SetBranchAddress("RmsColWires", RmsColWires, &b_RmsColWires);
   fChain->SetBranchAddress("RmsIndWires", RmsIndWires, &b_RmsIndWires);
   //fChain->SetBranchAddress("ColHist", &ColHist, &b_ColHist);
   //fChain->SetBranchAddress("IndHist", &IndHist, &b_IndHist);
   fChain->SetBranchAddress("FebEvent", &FebEvent, &b_FebEvent);
   fChain->SetBranchAddress("FebTime", &FebTime, &b_FebTime);
   fChain->SetBranchAddress("FebTopAmp", FebTopAmp, &b_FebTopAmp);
   fChain->SetBranchAddress("FebTopTotAmp", &FebTopTotAmp, &b_FebTopTotAmp);
   fChain->SetBranchAddress("FebBotAmp", FebBotAmp, &b_FebBotAmp);
   fChain->SetBranchAddress("FebBotTotAmp", &FebBotTotAmp, &b_FebBotTotAmp);
   fChain->SetBranchAddress("FebTotAmp", &FebTotAmp, &b_FebTotAmp);
   Notify();
   
   fInfile = fChain->GetCurrentFile();
   
   return true;
}

Bool_t RSTPC_T1wrapper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RSTPC_T1wrapper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RSTPC_T1wrapper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void RSTPC_T1wrapper::Loop()
{
//   In a ROOT session, you can do:
//      root> .L RSTPC_T1wrapper.C
//      root> RSTPC_T1wrapper t
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


#endif /* RSTPC_T1WRAPPER_CC */