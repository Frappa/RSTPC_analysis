#include "RSTPC_Globals.hh"
#include "RSTPC_RunProcessor.hh"
#include "RSTPC_Analyser.hh"
#include "MppcTreeWrapper.hh"

#include "TChain.h"
#include "TFile.h"
#include "TEventList.h"

void RSTPC_FebData()
{
	TChain *T1 = new TChain("T1");
	
	T1->Add("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002031_Merged.root");
	T1->Add("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root");
	T1->Add("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002033_Merged.root");
	
	
	TCanvas* c1 = new TCanvas("canvFeb","Feb data",800,800);
	T1->Draw("FebTopTotAmp:FebBotTotAmp>>hFEB(100,0,12500,100,0,12500)","","colz");
}


void DumpEventLists()
{
	TChain *T1 = new TChain("T1");
	
	T1->Add("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root");
	
	T1->Draw(">>FebSatList1","(FebTopAmp>4070)||(FebBotAmp>4070)");
	T1->Draw(">>FebSatList2","(FebTopTotAmp>12250)||(FebBotTotAmp>12250)");
	T1->Draw(">>FebSmallSigList","(FebTopTotAmp<1000)&&(FebBotTotAmp<1000)");
	
	TEventList *FebSatList1 = (TEventList*)gROOT->FindObject("FebSatList1");
	TEventList *FebSatList2 = (TEventList*)gROOT->FindObject("FebSatList2");
	TEventList *FebSmallSigList = (TEventList*)gROOT->FindObject("FebSmallSigList");
	
	TFile* outfile = TFile::Open("FebListRun2032.root","recreate");
	outfile->WriteTObject(FebSatList1);
	outfile->WriteTObject(FebSatList2);
	outfile->WriteTObject(FebSmallSigList);
}
