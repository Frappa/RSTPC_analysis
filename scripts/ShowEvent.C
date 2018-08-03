#include "RSTPC_Globals.hh"
#include "RSTPC_Analyser.hh"
#include "HistoManipulators.hh"
#include "DigitalFilters.hh"

#include "TChain.h"
#include "TFile.h"
#include "TEventList.h"

#include <sstream>
//R__LOAD_LIBRARY(RSTPCAnalyser.so)
//R__LOAD_LIBRARY(HistoManipulators.so)
//R__LOAD_LIBRARY(DigitalFilters.so)
Int_t iEvent = 0;


RSTPC_Analyser *an = NULL;
TEventList* evList = NULL;

void ShowEvent(Int_t iEv){
	
	
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell");
	
	if(!an)
	{
		an = new RSTPC_Analyser;
		an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
		an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
		an->OpenRun(2032);
	}
	
	//an->OpenRun(2002); //Calibration Run
	//an->OpenRun(2004); //Calibration Run
	an->SetBaselineROI(1000, 4000);
	an->LoadEvent(iEv);
	an->ApplyWiresMaps(an->hC0,an->hI0);
	
	an->hC0->GetXaxis()->SetTitle("Drift time [samples]");
	an->hC0->GetYaxis()->SetTitle("Wire number");
	
	an->hI0->GetXaxis()->SetTitle("Drift time [samples]");
	an->hI0->GetYaxis()->SetTitle("Wire number");
	
	cout << "\nShowing event " << iEv <<endl;
	
	TCanvas* c0 = (TCanvas*)gROOT->FindObject("c0");
	if(!c0)
	{
		c0 = new TCanvas("c0","Before any correction", 1020,640);
		c0->Divide(1,2);
	}
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	
	TH2D* hC0_0 = (TH2D*)an->hC0->Clone("hC0_0");
	
	TH2D* hI0_0 = (TH2D*)an->hI0->Clone("hI0_0");
	
	
	c0->cd(1);
	hC0_0->Draw("colz");
	
	c0->cd(2);
	hI0_0->Draw("colz");
	
	
	
	TCanvas* c1 = (TCanvas*)gROOT->FindObject("c1");
	if(!c1)
	{
		c1 = new TCanvas("c1","After CM rejection", 1020,640);
		c1->Divide(1,2);
	}
	
	an->CMrej(an->hC0, an->fColRMS, 5);
	an->CMrej(an->hI0, an->fIndRMS, 5);
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH2D* hC0_1 = (TH2D*)an->hC0->Clone("hC0_1");
	
	TH2D* hI0_1 = (TH2D*)an->hI0->Clone("hI0_1");
	
	c1->cd(1);
	hC0_1->Draw("colz");
	
	c1->cd(2);
	hI0_1->Draw("colz");
	
	
	
	TCanvas* c2 = (TCanvas*)gROOT->FindObject("c2");
	if(!c2)
	{
		c2 = new TCanvas("c2","Overlapped waveforms = before any correction", 1020,640);
		c2->Divide(1,2);
	}
	
	
	
	
	
	Int_t nChs = hC0_0->GetNbinsY();
	
	TLegend *leg = new TLegend(0.7,0.9,0.9,0.65);
	
	stringstream hname;
	for(Int_t iCh=0; iCh<nChs; iCh++)
	{
		hname.str(""); hname << "colWf0_ch" << iCh;
		TH1D *colWf_0 = hC0_0->ProjectionX( hname.str().c_str(), iCh+1, iCh+1 ); colWf_0->SetLineColor(kBlue);
		hname.str(""); hname << "colWf1_ch" << iCh;
		TH1D *colWf_1 = hC0_1->ProjectionX( hname.str().c_str(), iCh+1, iCh+1 ); colWf_1->SetLineColor(kRed);
		
		c2->cd(1);
		if(iCh==0)
		{
			colWf_0->Draw();
			leg->AddEntry(colWf_0, "Before CM rejection", "L");
			leg->AddEntry(colWf_1, "After CM rejection", "L");
		}
		else
		{
			colWf_0->Draw("same");
		}
		colWf_1->Draw("same");
		
		
		hname.str(""); hname << "indWf0_ch" << iCh;
		TH1D *indWf_0 = hI0_0->ProjectionX( hname.str().c_str(), iCh+1, iCh+1 ); ; indWf_0->SetLineColor(kBlue);
		hname.str(""); hname << "indWf1_ch" << iCh;
		TH1D *indWf_1 = hI0_1->ProjectionX( hname.str().c_str(), iCh+1, iCh+1 ); ; indWf_1->SetLineColor(kRed);
		
		c2->cd(2);
		if(iCh==0)
		{
			indWf_0->Draw();
		}
		else
		{
			indWf_0->Draw("same");
		}
		indWf_1->Draw("same");
	}
	
	c2->cd(1);
	leg->Draw();
	
	return;
}


void NextEvent()
{
	if(!evList)
	{
		ShowEvent(iEvent);
	}
	else
	{
		ShowEvent( evList->GetEntry(iEvent) );
	}
	iEvent++;
}




void SumAllEvents(Int_t RunNumber)
{//This is mostly useful for calibration runs
	
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell");
	
	if(!an)
	{
		an = new RSTPC_Analyser;
		an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
		an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
	}
	
	an->Set_CMnoiseRej(true, 5);
	an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
	an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
	
	an->OpenRun(RunNumber);
	
	//an->OpenRun(2001); //Calibration Run
	//an->OpenRun(2002); //Calibration Run
	//an->OpenRun(2004); //Calibration Run
	
	//an->OpenRun(2032);
	
	an->SetBaselineROI(1000, 4000);
	
	
	
	
	Int_t nEvs = an->fTrigTree->GetEntries();
	
	stringstream hname;
	TH1D *hColSum = NULL, *hIndSum = NULL;
	
	Int_t nChs;
	
	for(Int_t iEv=0; iEv<nEvs; iEv++)
	{
		an->LoadEvent(iEv);
		an->ApplyWiresMaps(an->hC0,an->hI0);
		
		an->BaselineCorr(an->hC0);
		an->BaselineCorr(an->hI0);
		
		nChs = an->hC0->GetNbinsY();
		if(!hColSum)
		{
			hColSum = an->hC0->ProjectionX("hColSum", 1, nChs);
		}
		else
		{
			hColSum->Add( an->hC0->ProjectionX("hColSum_tmp", 1, nChs) );
		}
		
		nChs = an->hI0->GetNbinsY();
		if(!hIndSum)
		{
			hIndSum = an->hI0->ProjectionX("hIndSum", 1, nChs);
		}
		else
		{
			hIndSum->Add( an->hI0->ProjectionX("hIndSum_tmp", 1, nChs) );
		}
	}
	
	
	TCanvas* c3 = (TCanvas*)gROOT->FindObject("c3");
	if(!c3)
	{
		c3 = new TCanvas("c3","Summed waveforms = before any correction", 1020,640);
		c3->Divide(1,2);
	}
	
	c3->cd(1);
	hColSum->Scale(1./nChs/nEvs);
	hColSum->Draw("hist");
	
	c3->cd(2);
	hIndSum->Scale(1./nChs/nEvs);
	hIndSum->Draw("hist");
	
	if(gROOT->FindObject("hColSum_tmp")) delete gROOT->FindObject("hColSum_tmp");
	if(gROOT->FindObject("hIndSum_tmp")) delete gROOT->FindObject("hIndSum_tmp");
	
	return;
}

