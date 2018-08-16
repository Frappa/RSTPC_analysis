#include "HistoManipulators.hh"
#include "DigitalFilters.hh"

#include "RSTPC_Analyser.hh"
#include "RSTPC_T1wrapper.hh"
#include "RSTPC_T2wrapper.hh"
#include "RSTPC_Hits.hh"

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

// Header file for the classes stored in the TTree if any.
#include "TObjArray.h"
#include "TClassTable.h"


using namespace std;



RSTPC_T2wrapper *t2w = NULL;


void DrawPulsesLenghts()
{
	if(!t2w)
	{
		t2w = new RSTPC_T2wrapper("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root", true);
	}
	
	if( !t2w->IsInit() ) return;
	
	
	
	if (t2w->fChain == 0) return;
	
	Long64_t nEvs = t2w->fChain->GetEntriesFast();
	
	
	TH1D *hColPulsesLength = (TH1D*)gROOT->FindObject("hColPulsesLength");
	if(!hColPulsesLength)
	{
		hColPulsesLength = new TH1D("hColPulsesLength","Collection pulses;Pulse lenght (samples);Counts", 1000, -0.5, 1000-0.5);
	}
	else
	{
		hColPulsesLength->Reset();
	}
	
	TH1D *hIndPulsesLength = (TH1D*)gROOT->FindObject("hIndPulsesLength");
	if(!hIndPulsesLength)
	{
		hIndPulsesLength = new TH1D("hIndPulsesLength","Induction pulses;Pulse lenght (samples);Counts", 1000, -0.5, 1000-0.5);
	}
	else
	{
		hIndPulsesLength->Reset();
	}
	
	
	TH1D *hColPulsesLedges = (TH1D*)gROOT->FindObject("hColPulsesLedges");
	if(!hColPulsesLedges)
	{
		hColPulsesLedges = new TH1D("hColPulsesLedges","Collection pulses;Pulse left edge (samples);Counts", 5000, -0.5, 5000-0.5);
	}
	else
	{
		hColPulsesLedges->Reset();
	}
	
	TH1D *hIndPulsesLedges = (TH1D*)gROOT->FindObject("hIndPulsesLedges");
	if(!hIndPulsesLedges)
	{
		hIndPulsesLedges = new TH1D("hIndPulsesLedges","Induction pulses;Pulse left edge (samples);Counts", 5000, -0.5, 5000-0.5);
	}
	else
	{
		hIndPulsesLedges->Reset();
	}
	
	
	
	RSTPC_Pulse *ColPulse, *IndPulse;
	
	for(Int_t iEv=0; iEv<nEvs; iEv++)
	{
		t2w->GetEntry(iEv);
		
		TIter ColPulsesIt(t2w->ColPulses);
		TIter IndPulsesIt(t2w->IndPulses);
		
		while( (ColPulse = (RSTPC_Pulse*)ColPulsesIt.Next()) )
		{
			Double_t len = (Double_t)( ColPulse->fRedge-ColPulse->fLedge );
			hColPulsesLength->Fill(len);
			hColPulsesLedges->Fill((Double_t)ColPulse->fLedge);
		}
		
		while( (IndPulse = (RSTPC_Pulse*)IndPulsesIt.Next()) )
		{
			Double_t len = (Double_t)( IndPulse->fRedge-IndPulse->fLedge );
			hIndPulsesLength->Fill(len);
			hIndPulsesLedges->Fill((Double_t)IndPulse->fLedge);
		}
	}
	
	
	TCanvas *c1 = (TCanvas*)gROOT->FindObject("canvPulsesLen");
	if(c1) delete c1;
	c1 = new TCanvas("canvPulsesLen", "Pulses lenght", 900, 500);
	c1->Divide(2,1);
	
	c1->cd(1);
	hColPulsesLength->Draw();
	
	c1->cd(2);
	hIndPulsesLength->Draw();
	
	
	TCanvas *c2 = (TCanvas*)gROOT->FindObject("canvPulsesLedges");
	if(c2) delete c2;
	c2 = new TCanvas("canvPulsesLedges", "Pulses Left Edges", 900, 500);
	c2->Divide(2,1);
	
	c2->cd(1);
	hColPulsesLedges->Draw();
	
	c2->cd(2);
	hIndPulsesLedges->Draw();
	
}


void DrawPulsesWfs(Int_t iEv)
{
	if(!t2w)
	{
		t2w = new RSTPC_T2wrapper("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root", true);
	}
	
	if( !t2w->IsInit() ) return;
	
	if (!t2w->fChain) return;
	
	if( !((t2w->fT1wr) && (t2w->fT1wr->IsInit())) ) return;
	
	
	t2w->GetEntry(iEv);
	
	
	Int_t nChs = t2w->fT1wr->ColHist->GetNbinsY();
	Int_t nBins = t2w->fT1wr->ColHist->GetNbinsX();
	Double_t xlow = -0.5;
	Double_t xup = (nBins+1)-0.5;
	
	vector<TH1D*> colWfs(nChs), indWfs(nChs);
	vector<TH1D*> colWfsOrig(nChs), indWfsOrig(nChs);
	
	stringstream hname;
	
	//Initialise the histograms
	for(Int_t iCh=0; iCh<nChs; iCh++)
	{
		hname.str(""); hname << t2w->fT1wr->ColHist->GetName() << "_ch" << iCh;
		colWfsOrig.at(iCh) = t2w->fT1wr->ColHist->ProjectionX(hname.str().c_str(), iCh+1, iCh+1);
		
		hname.str(""); hname << "hColPulsesWfs_" << iCh;
		TH1D* hColPulsesWfs = (TH1D*)gROOT->FindObject(hname.str().c_str());
		if(hColPulsesWfs)
		{
			colWfs.at(iCh) = hColPulsesWfs;
			hColPulsesWfs->Reset();
		}
		else
		{
			colWfs.at(iCh) = (TH1D*)colWfsOrig.at(iCh)->Clone(hname.str().c_str());
			colWfs.at(iCh)->Reset();
			colWfs.at(iCh)->SetLineColor(kBlue);
		}
		
		
		hname.str(""); hname << t2w->fT1wr->IndHist->GetName() << "_ch" << iCh;
		indWfsOrig.at(iCh) = t2w->fT1wr->IndHist->ProjectionX(hname.str().c_str(), iCh+1, iCh+1);
		
		hname.str(""); hname << "hIndPulsesWfs_" << iCh;
		TH1D* hIndPulsesWfs = (TH1D*)gROOT->FindObject(hname.str().c_str());
		if(hIndPulsesWfs)
		{
			indWfs.at(iCh) = hIndPulsesWfs;
			hIndPulsesWfs->Reset();
		}
		else
		{
			indWfs.at(iCh) = (TH1D*)indWfsOrig.at(iCh)->Clone(hname.str().c_str());
			indWfs.at(iCh)->Reset();
			indWfs.at(iCh)->SetLineColor(kRed);
		}
	}
	
	
	
	TIter ColPulsesIt(t2w->ColPulses);
	TIter IndPulsesIt(t2w->IndPulses);
	RSTPC_Pulse *ColPulse, *IndPulse;
	
	//Now fill the histograms
	Int_t iPulse = 0;
	Double_t ymax, ymin;
	while( (ColPulse = (RSTPC_Pulse*)ColPulsesIt.Next()) )
	{
		UInt_t ch = ColPulse->fWireNum;
		for(Int_t iSamp=ColPulse->fLedge; iSamp<=ColPulse->fRedge; iSamp++)
		{
			colWfs.at(ch)->SetBinContent(iSamp+1, colWfsOrig.at(ch)->GetBinContent(iSamp+1));
		}
		
		if(iPulse == 0)
		{
			ymax = indWfs.at(ch)->GetMaximum();
			ymin = indWfs.at(ch)->GetMinimum();
		}
		else
		{
			if(ymax < indWfs.at(ch)->GetMaximum()) ymax = indWfs.at(ch)->GetMaximum();
			if(ymin > indWfs.at(ch)->GetMinimum()) ymin = indWfs.at(ch)->GetMinimum();
		}
		
		iPulse++;
	}
	
	
	while( (IndPulse = (RSTPC_Pulse*)IndPulsesIt.Next()) )
	{
		UInt_t ch = IndPulse->fWireNum;
		
		for(Int_t iSamp=IndPulse->fLedge; iSamp<=IndPulse->fRedge; iSamp++)
		{
			indWfs.at(ch)->SetBinContent(iSamp+1, indWfsOrig.at(ch)->GetBinContent(iSamp+1));
		}
		
		if(iPulse == 0)
		{
			ymax = indWfs.at(ch)->GetMaximum();
			ymin = indWfs.at(ch)->GetMinimum();
		}
		else
		{
			if(ymax < indWfs.at(ch)->GetMaximum()) ymax = indWfs.at(ch)->GetMaximum();
			if(ymin > indWfs.at(ch)->GetMinimum()) ymin = indWfs.at(ch)->GetMinimum();
		}
		
		iPulse++;
	}
	
	Double_t deltay = ymax - ymin;
	
	//Make the frame histogram
	TH2D *frame = new TH2D("hPulsesWfs_frame", ";Time (samples);Amplitude [ADC]", nBins, xlow, xup, 1000, ymin-0.1*deltay, ymax+0.1*deltay);
	
	
	TCanvas *c1 = (TCanvas*)gROOT->FindObject("canvPulsesWfs");
	if(c1) delete c1;
	c1 = new TCanvas("canvPulsesWfs", "All pulses waveforms", 1100, 600);
	
	frame->Draw();
	
	vector<TH1D*>::iterator vecIt;
	for(vecIt=colWfs.begin(); vecIt!=colWfs.end(); vecIt++)
	{
		(*vecIt)->Draw("same");
	}
	for(vecIt=indWfs.begin(); vecIt!=indWfs.end(); vecIt++)
	{
		(*vecIt)->Draw("same");
	}
}


void TestPulsesOverlap(Int_t event, UInt_t colwire, UInt_t indwire)
{//This is a routine mainly used to debug the RSTPC_RunProcessor::CombinePulses(...) routine
	
	if(!t2w)
	{
		t2w = new RSTPC_T2wrapper("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root", true);
	}
	
	if( !t2w->IsInit() ) return;
	
	if (!t2w->fChain) return;
	
	if( !((t2w->fT1wr) && (t2w->fT1wr->IsInit())) ) return;
	
	
	t2w->GetEntry(event);
	
	
	Int_t nChs = t2w->fT1wr->ColHist->GetNbinsY();
	Int_t nBins = t2w->fT1wr->ColHist->GetNbinsX();
	Double_t xlow = -0.5;
	Double_t xup = (nBins+1)-0.5;
	
	
	
	//===== Initialisation of the histograms =====//
	
	TH1D *ColWfOrig, *IndWfOrig;//Here the entire wf of the channel is present
	TH1D *ColWfAll, *IndWfAll; //Here there are all the pulses of the corresponding channels
	TH1D *ColWfCoin, *IndWfCoin;//Here only the overlapping pulses will be inserted
	
	
	stringstream hname;
	
	hname.str(""); hname << t2w->fT1wr->ColHist->GetName() << "_ch" << colwire;
	ColWfOrig = t2w->fT1wr->ColHist->ProjectionX(hname.str().c_str(), colwire+1, colwire+1);
	
	hname.str(""); hname << t2w->fT1wr->IndHist->GetName() << "_ch" << indwire;
	IndWfOrig = t2w->fT1wr->IndHist->ProjectionX(hname.str().c_str(), indwire+1, indwire+1);
	
	
	hname.str(""); hname << "hColPulCoin_ev" << event << "_ch" << colwire;
	ColWfCoin = (TH1D*)gROOT->FindObject( hname.str().c_str() );
	if(ColWfCoin)
	{
		ColWfCoin->Reset();
	}
	else
	{
		ColWfCoin = (TH1D*)ColWfOrig->Clone( hname.str().c_str() );
		ColWfCoin->Reset();
		ColWfCoin->SetLineColor(kBlue);
	}
	
	hname.str(""); hname << "hColPulAll_ev" << event << "_ch" << colwire;
	ColWfAll = (TH1D*)gROOT->FindObject( hname.str().c_str() );
	if(ColWfAll)
	{
		ColWfAll->Reset();
	}
	else
	{
		ColWfAll = (TH1D*)ColWfCoin->Clone( hname.str().c_str() );
		ColWfAll->Reset();
		ColWfAll->SetLineColor(kBlue);
	}
	
	
	
	hname.str(""); hname << "hIndPulCoin_ev" << event << "_ch" << indwire;
	IndWfCoin = (TH1D*)gROOT->FindObject( hname.str().c_str() );
	if(IndWfCoin)
	{
		IndWfCoin->Reset();
	}
	else
	{
		IndWfCoin = (TH1D*)IndWfOrig->Clone( hname.str().c_str() );
		IndWfCoin->Reset();
		IndWfCoin->SetLineColor(kRed);
	}
	
	hname.str(""); hname << "hIndPulAll_ev" << event << "_ch" << indwire;
	IndWfAll = (TH1D*)gROOT->FindObject( hname.str().c_str() );
	if(IndWfAll)
	{
		IndWfAll->Reset();
	}
	else
	{
		IndWfAll = (TH1D*)IndWfCoin->Clone( hname.str().c_str() );
		IndWfAll->Reset();
		IndWfAll->SetLineColor(kRed);
	}
	
	
	
	
	//===== Coincidence algorithm =====//
	
	TIter ColPulsesIt(t2w->ColPulses);
	TIter IndPulsesIt(t2w->IndPulses);
	RSTPC_Pulse *ColPulse, *IndPulse;
	
	
	map<RSTPC_Pulse*, vector<RSTPC_Pulse*> > CoinMap;
	Int_t nMatches = 0;
	while( (ColPulse = (RSTPC_Pulse*)ColPulsesIt.Next()) )
	{
		Bool_t match = false;
		
		if(!(ColPulse->fWireNum==colwire)) continue;
		
		for(Int_t cSamp=ColPulse->fLedge; cSamp<=ColPulse->fRedge; cSamp++)
		{
			ColWfAll->SetBinContent( cSamp+1, ColWfOrig->GetBinContent(cSamp+1) );
		}
		
		vector<RSTPC_Pulse*> IndPulsesVec;
		
		while( (IndPulse = (RSTPC_Pulse*)IndPulsesIt.Next()) )
		{
			if(!(IndPulse->fWireNum==indwire)) continue;
			
			for(Int_t iSamp=IndPulse->fLedge; iSamp<=IndPulse->fRedge; iSamp++)
			{
				IndWfAll->SetBinContent( iSamp+1, IndWfOrig->GetBinContent(iSamp+1) );
			}
			
			if( !( (ColPulse->fLedge>IndPulse->fRedge)||(ColPulse->fRedge<IndPulse->fLedge) ) )
			{//This is the overlapping condition
				IndPulsesVec.push_back(IndPulse);
				match = true;
				nMatches++;
			}
		}
		if(match) CoinMap[ColPulse] = IndPulsesVec;
	}
	
	
	cout << "\nFound in total " << nMatches << " matching conditions:" << endl;
	
	map<RSTPC_Pulse*, vector<RSTPC_Pulse*> >::iterator mapIt;
	for(mapIt=CoinMap.begin(); mapIt!=CoinMap.end(); mapIt++)
	{
		ColPulse = mapIt->first;
		vector<RSTPC_Pulse*> IndPulsesVec = mapIt->second;
		vector<RSTPC_Pulse*>::iterator vecIt;
		for(vecIt=IndPulsesVec.begin(); vecIt!=IndPulsesVec.end(); vecIt++)
		{
			IndPulse = (*vecIt);
			cout << "Coin: Cpulse (" << ColPulse->fLedge << "," << ColPulse->fRedge << "); Ipulse (" << IndPulse->fLedge << "," << IndPulse->fRedge << ")." << endl;
			
			//Draw also the pulses in the respective histograms
			for(Int_t cSamp=ColPulse->fLedge; cSamp<=ColPulse->fRedge; cSamp++)
			{
				ColWfCoin->SetBinContent( cSamp+1, ColWfOrig->GetBinContent(cSamp+1) );
			}
			
			for(Int_t iSamp=IndPulse->fLedge; iSamp<=IndPulse->fRedge; iSamp++)
			{
				IndWfCoin->SetBinContent( iSamp+1, IndWfOrig->GetBinContent(iSamp+1) );
			}
		}
	}
	
	
	TCanvas *canvTestPulsesOverlap = (TCanvas*)gROOT->FindObject("canvTestPulsesOverlap");
	if(canvTestPulsesOverlap) delete canvTestPulsesOverlap;
	canvTestPulsesOverlap = new TCanvas("canvTestPulsesOverlap", "Overlapping pulses",1020,640);
	canvTestPulsesOverlap->Divide(1,2);
	
	stringstream htitle;
	
	
	canvTestPulsesOverlap->cd(1);
	
	htitle.str(""); htitle << "Col wire " << colwire << ", Ind wire " << indwire << " - All pulses; Time [samples]; Amplitude [A.U.]";
	ColWfAll->SetTitle( htitle.str().c_str() ); ColWfAll->SetLineColor(kBlue);
	ColWfAll->Draw();
	
	IndWfAll->SetTitle( htitle.str().c_str() ); IndWfAll->SetLineColor(kRed);
	IndWfAll->Draw("same");
	
	
	canvTestPulsesOverlap->cd(2);
	
	htitle.str(""); htitle << "Col wire " << colwire << ", Ind wire " << indwire << " - Overlapping pulses; Time [samples]; Amplitude [A.U.]";
	ColWfCoin->SetTitle( htitle.str().c_str() ); ColWfCoin->SetLineColor(kBlue);
	ColWfCoin->Draw();
	
	IndWfCoin->SetTitle( htitle.str().c_str() ); IndWfCoin->SetLineColor(kRed);
	IndWfCoin->Draw("same");
}


void TestPulsesOverlap(Int_t event)
{//This makes just the print out of all the coincidences
	
	if(!t2w)
	{
		t2w = new RSTPC_T2wrapper("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root", true);
	}
	
	if( !t2w->IsInit() ) return;
	
	if (!t2w->fChain) return;
	
	if( !((t2w->fT1wr) && (t2w->fT1wr->IsInit())) ) return;
	
	
	t2w->GetEntry(event);
	
	
	Int_t nChs = t2w->fT1wr->ColHist->GetNbinsY();
	Int_t nBins = t2w->fT1wr->ColHist->GetNbinsX();
	Double_t xlow = -0.5;
	Double_t xup = (nBins+1)-0.5;
	
	
	
	//===== Coincidence algorithm =====//
	
	TIter ColPulsesIt(t2w->ColPulses);
	TIter IndPulsesIt(t2w->IndPulses);
	RSTPC_Pulse *ColPulse, *IndPulse;
	
	
	map<RSTPC_Pulse*, vector<RSTPC_Pulse*> > CoinMap;
	Int_t nMatches = 0;
	while( (ColPulse = (RSTPC_Pulse*)ColPulsesIt.Next()) )
	{
		Bool_t match = false;
		
		vector<RSTPC_Pulse*> IndPulsesVec;
		
		while( (IndPulse = (RSTPC_Pulse*)IndPulsesIt.Next()) )
		{
			if( !( (ColPulse->fLedge>IndPulse->fRedge)||(ColPulse->fRedge<IndPulse->fLedge) ) )
			{//This is the overlapping condition
				IndPulsesVec.push_back(IndPulse);
				match = true;
				nMatches++;
			}
		}
		if(match) CoinMap[ColPulse] = IndPulsesVec;
	}
	
	
	cout << "\nFound in total " << nMatches << " matching conditions:" << endl;
	
	map<RSTPC_Pulse*, vector<RSTPC_Pulse*> >::iterator mapIt;
	for(mapIt=CoinMap.begin(); mapIt!=CoinMap.end(); mapIt++)
	{
		ColPulse = mapIt->first;
		vector<RSTPC_Pulse*> IndPulsesVec = mapIt->second;
		vector<RSTPC_Pulse*>::iterator vecIt;
		for(vecIt=IndPulsesVec.begin(); vecIt!=IndPulsesVec.end(); vecIt++)
		{
			IndPulse = (*vecIt);
			cout << "  Cpulse wire " << ColPulse->fWireNum << " (" << ColPulse->fLedge << "," << ColPulse->fRedge << "); Ipulse wire " << IndPulse->fWireNum << " (" << IndPulse->fLedge << "," << IndPulse->fRedge << ")." << endl;
		}
	}
	
}


void PlotColPulses()
{
	if(!t2w)
	{
		t2w = new RSTPC_T2wrapper("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root", true);
	}
	
	if( !t2w->IsInit() ) return;
	
	if( !t2w->fChain ) return;
	
	if( !((t2w->fT1wr) && (t2w->fT1wr->IsInit())) ) return;
	
	
	Int_t nEvs = t2w->fChain->GetEntries();
	
	
	//Plot the ampitude and the widths of the noise collection pulses
	const Int_t mintimenoise = 2500;
	
	vector<Int_t> ledgesVec, widthsVec, maxposVec;
	vector<Double_t> maxampVec;
	
	
	for(Int_t iEv=0; iEv<nEvs; iEv++)
	{
		t2w->GetEntry(iEv);
		
		//Iterate over all the pulses but select only those of the collection wires
		RSTPC_Pulse *ColPulse;
		TIter ColPulsesIt(t2w->ColPulses);
		while( (ColPulse = (RSTPC_Pulse*)ColPulsesIt.Next()) )
		{
			if(ColPulse->fWireType!=kCol) continue;
			if(ColPulse->fLedge<mintimenoise) continue;
			
			Int_t width = ColPulse->fRedge - ColPulse->fLedge;
			maxposVec.push_back( ColPulse->fMaxPos );
			widthsVec.push_back( width );
			maxampVec.push_back( ColPulse->fMax );
		}
	}
	
	Double_t minWidth = (Double_t)TMath::MinElement( widthsVec.size(), &widthsVec.at(0) );
	Double_t maxWidth = (Double_t)TMath::MinElement( widthsVec.size(), &widthsVec.at(0) );
	Double_t deltaX = maxWidth - minWidth;
	minWidth -= 0.1*deltaX;
	maxWidth += 0.1*deltaX;
	
	Double_t minAmpl = TMath::MinElement( maxampVec.size(), &maxampVec.at(0) );
	Double_t maxAmpl = TMath::MinElement( maxampVec.size(), &maxampVec.at(0) );
	Double_t deltaY = maxAmpl - minAmpl;
	minAmpl -= 0.1*deltaY;
	maxAmpl += 0.1*deltaY;
	
	TH2D* hNoiseColPulses = (TH2D*)gROOT->FindObject("hNoiseColPulses");
	if(hNoiseColPulses) delete hNoiseColPulses;
	hNoiseColPulses = new TH2D("hNoiseColPulses","Collection pulses - Only noise region;Width [samples];Amplitude [AU]", 50, minWidth, maxWidth, 50, minAmpl, maxAmpl);
	
	for(Int_t iPulse=0; iPulse<maxposVec.size(); iPulse++)
	{
		hNoiseColPulses->Fill( widthsVec.at(iPulse), maxampVec.at(iPulse) );
	}
	
	TCanvas *canvPlotColPulses1 = (TCanvas*)gROOT->FindObject("canvPlotColPulses1");
	if(!canvPlotColPulses1)
	{
		canvPlotColPulses1 = new TCanvas("canvPlotColPulses1","Noise collection pulses", 800, 600);
	}
	canvPlotColPulses1->cd()->SetLogz();
	
	hNoiseColPulses->Draw("colz");
	
	
	
	
	//Now plot the same for the signal region
	const Double_t minSignTime = 200; //Time in samples units
	const Double_t maxSignTime = 1300; //Time in samples number
	
	maxposVec.clear();
	widthsVec.clear();
	maxampVec.clear();
	
	for(Int_t iEv=0; iEv<nEvs; iEv++)
	{
		t2w->GetEntry(iEv);
		
		//Iterate over all the pulses but select only those of the collection wires
		RSTPC_Pulse *ColPulse;
		TIter ColPulsesIt(t2w->ColPulses);
		while( (ColPulse = (RSTPC_Pulse*)ColPulsesIt.Next()) )
		{
			if(ColPulse->fWireType!=kCol) continue;
			if(ColPulse->fLedge<minSignTime) continue;
			if(ColPulse->fLedge>maxSignTime) continue;
			
			Int_t width = ColPulse->fRedge - ColPulse->fLedge;
			maxposVec.push_back( ColPulse->fMaxPos );
			widthsVec.push_back( width );
			maxampVec.push_back( ColPulse->fMax );
		}
	}
	
	minWidth = (Double_t)TMath::MinElement( widthsVec.size(), &widthsVec.at(0) );
	maxWidth = (Double_t)TMath::MinElement( widthsVec.size(), &widthsVec.at(0) );
	deltaX = maxWidth - minWidth;
	minWidth -= 0.1*deltaX;
	maxWidth += 0.1*deltaX;
	
	minAmpl = TMath::MinElement( maxampVec.size(), &maxampVec.at(0) );
	maxAmpl = TMath::MinElement( maxampVec.size(), &maxampVec.at(0) );
	deltaY = maxAmpl - minAmpl;
	minAmpl -= 0.1*deltaY;
	maxAmpl += 0.1*deltaY;
	
	TH2D* hSignalColPulses = (TH2D*)gROOT->FindObject("hSignalColPulses");
	if(hSignalColPulses) delete hSignalColPulses;
	hSignalColPulses = new TH2D("hSignalColPulses","Collection pulses - Only signal region;Width [samples];Amplitude [AU]", 500, minWidth, maxWidth, 500, minAmpl, maxAmpl);
	
	for(Int_t iPulse=0; iPulse<maxposVec.size(); iPulse++)
	{
		hSignalColPulses->Fill( widthsVec.at(iPulse), maxampVec.at(iPulse) );
	}
	
	TCanvas *canvPlotColPulses2 = (TCanvas*)gROOT->FindObject("canvPlotColPulses2");
	if(!canvPlotColPulses2)
	{
		canvPlotColPulses2 = new TCanvas("canvPlotColPulses2","Signal collection pulses", 800, 600);
	}
	canvPlotColPulses2->cd()->SetLogz();
	
	hSignalColPulses->Draw("colz");
	
}