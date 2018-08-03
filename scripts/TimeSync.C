#include "MppcTreeWrapper.hh"
#include "RSTPC_Analyser.hh"
#include "HistoManipulators.hh"
#include "DigitalFilters.hh"

#include <set>


RSTPC_Analyser *an = NULL;
MppcTreeWrapper *tw = NULL;

//Stuff that goes in the c1 canvas
TGraph *gr1 = NULL;
TH1D *hTimeDiffs = NULL;

//Stuff that goes in the c2 canvas
TGraph *gr2 = NULL;
TH1D* hDeltaTimesDiff = NULL;

//Stuff that goes in the c3 canvas
TGraph *gr3 = NULL;
TGraph *gr4 = NULL;


void TimeSync_1()
{//This is the first algorithm
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell/");
	
	if(!an)
	{
		an = new RSTPC_Analyser;
		an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
		an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
		//an->OpenRun(2026); //Peaking time 2 usec, 15 kV
		//an->OpenRun(2027); //Peaking time 2 usec, 3 kV
		//an->OpenRun(2028); //Peaking time 2 usec, 3 kV
		//an->OpenRun(2031); //Peaking time 1 usec, 15 kV
		an->OpenRun(2032); //Peaking time 1 usec, 20 kV
		//an->OpenRun(2033); //Peaking time 1 usec, 30 kV
	}
	
	
	
	if(!(an->fTrigTree)) return;
	ULong64_t EventTime;
	an->fTrigTree->SetBranchAddress("EventTime", &EventTime);
	Int_t nTriggers = an->fTrigTree->GetEntries();
	
	
	
	if(!tw)
	{
		tw = new MppcTreeWrapper("/home/francescop/data/ResistiveShell/July_run_2032_Cosmics.root");
		if(!tw->fChain)
		{
			delete tw;
			tw = NULL;
			return;
		}
	}
	
	Int_t nEvents = tw->fChain->GetEntries();
	
	
	
	vector<Double_t> TrigTimesVec(nTriggers);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		an->fTrigTree->GetEntry(iTrig);
		TrigTimesVec[iTrig] = ((Double_t)EventTime)/1e6;
	}
	
	vector<Double_t> FebTimesVec(nEvents);
	set<Int_t> FebFreeEvs;
	for(int iEv=0; iEv<nEvents; iEv++)
	{
		tw->fChain->GetEntry(iEv);
		FebTimesVec[iEv] = tw->GetEventTime();
		FebFreeEvs.insert(iEv);
	}
	
	
	
	map<Int_t, Int_t> FebEvsMap;
	
	//Cycle over the triggered events
	Int_t iEvSelected = -1;
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		Double_t mindiff, minabsdiff;
		
		set<Int_t>::iterator iEv, selEv;
		for(iEv=FebFreeEvs.begin(); iEv!=FebFreeEvs.end(); ++iEv)
		{
			Double_t diff = ( TrigTimesVec[iTrig]-FebTimesVec[*iEv] );
			if( iEv==FebFreeEvs.begin() )
			{
				mindiff = diff;
				minabsdiff = TMath::Abs(diff);
				selEv = iEv;
			}
			else if( TMath::Abs(diff)<mindiff )
			{
				mindiff = diff;
				minabsdiff = TMath::Abs(diff);
				selEv = iEv;
			}
		}//End cycle over the FEB events
		
		//TimeDiffVec[iTrig] = mindiff;
		FebEvsMap[iTrig] = *selEv;
		FebFreeEvs.erase(selEv);
		
	}//End of cycle over the triggered events
	
	
	vector<Double_t> SelFebTimesVec(nTriggers);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		SelFebTimesVec.at(iTrig) = FebTimesVec.at( FebEvsMap[iTrig] );
	}
	
	
	vector<Double_t> TimeDiffVec(nTriggers);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		TimeDiffVec.at(iTrig) = SelFebTimesVec.at(iTrig) - TrigTimesVec.at(iTrig);
	}
	
	
	
	
	TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
	if(!c1)
	{
		c1 = new TCanvas("c1","c1",1000,600);
		c1->Divide(2,1);
	}
	
	
	//Make and fill an histogram with the minimum times of the 
	Double_t minDiff = TMath::MinElement(nTriggers, &TimeDiffVec.at(0) );
	Double_t maxDiff = TMath::MaxElement(nTriggers, &TimeDiffVec.at(0) );
	Double_t deltaDiff = maxDiff-minDiff;
	minDiff = minDiff - 0.1*deltaDiff;
	maxDiff = maxDiff + 0.1*deltaDiff;
	
	
	c1->cd(1);
	
	if(gr1) delete gr1;
	gr1 = new TGraph;
	gr1->SetMarkerStyle(7);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		gr1->SetPoint(iTrig, iTrig, TimeDiffVec[iTrig]);
	}
	gr1->SetTitle(";Event number;Time_{FEB} - Time_{TPC} [s]");
	gr1->Draw("AP");
	
	
	c1->cd(2);
	
	hTimeDiffs = (TH1D*)gROOT->FindObject("hTimeDiffs");
	if(hTimeDiffs) delete hTimeDiffs;
	hTimeDiffs = new TH1D("hTimeDiffs", ";Time_{FEB} - Time_{TPC} [s]; Counts", 100, minDiff, maxDiff);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		hTimeDiffs->Fill(TimeDiffVec[iTrig]);
	}
	hTimeDiffs->Draw();
	
	
	
	TCanvas *c2 = (TCanvas*)gROOT->FindObject("c2");
	if(!c2)
	{
		c2 = new TCanvas("c2","c2",1000,600);
		c2->Divide(2,1);
	}
	
	
	//Make the vector of the FEB-TPC differences of the time differences between two consecutive TPC events
	vector<Double_t> DeltaTimesDiffVec;
	for(Int_t iTrig=0; iTrig<nTriggers-1; iTrig++)
	{
		Double_t val = (SelFebTimesVec[iTrig+1]-SelFebTimesVec[iTrig]) - (TrigTimesVec[iTrig+1]-TrigTimesVec[iTrig]) ;
		DeltaTimesDiffVec.push_back(val);
	}
	
	Double_t minDeltaTimeDiff = TMath::MinElement(DeltaTimesDiffVec.size(), &DeltaTimesDiffVec.at(0));
	Double_t maxDeltaTimeDiff = TMath::MaxElement(DeltaTimesDiffVec.size(), &DeltaTimesDiffVec.at(0));
	
	
	c2->cd(1);
	
	if(gr2) delete gr2;
	gr2 = new TGraph;
	gr2->SetMarkerStyle(7);
	for(Int_t iTrig=0; iTrig<DeltaTimesDiffVec.size(); iTrig++)
	{
		//gr->SetPoint(iTrig, iTrig, FEBevs[iTrig]);
		//gr->SetPoint(iTrig, iTrig, FEBTimes[FEBevs[iTrig]]);
		gr2->SetPoint(iTrig, iTrig, DeltaTimesDiffVec[iTrig]);
	}
	gr2->SetTitle(";Event number;#DeltaTime_{FEB} - #DeltaTime_{TPC} [s]");
	gr2->Draw("AP");
	
	
	c2->cd(2);
	
	hDeltaTimesDiff = (TH1D*)gROOT->FindObject("hDeltaTimesDiff");
	if(hDeltaTimesDiff) delete hDeltaTimesDiff;
	hDeltaTimesDiff = new TH1D("hDeltaTimesDiff", "TPC events;#DeltaTime_{FEB} - #DeltaTime_{TPC} [s];Counts", 100, minDeltaTimeDiff, maxDeltaTimeDiff);
	for(Int_t iTrig=0; iTrig<DeltaTimesDiffVec.size(); iTrig++)
	{
		hDeltaTimesDiff->Fill( DeltaTimesDiffVec[iTrig] );
	}
	hDeltaTimesDiff->SetLineColor(kBlue);
	hDeltaTimesDiff->Draw();
	
	
	
	TCanvas *c3 = new TCanvas("c3","c3",1000,800);
	c3->Divide(1,2);
	
	
	c3->cd(1);
	
	if(gr3) delete gr3;
	gr3 = new TGraph;
	gr3->SetMarkerStyle(7);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		gr3->SetPoint(iTrig, iTrig, FebEvsMap[iTrig]);
	}
	gr3->SetTitle(";TPC event number;FEB event number");
	gr3->Draw("AP");
	
	
	c3->cd(2);
	
	if(gr4) delete gr4;
	gr4 = new TGraph;
	gr4->SetMarkerStyle(7);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		gr4->SetPoint(iTrig, iTrig, SelFebTimesVec.at(iTrig));
	}
	gr4->SetTitle(";TPC event number;FEB time [s]");
	gr4->Draw("AP");
	
}


void TimeSync_2()
{//This is the first algorithm: here the selection is based on the time differences 
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell/");
	
	if(!an)
	{
		an = new RSTPC_Analyser;
		an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
		an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
		//an->OpenRun(2026); //Peaking time 2 usec, 15 kV
		//an->OpenRun(2027); //Peaking time 2 usec, 3 kV
		//an->OpenRun(2028); //Peaking time 2 usec, 3 kV
		//an->OpenRun(2031); //Peaking time 1 usec, 15 kV
		an->OpenRun(2032); //Peaking time 1 usec, 20 kV
		//an->OpenRun(2033); //Peaking time 1 usec, 30 kV
	}
	
	
	//Set the max tollerated time difference between FEB and TPC events
	const Double_t MaxTimeDiff = 1.0; //In seconds
	
	
	if(!(an->fTrigTree)) return;
	ULong64_t EventTime;
	an->fTrigTree->SetBranchAddress("EventTime", &EventTime);
	Int_t nTriggers = an->fTrigTree->GetEntries();
	
	
	
	if(!tw)
	{
		tw = new MppcTreeWrapper("/home/francescop/data/ResistiveShell/July_run_2032_Cosmics.root");
		if(!tw->fChain)
		{
			delete tw;
			tw = NULL;
			return;
		}
	}
	
	Int_t nEvents = tw->fChain->GetEntries();
	
	
	vector<Double_t> TrigTimesVec(nTriggers);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		an->fTrigTree->GetEntry(iTrig);
		TrigTimesVec.at(iTrig) = ((Double_t)EventTime)/1e6;
	}
	
	//This will contain the time differences between 2 consecutive events. It has 1 element less than the vector of the times
	map<Int_t, Double_t> TrigTimesDiffMap; 
	for(Int_t iTrig=0; iTrig<nTriggers-1; iTrig++)
	{
		TrigTimesDiffMap[iTrig] = TrigTimesVec[iTrig+1]-TrigTimesVec[iTrig];
	}
	
	
	set<Int_t> FebFreeEvs;
	vector<Double_t> FebTimesVec(nEvents);
	for(int iEv=0; iEv<nEvents; iEv++)
	{
		tw->fChain->GetEntry(iEv);
		FebTimesVec[iEv] = tw->GetEventTime();
		FebFreeEvs.insert(iEv);
	}
	
	map<Int_t,Int_t> FebEvsMap; //The key corresponds to the TPC event the second to the corresponding FEB event
	
	//Cycle over the triggered events
	Int_t iEvSelected = -1;
	for(Int_t iTrig=0; iTrig<nTriggers-1; iTrig++)
	{
		
		Double_t dTpcDiffTime = TrigTimesDiffMap[iTrig]; //This is the goal time difference between two consecutive Feb events
		Double_t dMinTime = TrigTimesVec[iTrig] - MaxTimeDiff;
		Double_t dMaxTime = TrigTimesVec[iTrig+1] + MaxTimeDiff;
		
		vector<Int_t> candEvs; candEvs.clear();
		
		set<Int_t>::iterator iEvIt;
		
		Int_t lowIdx, upIdx; //Indexes of the Feb events
		if(iTrig==0)
		{//Here I have to determine both the low index and the up index
			
			//Select all the candidate FEB events that have their time within the max allowed time difference
			for(iEvIt=FebFreeEvs.begin(); iEvIt!=FebFreeEvs.end(); iEvIt++)
			{
				Double_t dFebTime = FebTimesVec[*iEvIt];
				if( (dFebTime>=dMinTime)&&(dFebTime<=dMaxTime) ) candEvs.push_back(*iEvIt);
			}
			
			Int_t nCandEvs = candEvs.size();
			
			Double_t deltadiff; //Absolute difference of time differences
			for(Int_t iEv=0; iEv<nCandEvs; iEv++)
			{
				for(Int_t jEv=iEv+1; jEv<nCandEvs; jEv++)
				{
					if( (iEv==0)&&(jEv==1) )
					{
						lowIdx = candEvs[0];
						upIdx = candEvs[1];
						deltadiff = abs( dTpcDiffTime - (FebTimesVec[candEvs[1]]-FebTimesVec[candEvs[0]]) );
					}
					else if( deltadiff>abs( dTpcDiffTime - (FebTimesVec[candEvs[jEv]]-FebTimesVec[candEvs[iEv]]) ) )
					{
						if( (FebTimesVec[candEvs[jEv]]-FebTimesVec[candEvs[iEv]])>0. )
						{//This is redundant but better to keep
							lowIdx = candEvs[iEv];
							upIdx = candEvs[jEv];
							deltadiff = abs( dTpcDiffTime - (FebTimesVec[candEvs[jEv]]-FebTimesVec[candEvs[iEv]]) );
						}
					}
				}
			}
		
			//Here I have selected the two events FEB events that can correspond to the TPC event iTrig and iTrig+1
			FebEvsMap[iTrig] = lowIdx;
			FebEvsMap[iTrig+1] = upIdx;
		}
		else
		{//Here the low index is the up index of the previous TPC event and only the up index has to searched for
			lowIdx = upIdx;
			
			//Select all the candidate FEB events that have their time within the max allowed time difference
			for(iEvIt=FebFreeEvs.begin(); iEvIt!=FebFreeEvs.end(); iEvIt++)
			{
				Double_t dFebTime = FebTimesVec[*iEvIt];
				if( (dFebTime>=dMinTime) && (dFebTime<=dMaxTime) ) candEvs.push_back(*iEvIt);
			}
			
			Int_t nCandEvs = candEvs.size();
			
			Double_t deltadiff; //Absolute difference of time differences
			
			for(Int_t iEv=0; iEv<nCandEvs; iEv++)
			{
				if(iEv==0)
				{
					upIdx = candEvs[iEv];
					deltadiff = abs( dTpcDiffTime - (FebTimesVec[candEvs[iEv]]-FebTimesVec[lowIdx]) );
				}
				else if(deltadiff < abs( dTpcDiffTime - (FebTimesVec[candEvs[iEv]]-FebTimesVec[lowIdx]) ))
				{
					if( (FebTimesVec[candEvs[iEv]]-FebTimesVec[lowIdx])>0. )
					{//This is redundant but better to keep
						upIdx = candEvs[iEv];
						deltadiff = abs( dTpcDiffTime - (FebTimesVec[candEvs[iEv]]-FebTimesVec[lowIdx]) );
					}
					
				}
			}
			
			FebEvsMap[iTrig+1] = upIdx;
		}
		
		
		//Remove from the FEB free events all those that have a time lower than FEBTimesVec[upIdx]
		//Must be done in 2 steps as I cannot change the elements in the FEBfreeEvs while iterating over it
		set<Int_t> delList; delList.clear();
		FebFreeEvs.erase(lowIdx);
		FebFreeEvs.erase(upIdx);
		for(iEvIt=FebFreeEvs.begin(); iEvIt!=FebFreeEvs.end(); iEvIt++)
		{
			Double_t dFebTime = FebTimesVec[*iEvIt];
			if( dFebTime<FebTimesVec[upIdx] ) delList.insert(*iEvIt);
		}
		
		for(iEvIt=delList.begin(); iEvIt!=delList.end(); iEvIt++)
		{
			FebFreeEvs.erase(*iEvIt);
		}
		
		
	}//End of cycle over the triggered events
	
	vector<Double_t> TimeDiffVec(nTriggers);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		Int_t febIdx = FebEvsMap[iTrig];
		TimeDiffVec.at(iTrig) = FebTimesVec.at(febIdx)-TrigTimesVec.at(iTrig);
	}
	
	
	
	TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
	if(!c1)
	{
		c1 = new TCanvas("c1","c1",1280,640);
		c1->Divide(2,1);
	}
	
	//Make and fill an histogram with the minimum times of the 
	Double_t minDiff = TMath::MinElement( nTriggers, &TimeDiffVec.at(0) );
	Double_t maxDiff = TMath::MaxElement( nTriggers, &TimeDiffVec.at(0) );
	Double_t deltaDiff = maxDiff-minDiff;
	minDiff = minDiff - 0.1*deltaDiff;
	maxDiff = maxDiff + 0.1*deltaDiff;
	
	
	c1->cd(1);
	
	if(gr1) delete gr1;
	gr1 = new TGraph;
	gr1->SetMarkerStyle(7);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		//gr->SetPoint(iTrig, iTrig, FEBevs[iTrig]);
		//gr->SetPoint(iTrig, iTrig, FEBTimes[FEBevs[iTrig]]);
		gr1->SetPoint(iTrig, iTrig, TimeDiffVec[iTrig]);
	}
	gr1->SetTitle(";Event number;Time_{FEB} - Time_{TPC} [s]");
	gr1->Draw("AP");
	
	
	c1->cd(2);
	
	hTimeDiffs = (TH1D*)gROOT->FindObject("hTimeDiffs");
	if(hTimeDiffs) delete hTimeDiffs;
	hTimeDiffs = new TH1D("hTimeDiffs", ";Time_{FEB} - Time_{TPC} [s]; Counts", 100, minDiff, maxDiff);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		hTimeDiffs->Fill( TimeDiffVec.at(iTrig) );
	}
	hTimeDiffs->Draw();
	
	
	
	TCanvas *c2 = (TCanvas*)gROOT->FindObject("c2");
	if(!c2)
	{
		c2 = new TCanvas("c2","c2",1280,640);
		c2->Divide(2,1);
	}
	
	//Make the vector of the FEB-TPC differences of the time differences between two consecutive TPC events
	vector<Double_t> DeltaTimesDiffVec;
	for(Int_t iTrig=0; iTrig<nTriggers-1; iTrig++)
	{
		Double_t val = (FebTimesVec[FebEvsMap[iTrig+1]]-FebTimesVec[FebEvsMap[iTrig]]) - (TrigTimesVec[iTrig+1]-TrigTimesVec[iTrig]) ;
		DeltaTimesDiffVec.push_back(val);
	}
	
	Double_t minDeltaTimeDiff = TMath::MinElement(DeltaTimesDiffVec.size(), &DeltaTimesDiffVec.at(0));
	Double_t maxDeltaTimeDiff = TMath::MaxElement(DeltaTimesDiffVec.size(), &DeltaTimesDiffVec.at(0));
	
	
	c2->cd(1);
	
	if(gr2) delete gr2;
	gr2 = new TGraph;
	gr2->SetMarkerStyle(7);
	for(Int_t iTrig=0; iTrig<DeltaTimesDiffVec.size(); iTrig++)
	{
		//gr->SetPoint(iTrig, iTrig, FEBevs[iTrig]);
		//gr->SetPoint(iTrig, iTrig, FEBTimes[FEBevs[iTrig]]);
		gr2->SetPoint(iTrig, iTrig, DeltaTimesDiffVec[iTrig]);
	}
	gr2->SetTitle(";Event number;#DeltaTime_{FEB} - #DeltaTime_{TPC} [s]");
	gr2->Draw("AP");
	
	
	c2->cd(2);
	
	TH1D* hDeltaTimesDiff = (TH1D*)gROOT->FindObject("hDeltaTimesDiff");
	if(hDeltaTimesDiff) delete hDeltaTimesDiff;
	hDeltaTimesDiff = new TH1D("hDeltaTimesDiff", ";#DeltaTime_{FEB} - #DeltaTime_{TPC} [s];Counts", 100, minDeltaTimeDiff, maxDeltaTimeDiff);
	for(Int_t iTrig=0; iTrig<DeltaTimesDiffVec.size(); iTrig++)
	{
		hDeltaTimesDiff->Fill( DeltaTimesDiffVec[iTrig] );
	}
	hDeltaTimesDiff->SetLineColor(kBlue);
	hDeltaTimesDiff->Draw();
	
	
	
	TCanvas *c3 = (TCanvas*)gROOT->FindObject("c3");
	if(!c3)
	{
		c3 = new TCanvas("c3","c3",1280,640);
		c3->Divide(1,2);
	}
	
	
	c3->cd(1);
	
	if(gr3) delete gr3;
	gr3 = new TGraph;
	gr3->SetMarkerStyle(7);
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		gr3->SetPoint(iTrig, iTrig, FebEvsMap[iTrig]);
	}
	gr3->SetTitle(";TPC event number;FEB event number");
	gr3->Draw("AP");
	
	
	c3->cd(2);
	
	if(gr4) delete gr4;
	gr4 = new TGraph;
	gr4->SetMarkerStyle(7);
	
	for(Int_t iTrig=0; iTrig<nTriggers; iTrig++)
	{
		gr4->SetPoint(iTrig, iTrig, FebTimesVec[FebEvsMap[iTrig]]);
	}
	gr4->SetTitle(";TPC event number;FEB time [s]");
	gr4->Draw("AP");
	
}


void ComparePlots()
{
	TFile *infile = TFile::Open("SyncHists.root");
	
	
	
	
	TGraph *gr1_alg1 = (TGraph*)infile->Get("gr1_alg1");
	if(gr1_alg1) gr1_alg1->SetMarkerColor(kBlue);
	TGraph *gr1_alg2 = (TGraph*)infile->Get("gr1_alg2");
	if(gr1_alg2) gr1_alg2->SetMarkerColor(kRed);
	
	TMultiGraph *mgr1 = new TMultiGraph;
	mgr1->Add(gr1_alg1, "P");
	mgr1->Add(gr1_alg2, "P");
	mgr1->SetTitle(";Event number;Time_{FEB} - Time_{TPC} [s]");
	
	TH1D* hTimeDiffs_alg1 = (TH1D*)infile->Get("hTimeDiffs_alg1");
	if(hTimeDiffs_alg1) hTimeDiffs_alg1->SetLineColor(kBlue);
	TH1D* hTimeDiffs_alg2 = (TH1D*)infile->Get("hTimeDiffs_alg2");
	if(hTimeDiffs_alg2) hTimeDiffs_alg2->SetLineColor(kRed);
	
	TCanvas *c1 = new TCanvas("c1","c1", 1100, 500);
	c1->Divide(2,1);
	
	c1->cd(1);
	mgr1->Draw("AP");
	
	c1->cd(2);
	hTimeDiffs_alg1->Draw();
	hTimeDiffs_alg2->Draw("same");
	
	TLegend *leg = new TLegend(0.9,0.65,0.75,0.9);
	leg->AddEntry(hTimeDiffs_alg1, "Algorithm 1", "L");
	leg->AddEntry(hTimeDiffs_alg2, "Algorithm 2", "L");
	leg->Draw();
	
	
	
	TGraph *gr2_alg1 = (TGraph*)infile->Get("gr2_alg1");
	if(gr2_alg1) gr2_alg1->SetMarkerColor(kBlue);
	TGraph *gr2_alg2 = (TGraph*)infile->Get("gr2_alg2");
	if(gr2_alg2) gr2_alg2->SetMarkerColor(kRed);
	
	TMultiGraph *mgr2 = new TMultiGraph;
	mgr2->Add(gr2_alg1, "P");
	mgr2->Add(gr2_alg2, "P");
	mgr2->SetTitle(";Event number;#DeltaTime_{FEB} - #DeltaTime_{TPC} [s]");
	
	TH1D* hDeltaTimesDiff_alg1 = (TH1D*)infile->Get("hDeltaTimesDiff_alg1");
	if(hDeltaTimesDiff_alg1) hDeltaTimesDiff_alg1->SetLineColor(kBlue);
	TH1D* hDeltaTimesDiff_alg2 = (TH1D*)infile->Get("hDeltaTimesDiff_alg2");
	if(hDeltaTimesDiff_alg2) hDeltaTimesDiff_alg2->SetLineColor(kRed);
	
	TCanvas *c2 = new TCanvas("c2","c2", 1100, 500);
	c2->Divide(2,1);
	
	c2->cd(1);
	mgr2->Draw("AP");
	
	c2->cd(2);
	hDeltaTimesDiff_alg1->Draw();
	hDeltaTimesDiff_alg2->Draw("same");
	
	leg->Draw();
	
	
	
	TGraph *gr3_alg1 = (TGraph*)infile->Get("gr3_alg1");
	if(gr3_alg1) gr3_alg1->SetMarkerColor(kBlue);
	TGraph *gr3_alg2 = (TGraph*)infile->Get("gr3_alg2");
	if(gr3_alg2) gr3_alg2->SetMarkerColor(kRed);
	
	TMultiGraph *mgr3 = new TMultiGraph;
	mgr3->Add(gr3_alg1, "P");
	mgr3->Add(gr3_alg2, "P");
	mgr3->SetTitle(";TPC event number;FEB event number");
	
	TGraph *gr4_alg1 = (TGraph*)infile->Get("gr4_alg1");
	if(gr4_alg1) gr4_alg1->SetMarkerColor(kBlue);
	TGraph *gr4_alg2 = (TGraph*)infile->Get("gr4_alg2");
	if(gr4_alg2) gr4_alg2->SetMarkerColor(kRed);
	
	TMultiGraph *mgr4 = new TMultiGraph;
	mgr4->Add(gr4_alg1, "P");
	mgr4->Add(gr4_alg2, "P");
	mgr4->SetTitle(";TPC event number;FEB time [s]");
	
	TCanvas *c3 = new TCanvas("c3","c3",1000,800);
	c3->Divide(1,2);
	
	c3->cd(1);
	mgr3->Draw("AP");
	
	c3->cd(2);
	mgr4->Draw("AP");
}



Double_t *febtimes_arr = NULL;
Int_t CheckFebTimesConsistency()
{
	Int_t ret = 0;
	
	MppcTreeWrapper *twx = new MppcTreeWrapper("/home/francescop/data/ResistiveShell/July_run_2031_Cosmics.root");
	
	if(!twx->fChain) return -1;
	Int_t nEvents = twx->fChain->GetEntriesFast();
	
	if(febtimes_arr)
	{
		delete [] febtimes_arr;
	}
	febtimes_arr = new Double_t[nEvents];
	
	for(Int_t iEv=0; iEv<nEvents; iEv++)
	{
		twx->GetEntry(iEv);
		febtimes_arr[iEv] = twx->GetEventTime();
	}
	delete twx;
	
	
	for(Int_t iEv=0; iEv<nEvents-1; iEv++)
	{
		if(febtimes_arr[iEv]>febtimes_arr[iEv+1])
		{
			std::cout << "Suspicious events: (" << iEv << ": " << febtimes_arr[iEv] << ") and (" << iEv+1 << ": " << febtimes_arr[iEv+1] << ")" << endl;
			ret++;
		}
	}
	
	return ret;
}



ULong64_t *tpctimes_arr = NULL;
Int_t CheckTpcTimesConsistency()
{
	Int_t ret = 0;
	
	RSTPC_Analyser *anx = new RSTPC_Analyser;
	anx->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
	anx->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
	anx->OpenRun(2031);
	
	if(!anx->fTrigTree)
	{
		delete anx;
		return -1;
	}
	ULong64_t EventTime;
	anx->fTrigTree->SetBranchAddress("EventTime", &EventTime);
	
	Int_t nEvents = anx->fTrigTree->GetEntries();
	
	if(tpctimes_arr)
	{
		delete [] tpctimes_arr;
	}
	tpctimes_arr = new ULong64_t[nEvents];
	
	for(Int_t iEv=0; iEv<nEvents; iEv++)
	{
		anx->fTrigTree->GetEntry(iEv);
		tpctimes_arr[iEv] = EventTime;
	}
	if(anx) delete anx;
	
	
	for(Int_t iEv=0; iEv<nEvents-1; iEv++)
	{
		if(tpctimes_arr[iEv]>tpctimes_arr[iEv+1])
		{
			std::cout << "Suspicious events: (" << iEv << ": " << tpctimes_arr[iEv] << ") and (" << iEv+1 << ": " << tpctimes_arr[iEv+1] << ")" << endl;
			ret++;
		}
	}
	
	return ret;
}