#include "HistoManipulators.hh"
#include "DigitalFilters.hh"

#include "RSTPC_RunProcessor.hh"
#include "RSTPC_Analyser.hh"
#include "MppcTreeWrapper.hh"
#include "RSTPC_T1wrapper.hh"
#include "RSTPC_T2wrapper.hh"
#include "RSTPC_Hits.hh"

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

// Header file for the classes stored in the TTree if any.
#include "TObjArray.h"
#include "TClassTable.h"

#include <vector>
#include <map>

using namespace std;



RSTPC_T2wrapper *t2w = NULL;


TClonesArray *gHitsArr = NULL;


class CrossCorrRes
{
public:
	Double_t max;
	Double_t taumax;
	Bool_t bad;
};

CrossCorrRes CalculatePulsesCrossCorrelation( RSTPC_Pulse* ColPulse, TH2D* ColHist, RSTPC_Pulse* IndPulse, TH2D* IndHist );

/*
Int_t HitsFinderSlices(Int_t event)
{
	if(!t2w)
	{
		t2w = new RSTPC_T2wrapper("/home/francescop/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root", true);
	}
	
	if( !t2w->IsInit() ) return -1;
	
	
	RSTPC_RunProcessor::SetPeakingTime(1.0);
	RSTPC_RunProcessor::SetDriftVel(2.06);
	
	t2w->GetEntry(event);
	
	vector<RSTPC_Pulse*> ColPulsesVec;
	vector<RSTPC_Pulse*> IndPulsesVec;
	
	TIter ColPulsesIt(t2w->ColPulses);
	TIter IndPulsesIt(t2w->IndPulses);
	
	RSTPC_Pulse *ColPulse=NULL, *IndPulse=NULL;
	
	while( (ColPulse = (RSTPC_Pulse*)ColPulsesIt.Next()) )
	{
		ColPulsesVec.push_back(ColPulse);
	}
	
	while( (IndPulse = (RSTPC_Pulse*)IndPulsesIt.Next()) )
	{
		IndPulsesVec.push_back(IndPulse);
	}
	
	
	map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >* fPulseMap = RSTPC_RunProcessor::CombinePulses( &ColPulsesVec, &IndPulsesVec );
	
	
	
	Int_t nHits = 0;
	
	//Make time slices and find the hits
	Double_t gPeakTime = RSTPC_RunProcessor::GetPeakingTime();
	Double_t gSamplFreq = RSTPC_RunProcessor::GetSamplingFreq();
	Double_t gDriftLen = RSTPC_RunProcessor::GetDriftLenght();
	Double_t gDriftVel = RSTPC_RunProcessor::GetDriftVel();
	
	Int_t timeSlice = (Int_t)( (gPeakTime*gSamplFreq) + 0.5 ); //Trick to round to the closest integer
	Int_t nSlices = ceil( gSamplFreq*gDriftLen/gDriftVel/timeSlice );
	
	cout << "Debug --> Using " << nSlices << " slices to determine the hits." << endl;
	cout << "Debug --> Time slice: " << timeSlice << " samples = " << timeSlice/gSamplFreq << " usecs = " << gDriftVel*timeSlice/gSamplFreq << " mm." << endl;
	
	
	vector<Int_t> TimeEdges(nSlices+1);
	for(Int_t iSl=0; iSl<=nSlices; iSl++)
	{
		TimeEdges.at(iSl) = iSl*timeSlice;
	}
	
	
	Int_t nChs = t2w->fT1wr->ColHist->GetNbinsY();
	
	vector<vector<Double_t> >gColWfsVect(nChs);
	for(Int_t iCh=0; iCh<nChs; iCh++)
	{
		TH1D* hChWf = t2w->fT1wr->ColHist->ProjectionX("hChWf",iCh+1,iCh+1);
		Int_t nSamps = hChWf->GetNbinsX();
		vector<Double_t> wf(nSamps);
		for(Int_t iSamp=0; iSamp<nSamps; iSamp++)
		{
			wf.at(iSamp) = hChWf->GetBinContent(iSamp+1);
		}
		gColWfsVect.at(iCh) = wf;
	}
	
	
	
	vector<RSTPC_Hit*> fHitsVec;
	
	for(Int_t iSl=0; iSl<nSlices; iSl++)
	{
		Int_t ledge = TimeEdges.at(iSl);
		Int_t redge = TimeEdges.at(iSl+1);
		Int_t nsamps = redge-ledge+1;
		
		cout << "Debug --> Slice " << iSl << ": ledge=" << ledge << "\tredge=" << redge << endl;
		
		
		map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >::iterator mapIt;
		for(mapIt=fPulseMap->begin(); mapIt!=fPulseMap->end(); mapIt++)
		{
			ColPulse = mapIt->first;
			vector<RSTPC_Pulse*> *IndPulsesVec = mapIt->second;
			
			
			if( ( (ColPulse->fLedge>redge) && (ColPulse->fRedge<ledge) ) )
			{//There is no overlapping with the time slice
				continue;
			}
			
			
			
			//Here there is overlapping condition with the time slice check which of the induction pulses also overlap
			vector<RSTPC_Pulse*>::iterator vecIt;
			for(vecIt=IndPulsesVec->begin(); vecIt!=IndPulsesVec->end(); vecIt++)
			{
				RSTPC_Pulse* IndPulse = (*vecIt);
				
				if( ( (IndPulse->fLedge>redge) && (IndPulse->fRedge<ledge) ) )
				{//There is no overlapping with the time slice
					continue;
				}
				
				//Here there are both collection and induction pulses overlapping the time slice ==> Can make a hit
				RSTPC_Hit *hit = new RSTPC_Hit(ColPulse, IndPulse);
				hit->fLedge = (Double_t)ledge;
				hit->fRedge = (Double_t)redge;
				hit->SetCentreTime( (hit->fLedge+hit->fRedge)/2 );
				
				
				//Find the meaen (weighted) quantities
				Double_t meantime=0, meanheight=0, sum2=0;
				for(Int_t kSamp=ledge; kSamp<=redge; kSamp++)
				{//Using the square as a measure of the pulse chunk on the specific interval
					sum2 += pow((gColWfsVect.at(ColPulse->fWireNum)).at(kSamp) ,2);
					meantime += kSamp* pow((gColWfsVect.at(ColPulse->fWireNum)).at(kSamp) ,2); 
					meanheight += (gColWfsVect.at(ColPulse->fWireNum)).at(kSamp)/nsamps;
				}
				
				meantime = meantime/sum2;
				
				hit->fMeanTime = meantime;
				hit->fMeanHeight = meanheight;
				
				fHitsVec.push_back(hit);
				
				if(hit->fCentreTime <= 0)
				{
					cout << "WARNING -->HitsFinder(...): Non positive centre time!" << endl;
					cout << "           Hit ID: " << hit->fHitID << "->\tledge=" << hit->fLedge << "\tredge=" << hit->fRedge << "\tcentre time="<< hit->fCentreTime <<"\tmeantime=" << meantime << "\tmeanheight=" << meanheight << endl;
				}
				else if(meantime<=0)
				{
					cout << "WARNING ---> RSTPC_RunProcessor::HitsFinder(...): Non positive mean time!" << endl;
					cout << "           Hit ID: " << hit->fHitID << "->\tledge=" << hit->fLedge << "\tredge=" << hit->fRedge << "\tcentre time="<< hit->fCentreTime <<"\tmeantime=" << meantime << "\tmeanheight=" << meanheight << endl;
				}
			} //End of cycling over the induction pulses corresponding to a collection pulse
		} //End of cycling over the pulses map
	} //End of cycling over time slices
	
	
	//Avoid memory leaks
	if(gROOT->FindObject("hChWf")) delete gROOT->FindObject("hChWf");
	
	
	if(gHitsArr)
	{
		delete gHitsArr;
		gHitsArr = NULL;
	}
	
	gHitsArr = new TClonesArray("RSTPC_Hit", fHitsVec.size());
	
	TClonesArray &arr = *gHitsArr;
	
	for(Int_t iHit=0; iHit<fHitsVec.size(); iHit++)
	{
		new(arr[iHit]) RSTPC_Hit( *fHitsVec.at(iHit) );
	}
	
	return fHitsVec.size();
}
*/


void HitsFinderNoSlices(Int_t event)
{
	if(!t2w)
	{
		t2w = new RSTPC_T2wrapper("/home/francescop/data/ResistiveShell/merged/test/RSTPC_Run000002032_Merged.root", true);
	}
	
	if( !t2w->IsInit() ) return -1;
	
	RSTPC_RunProcessor::SetPitchSize( 52.5/31 ); //Measured
	RSTPC_RunProcessor::SetPeakingTime(1.0);
	RSTPC_RunProcessor::SetDriftVel(2.06);
	
	t2w->GetEntry(event);
	
	
	vector<RSTPC_Pulse*> ColPulsesVec(t2w->ColPulses->GetSize());
	vector<RSTPC_Pulse*> IndPulsesVec(t2w->IndPulses->GetSize());
	vector<RSTPC_Hit*> HitsVec;
	
	//This two guys are needed to make a fast search later on
	map<UInt_t, RSTPC_Pulse*> ColPulsesMap;
	map<UInt_t, RSTPC_Pulse*> IndPulsesMap;
	
	for(Int_t iPulse=0; iPulse<t2w->ColPulses->GetSize(); iPulse++)
	{
		RSTPC_Pulse* pulse = (RSTPC_Pulse*)t2w->ColPulses->At(iPulse);
		ColPulsesVec.at(iPulse) = pulse;
		ColPulsesMap[pulse->fPulseID] = pulse;
	}
	
	for(Int_t iPulse=0; iPulse<t2w->IndPulses->GetSize(); iPulse++)
	{
		RSTPC_Pulse* pulse = (RSTPC_Pulse*)t2w->IndPulses->At(iPulse);
		IndPulsesVec.at(iPulse) = pulse;
		IndPulsesMap[pulse->fPulseID] = pulse;
	}
	
	
	{
		vector<RSTPC_Pulse*>::iterator cVecIt;
		for(cVecIt=ColPulsesVec.begin(); cVecIt!=ColPulsesVec.end(); cVecIt++ )
		{
			RSTPC_Pulse* Pulse = (*cVecIt);
			if( !Pulse->fColCoinIDs )
			{
				Pulse->fColCoinIDs = new set<UInt_t>;
			}
			else
			{
				Pulse->fColCoinIDs->clear();
			}
			Pulse->fColCoinNum = 0;
			
			if( !Pulse->fIndCoinIDs )
			{
				Pulse->fIndCoinIDs = new set<UInt_t>;
			}
			else
			{
				Pulse->fIndCoinIDs->clear();
			}
			Pulse->fIndCoinNum = 0;
		}
		
		vector<RSTPC_Pulse*>::iterator iVecIt;
		for(iVecIt=IndPulsesVec.begin(); iVecIt!=IndPulsesVec.end(); iVecIt++ )
		{
			RSTPC_Pulse* Pulse = (*iVecIt);
			if( !Pulse->fColCoinIDs )
			{
				Pulse->fColCoinIDs = new set<UInt_t>;
			}
			else
			{
				Pulse->fColCoinIDs->clear();
			}
			Pulse->fColCoinNum = 0;
			
			if( !Pulse->fIndCoinIDs )
			{
				Pulse->fIndCoinIDs = new set<UInt_t>;
			}
			else
			{
				Pulse->fIndCoinIDs->clear();
			}
			Pulse->fIndCoinNum = 0;
		}
	}
	
	
	vector<RSTPC_Pulse*>::iterator cVecIt, iVecIt;
	
	//Collection pulses with amplitudes smaller than 200 ADC counts are considered noise
	set<UInt_t> ColPulsesNoiseIDs;
	
	for(cVecIt=ColPulsesVec.begin(); cVecIt!=ColPulsesVec.end(); cVecIt++ )
	{
		if( (*cVecIt)->fMax < 300 ) ColPulsesNoiseIDs.insert( (*cVecIt)->fPulseID );
	}
	
	set<UInt_t> IndPulsesNoiseIDs;
	
	for(iVecIt=IndPulsesVec.begin(); iVecIt!=IndPulsesVec.end(); iVecIt++ )
	{
		if( (*iVecIt)->fMax < 150 ) IndPulsesNoiseIDs.insert( (*iVecIt)->fPulseID );
	}
	
	
	
	for(cVecIt=ColPulsesVec.begin(); cVecIt!=ColPulsesVec.end(); cVecIt++ )
	{//Cycling over the collection pulses
		RSTPC_Pulse *ColPulse = (*cVecIt);
		
		if( ColPulsesNoiseIDs.find(ColPulse->fPulseID)!=ColPulsesNoiseIDs.end() ) continue;
		
		//Determine how many collection pulses are in coincidence with this pulse (inefficient algorithm)
		vector<RSTPC_Pulse*>::iterator cVecIt2;
		for(cVecIt2=ColPulsesVec.begin(); cVecIt2!=ColPulsesVec.end(); cVecIt2++ )
		{
			if(cVecIt2!=cVecIt)
			{
				if( ColPulsesNoiseIDs.find((*cVecIt2)->fPulseID)==ColPulsesNoiseIDs.end() )
				{
					if( !( (ColPulse->fLedge>(*cVecIt2)->fRedge)||(ColPulse->fRedge<(*cVecIt2)->fLedge) ) )
					{//This is the overlapping condition
						ColPulse->fColCoinIDs->insert( (*cVecIt2)->fPulseID );
					}
				}
			}
		}
		ColPulse->fColCoinNum = ColPulse->fColCoinIDs->size();
		
		
		//Fill the vector of the induction pulses in coincidence with this collection pulse (hits finder)
		for(iVecIt=IndPulsesVec.begin(); iVecIt!=IndPulsesVec.end(); iVecIt++ )
		{
			RSTPC_Pulse *IndPulse = (*iVecIt);
			
			if( IndPulsesNoiseIDs.find(IndPulse->fPulseID)!=IndPulsesNoiseIDs.end() ) continue;
			
			if( !( (ColPulse->fLedge>IndPulse->fRedge)||(ColPulse->fRedge<IndPulse->fLedge) ) )
			{//This is the overlapping condition
				ColPulse->fIndCoinIDs->insert( IndPulse->fPulseID );
				IndPulse->fColCoinIDs->insert( ColPulse->fPulseID );
				
			}
		}//End of cycling over all the induction pulses (to find the hits)
		ColPulse->fIndCoinNum = ColPulse->fIndCoinIDs->size();
		
	}//End of cycle over the collection pulses
	
	
	//Cycle over the induction pulses only to determine the coincidences with other induction pulses
	for(iVecIt=IndPulsesVec.begin(); iVecIt!=IndPulsesVec.end(); iVecIt++ )
	{
		RSTPC_Pulse *IndPulse = (*iVecIt);
		IndPulse->fIndCoinIDs->clear();//Should not be necessary
		
		if( IndPulsesNoiseIDs.find(IndPulse->fPulseID)!=IndPulsesNoiseIDs.end() ) continue;
		
		vector<RSTPC_Pulse*>::iterator iVecIt2;
		for(iVecIt2=IndPulsesVec.begin(); iVecIt2!=IndPulsesVec.end(); iVecIt2++ )
		{
			if(iVecIt2!=iVecIt)
			{
				if( IndPulsesNoiseIDs.find((*iVecIt2)->fPulseID)==IndPulsesNoiseIDs.end() )
				{
					if( !( (IndPulse->fLedge>(*iVecIt2)->fRedge) || (IndPulse->fRedge<(*iVecIt2)->fLedge) ) )
					{//This is the overlapping condition
						IndPulse->fIndCoinIDs->insert((*iVecIt2)->fPulseID);
					}
				}
			}
		}
		
		IndPulse->fColCoinNum = IndPulse->fColCoinIDs->size();
		IndPulse->fIndCoinNum = IndPulse->fIndCoinIDs->size();
	}
	
	
	
	///////////////////////////////////////////////////////////////////
	// Finished with finding all the coincidences, now find the hits //
	///////////////////////////////////////////////////////////////////
	
	//Cycle over the collection pulses
	for(cVecIt=ColPulsesVec.begin(); cVecIt!=ColPulsesVec.end(); cVecIt++ )
	{
		RSTPC_Pulse *ColPulse = (*cVecIt);
		
		if(ColPulse->fIndCoinIDs->size()==0) continue;
		
		if(ColPulse->fIndCoinIDs->size()==1)
		{//I consider the coincidence an hit only if the mean of the ind pulse is close enough to the maximum of the col pulse
			RSTPC_Pulse *IndPulse = IndPulsesMap[ (*ColPulse->fIndCoinIDs->begin()) ];
			if( abs( ((Double_t)ColPulse->fMaxPos)-IndPulse->fMeanTime )<=100 )
			{
				//Make a hit here
				RSTPC_Hit *hit = new RSTPC_Hit(ColPulse, IndPulse);
				
				//Copy the times amplitude and all pulse shape info from the collection pulse
				hit->fLedge = (Double_t)ColPulse->fLedge;
				hit->fRedge = (Double_t)ColPulse->fRedge;
				hit->SetCentreTime( 1.*(hit->fRedge+hit->fLedge)/2 );
				
				hit->fZ = -RSTPC_RunProcessor::GetDriftVel()*ColPulse->fMaxPos/RSTPC_RunProcessor::GetSamplingFreq();
				
				HitsVec.push_back(hit);
			}
		}
		else
		{//There several candidate induction pulses and therefore I need to select the one with best overlapping
			map<RSTPC_Pulse*, CrossCorrRes> CrossCorrMap; //Here there are the candidate induction pulses that can form an hit together with the collection pulse
			set<UInt_t>::iterator setIt;
			for(setIt=ColPulse->fIndCoinIDs->begin(); setIt!=ColPulse->fIndCoinIDs->end(); setIt++)
			{
				RSTPC_Pulse *IndPulse = IndPulsesMap[(*setIt)];
				
				if( abs( ((Double_t)ColPulse->fMaxPos)-IndPulse->fMeanTime )>100 )
				{//It is too off with respect to the collection pulse maximum
					continue;
				}
				CrossCorrRes cc_par = CalculatePulsesCrossCorrelation( ColPulse, t2w->fT1wr->ColHist, IndPulse, t2w->fT1wr->IndHist );
				if(!cc_par.bad) CrossCorrMap[IndPulse] = cc_par;
			}
			
			
			if(CrossCorrMap.size()>0)
			{
				//Iterate over the map to estabilish which induction pulse is the better for this
				map<RSTPC_Pulse*, CrossCorrRes>::iterator mapIt;
				RSTPC_Pulse* selIndPulse = NULL;
				for(mapIt=CrossCorrMap.begin(); mapIt!=CrossCorrMap.end(); mapIt++)
				{
					if(mapIt==CrossCorrMap.begin())
					{
						selIndPulse=mapIt->first;
					}
					else
					{
						CrossCorrRes cc_par = mapIt->second;
						/*
						if( abs(cc_par.taumax)<abs(CrossCorrMap[selIndPulse].taumax) )
						{
							selIndPulse=mapIt->first;
						}
						else
						{//The two induction pulses have the same tau of the cc maximum, cmpare the amplitudes
							if(cc_par.max>CrossCorrMap[selIndPulse].max) selIndPulse=mapIt->first;
						}
						*/
						if( (cc_par.max)>(CrossCorrMap[selIndPulse].max) )
						{
							selIndPulse=mapIt->first;
						}
					}
				}
				
				//With the selected pulse make an hit object
				RSTPC_Hit *hit = new RSTPC_Hit(ColPulse, selIndPulse);
				
				//Copy the times amplitude and all pulse shape info from the collection pulse
				hit->fLedge = (Double_t)ColPulse->fLedge;
				hit->fRedge = (Double_t)ColPulse->fRedge;
				hit->SetCentreTime( 1.*(hit->fRedge+hit->fLedge)/2 );
				
				hit->fZ = -RSTPC_RunProcessor::GetDriftVel()*ColPulse->fMaxPos/RSTPC_RunProcessor::GetSamplingFreq();
				
				HitsVec.push_back(hit);
			}
		}
	}//End cyclyng over the collection pulses for hit determination
	
	
	
	
	if(gHitsArr)
	{
		delete gHitsArr;
		gHitsArr = NULL;
	}
	
	gHitsArr = new TClonesArray("RSTPC_Hit", HitsVec.size());
	
	TClonesArray &arr = *gHitsArr;
	
	for(Int_t iHit=0; iHit<HitsVec.size(); iHit++)
	{
		new(arr[iHit]) RSTPC_Hit( *HitsVec.at(iHit) );
	}
	
	
	//////////////////////////////////////
	// Make plots of the hits positions //
	//////////////////////////////////////
	
	
	//The second set of bad events are those where one or more coincident collection pulses are above 300 ADC counts
	
	vector<RSTPC_Hit*>::iterator hitsIt;
	
	set<UInt_t> ColPulsesBadCoinIDs;
	
	for(hitsIt=HitsVec.begin(); hitsIt!=HitsVec.end(); hitsIt++ )
	{
		RSTPC_Hit* Hit = (*hitsIt);
		
		RSTPC_Pulse* CollPulse = ColPulsesMap[Hit->fColPulseID];
		
		vector<RSTPC_Hit*>::iterator hitsIt2;
		for(hitsIt2=HitsVec.begin(); hitsIt2!=HitsVec.end(); hitsIt2++ )
		{
			if(hitsIt==hitsIt2) continue;
			
			RSTPC_Hit* Hit2 = (*hitsIt2);
			RSTPC_Pulse* CollPulse2 = ColPulsesMap[Hit2->fColPulseID];
			
			
			if( abs(CollPulse->fMaxPos-CollPulse2->fMaxPos)<50 )
			{
				ColPulsesBadCoinIDs.insert(Hit->fColPulseID);
				ColPulsesBadCoinIDs.insert(Hit2->fColPulseID);
			}
			
		}
	}
	
	
	vector<RSTPC_Hit*> GoodHitsVec;
	for(hitsIt=HitsVec.begin(); hitsIt!=HitsVec.end(); hitsIt++)
	{
		RSTPC_Hit* hit = (*hitsIt);
		UInt_t ColPulseID = hit->fColPulseID;
		UInt_t IndPulseID = hit->fIndPulseID;
		
		if( (ColPulsesBadCoinIDs.find(ColPulseID)==ColPulsesBadCoinIDs.end()) )
		//if( (ColPulsesNoiseIDs.find(ColPulseID)==ColPulsesNoiseIDs.end()) && (IndPulsesNoiseIDs.find(IndPulseID)==IndPulsesNoiseIDs.end()) )
		//if( (ColPulsesNoiseIDs.find(ColPulseID)==ColPulsesNoiseIDs.end()) )
		//if( (IndPulsesNoiseIDs.find(IndPulseID)==IndPulsesNoiseIDs.end()) )
		//if(true)
		{
			GoodHitsVec.push_back(hit);
		}
	}
	
	cout << "\nHitsFinderNoSlices summary for event " << event << ":" << endl;
	cout << "\tTotal collection pulses: " << ColPulsesVec.size() << endl;
	cout << "\tNoisy pulses: " << ColPulsesNoiseIDs.size() << " (" << 100.*ColPulsesNoiseIDs.size()/ColPulsesVec.size() << "%)" << endl;
	cout << "\tMultiple coincidence pulses: " << ColPulsesBadCoinIDs.size() << " (" << 100.*ColPulsesBadCoinIDs.size()/ColPulsesVec.size() << "%)" << endl;
	cout << "\tGood hits: " << GoodHitsVec.size() << " out of " << HitsVec.size() << " (" << 100.*GoodHitsVec.size()/HitsVec.size() << "%)" << endl << endl;
	
	if(!GoodHitsVec.size()) return;
	
	
	TCanvas *c1 = (TCanvas*)gROOT->FindObject( "canvHitPatterns" );
	if(!c1)
	{
		c1 = new TCanvas("canvHitPatterns", "Hit Patterns", 800, 1000);
		c1->Divide(2,1);
	}
	
	TGraph *hitsXZ = new TGraph();
	hitsXZ->SetTitle(";X [mm];Z [mm]");
	hitsXZ->SetMarkerStyle(7);
	
	TGraphErrors *hitsXZwidths = new TGraphErrors();
	hitsXZwidths->SetMarkerStyle(1);
	hitsXZwidths->SetMarkerColor(kRed);
	hitsXZwidths->SetLineColor(kRed);
	
	TGraph *hitsYZ = new TGraph();
	hitsYZ->SetTitle(";Y [mm];Z [mm]");
	hitsYZ->SetMarkerStyle(7);
	
	TGraphErrors *hitsYZwidths = new TGraphErrors();
	hitsYZwidths->SetMarkerStyle(1);
	hitsYZwidths->SetMarkerColor(kRed);
	hitsYZwidths->SetLineColor(kRed);
	
	//cout << "\nCollection pulses of good hits:" << endl;
	//cout << "\nInduction pulses of good hits:" << endl;
	for(Int_t iPoint=0; iPoint<GoodHitsVec.size(); iPoint++)
	{
		RSTPC_Hit* hit = GoodHitsVec.at(iPoint);
		hitsXZ->SetPoint(iPoint, hit->fX, hit->fZ);
		hitsYZ->SetPoint(iPoint, hit->fY, hit->fZ);
		
		
		
		RSTPC_Pulse* iPulse = IndPulsesMap[hit->fIndPulseID];
		
		//double centreZ = -RSTPC_RunProcessor::GetDriftVel()*(iPulse->fRedge+iPulse->fLedge)/2/RSTPC_RunProcessor::GetSamplingFreq();
		//double width = RSTPC_RunProcessor::GetDriftVel()*(iPulse->fRedge-iPulse->fLedge)/RSTPC_RunProcessor::GetSamplingFreq();
		//cout << "Pulse "<< hit->fIndPulseID<< ": " << "Centre Z = " << centreZ << ", width = " << width << endl;
		
		double centreZ = -RSTPC_RunProcessor::GetDriftVel()*(iPulse->fMeanTime)/RSTPC_RunProcessor::GetSamplingFreq();
		double width = RSTPC_RunProcessor::GetDriftVel()*sqrt(iPulse->fSigma)/RSTPC_RunProcessor::GetSamplingFreq();
		//cout << "Pulse "<< hit->fIndPulseID<< ": " << "Mean Z = " << centreZ << ", width = " << width << endl;
		
		hitsXZwidths->SetPoint(iPoint, hit->fX, centreZ);
		hitsXZwidths->SetPointError(iPoint, 0, width);
		
		
		RSTPC_Pulse* cPulse = ColPulsesMap[hit->fColPulseID];
		
		//centreZ = -RSTPC_RunProcessor::GetDriftVel()*(cPulse->fRedge+cPulse->fLedge)/2/RSTPC_RunProcessor::GetSamplingFreq();
		//width = RSTPC_RunProcessor::GetDriftVel()*(cPulse->fRedge-cPulse->fLedge)/RSTPC_RunProcessor::GetSamplingFreq();
		
		centreZ = -RSTPC_RunProcessor::GetDriftVel()*(cPulse->fMeanTime)/RSTPC_RunProcessor::GetSamplingFreq();
		width = RSTPC_RunProcessor::GetDriftVel()*sqrt(cPulse->fSigma)/RSTPC_RunProcessor::GetSamplingFreq();
		
		//cout << "Pulse "<< hit->fColPulseID<< ": " << "Centre Z = " << centreZ << ", width = " << width << endl;
		
		hitsYZwidths->SetPoint(iPoint, hit->fY, centreZ);
		hitsYZwidths->SetPointError(iPoint, 0, width);
	}
	
	
	cout << "\nGood hits:" << endl;
	for(Int_t iPoint=0; iPoint<GoodHitsVec.size(); iPoint++)
	{
		RSTPC_Hit* hit = GoodHitsVec.at(iPoint);
		cout << "Hit "<< hit->fHitID << ": " << "Z = " << hit->fZ << ", X = " << hit->fX << ", Y = " << hit->fY << endl;
	}
	
	
	
	
	TH1F *frame1 = c1->cd(1)->DrawFrame(-5, -155, 57.5, 5);
	TMultiGraph *mg_XZhits = new TMultiGraph;
	mg_XZhits->Add(hitsXZ, "P");
	//mg_XZhits->Add(hitsXZwidths, "PE");
	//hitsXZ->Draw("AP");
	//hitsXZwidths->Draw("AP");
	mg_XZhits->Draw("P");
	
	TH1F *frame2 = c1->cd(2)->DrawFrame(-5, -155, 57.5, 5);
	TMultiGraph *mg_YZhits = new TMultiGraph;
	mg_YZhits->Add(hitsYZ, "P");
	//mg_YZhits->Add(hitsYZwidths, "PE");
	//hitsYZ->Draw("AP");
	//hitsYZwidths->Draw("AP");
	mg_YZhits->Draw("P");
}

/*
void PlotHitTimes()
{
	if(!t2w)
	{
		t2w = new RSTPC_T2wrapper("/home/francescop/data/ResistiveShell/merged/test/RSTPC_Run000002032_Merged.root", true);
	}
	
	if( !t2w->IsInit() ) return;
	
	
	
	if (t2w->fChain == 0) return;
	
	Long64_t nEvs = t2w->fChain->GetEntriesFast();
	
	
	TH1D *hHitCentreTimes = (TH1D*)gROOT->FindObject("hHitCentreTimes");
	if(!hHitCentreTimes)
	{
		hHitCentreTimes = new TH1D("hHitCentreTimes",";Hit centre time (samples);Counts", 4096, -0.5, 4096-0.5);
	}
	else
	{
		hHitCentreTimes->Reset();
	}
	
	TH1D *hHitMeanTimes = (TH1D*)gROOT->FindObject("hHitMeanTimes");
	if(!hHitMeanTimes)
	{
		hHitMeanTimes = new TH1D("hHitMeanTimes",";Hit mean time (samples);Counts", 4096, -0.5, 4096-0.5);
	}
	else
	{
		hHitMeanTimes->Reset();
	}
	
	
	
	RSTPC_Hit *hit;
	
	for(Int_t iEv=0; iEv<nEvs; iEv++)
	{
		t2w->GetEntry(iEv);
		
		TIter HitsIt(t2w->Hits);
		
		while( (hit = (RSTPC_Hit*)HitsIt.Next()) )
		{
			hHitCentreTimes->Fill( hit->fCentreTime );
			hHitMeanTimes->Fill( hit->fMeanTime );
		}
	}
	
	TCanvas *c1 = (TCanvas*)gROOT->FindObject("canvHitTimes");
	if(c1) delete c1;
	c1 = new TCanvas("canvPulsesLen", "Hits times", 1200, 800);
	c1->Divide(1,2);
	
	c1->cd(1);
	hHitCentreTimes->Draw("hist");
	
	c1->cd(2);
	hHitMeanTimes->Draw("hist");
}
*/


CrossCorrRes CalculatePulsesCrossCorrelation( RSTPC_Pulse* ColPulse, TH2D* ColHist, RSTPC_Pulse* IndPulse, TH2D* IndHist )
{
	CrossCorrRes result;
	result.bad = false;
	
	if(!(ColPulse && ColHist && IndPulse && IndHist))
	{
		result.bad=true;
		return result;
	}
	
	vector<Double_t> ColPulseWf(ColHist->GetNbinsX(), 0.);
	vector<Double_t> IndPulseWf(IndHist->GetNbinsX(), 0.);
	
	Int_t nSamps = ColPulseWf.size();
	if(IndPulseWf.size()!=ColPulseWf.size())
	{
		result.bad=true;
		return result;
	}
	
	for(Int_t iSamp=ColPulse->fLedge; iSamp<=ColPulse->fRedge; iSamp++)
	{
		ColPulseWf.at(iSamp) = ColHist->GetBinContent(iSamp+1, ColPulse->fWireNum+1);
	}
	
	for(Int_t iSamp=IndPulse->fLedge; iSamp<=IndPulse->fRedge; iSamp++)
	{
		IndPulseWf.at(iSamp) = IndHist->GetBinContent(iSamp+1, IndPulse->fWireNum+1);
	}
	
	
	//Minimum and maximum CC times (HARD CODED HERE)
	Double_t DeltaT = 1./RSTPC_RunProcessor::GetSamplingFreq();
	Int_t TauMaxShift = RSTPC_RunProcessor::GetCrosCorrsMaxAbsDelaySamps(); //In samples
	Int_t nTimes = 2*TauMaxShift + 1;
	
	for(Int_t iTime=0; iTime<nTimes; iTime++)
	{
		Int_t SampShift = iTime-TauMaxShift;
		Double_t tau = SampShift*DeltaT;
		Double_t croscorr = 0;
		
		for(Int_t iSamp=0; iSamp<ColPulseWf.size(); iSamp++)
		{
			if( (iSamp+SampShift>=0) && (iSamp+SampShift<nSamps) )
			{
				croscorr += ColPulseWf.at(iSamp)*IndPulseWf.at(iSamp+SampShift);
			}
		}
		
		if(iTime==0)
		{
			result.taumax = tau;
			result.max = croscorr;
		}
		else
		{
			if( result.max < croscorr )
			{
				result.taumax = tau;
				result.max = croscorr;
			}
		}
	}
	
	return result;
}