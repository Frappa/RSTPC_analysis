#ifndef RSTPC_RUNPROCESSOR_CC
#define RSTPC_RUNPROCESSOR_CC

#include "RSTPC_RunProcessor.hh"
#include "RSTPC_Hits.hh"

#include "TClassTable.h"

#include <vector>
#include <map>
#include <set>
#include <cmath>

using namespace std;





//#ifndef RSTPC_RUNPROC_STATICS
Bool_t RSTPC_RunProcessor::fgDebug = false;
Int_t RSTPC_RunProcessor::fgNevPrint = 100;
Double_t RSTPC_RunProcessor::fgSigmaThr = 3.0;
Double_t RSTPC_RunProcessor::fgPeakingTime = 0.0;
Double_t RSTPC_RunProcessor::fgSamplingFreq = 1./50e-3;
Double_t RSTPC_RunProcessor::fgPitchSize = 1.875;
Double_t RSTPC_RunProcessor::fgDriftLenght = 150;
Double_t RSTPC_RunProcessor::fgDriftVel = 0.0;
//#endif


void RSTPC_RunProcessor::SetDebug(Bool_t flag, Int_t nEvs)
{
	fgDebug = flag;
	fgNevPrint = nEvs;
}

void RSTPC_RunProcessor::SetSigmaThr(Double_t _val)
{
	fgSigmaThr = _val;
}
Double_t RSTPC_RunProcessor::GetSigmaThr()
{
	return fgSigmaThr;
}

void RSTPC_RunProcessor::SetPeakingTime(Double_t _val)
{
	fgPeakingTime = _val;
}
Double_t RSTPC_RunProcessor::GetPeakingTime()
{
	return fgPeakingTime;
}

void RSTPC_RunProcessor::SetSamplingFreq(Double_t _val)
{
	fgSamplingFreq = _val;
}
Double_t RSTPC_RunProcessor::GetSamplingFreq()
{
	return fgSamplingFreq;
}

void RSTPC_RunProcessor::SetPitchSize(Double_t _val)
{
	fgPitchSize = _val;
}
Double_t RSTPC_RunProcessor::GetPitchSize()
{
	return fgPitchSize;
}

void RSTPC_RunProcessor::SetDriftLenght(Double_t _val)
{
	fgDriftLenght = _val;
}
Double_t RSTPC_RunProcessor::GetDriftLenght()
{
	return fgDriftLenght;
}

void RSTPC_RunProcessor::SetDriftVel(Double_t _val)
{
	fgDriftVel = _val;
}
Double_t RSTPC_RunProcessor::GetDriftVel()
{
	return fgDriftVel;
}


RSTPC_RunProcessor::RSTPC_RunProcessor() : gColPulses(0), gIndPulses(0)
{
	fRunNumber = -1;
	
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fDataDir = RSTPC_Options::GetInstance()->GetDataDir();
		if(fDataDir.at(fDataDir.size()-1)!='/') fDataDir += "/";
		
		if( access(fDataDir.c_str(), X_OK|W_OK|R_OK) )
		{
			exit(-1);
		}
	}
	else
	{
		if(!RSTPC_Options::GetInstance()->IsOutDirSet())
		{
			exit(-2);
		}
	}
	
	
	if( RSTPC_Options::GetInstance()->IsOutDirSet() )
	{//Here I already have the output directory ==> if I want only T2/T3 processsing I do not need the data dir
		fOutDir = RSTPC_Options::GetInstance()->GetOutDir();
		if( access(fOutDir.c_str(), X_OK|W_OK|R_OK) )
		{
			gSystem->mkdir(fOutDir.c_str(),true);
		}
		if( access(fOutDir.c_str(), X_OK|W_OK|R_OK) )
		{
			exit(-3);
		}
	}
	else
	{//In this case the data dir name must be given and the out dir will be created inside it
		if( !RSTPC_Options::GetInstance()->IsDataDirSet() )
		{
			exit(-2);
		}
		
		fDataDir = RSTPC_Options::GetInstance()->GetDataDir();
		
		if(!fDataDir.size())
		{
			exit(-2);
		}
		
		if(fDataDir.at(fDataDir.size()-1)!='/') fDataDir += "/";
		
		fOutDir = fDataDir + string("merged/");
		if( access(fOutDir.c_str(), X_OK|W_OK|R_OK) )
		{
			gSystem->mkdir(fOutDir.c_str(),true);
		}
		if( access(fOutDir.c_str(), X_OK|W_OK|R_OK) )
		{
			exit(-4);
		}
	}
	
	fTpcMan = new RSTPC_Analyser;
	fFebMan = new MppcTreeWrapper;
	fTimeDiffVec = new vector<Double_t>;
	fTpcTimes = new vector<Double_t>;
	fFebTimes = new vector<Double_t>;
	fFebEvs = new vector<Int_t>;
	
	fOutFile = NULL;
	fOutT1 = NULL;
	fT1wr = NULL;
	fOutT2 = NULL;
	
	fRunOpenFlag = false;
	
	if( RSTPC_Options::GetInstance()->IsRunNumberSet() )
	{
		fRunOpenFlag = InitT1proc( RSTPC_Options::GetInstance()->GetRun() );
	}
	
	fProcT1 = false;
	fProcT2 = false;
	
	gEventData = NULL;//This is allocated in the processors
	
	gColWfsVect = NULL;
	gIndWfsVect = NULL;
	
	gColPulses = new vector<RSTPC_Pulse*>;
	gIndPulses = new vector<RSTPC_Pulse*>;
	gHits = new vector<RSTPC_Hit*>;
}


RSTPC_RunProcessor::RSTPC_RunProcessor(Int_t RunNumber) : gColPulses(0), gIndPulses(0)
{
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fDataDir = RSTPC_Options::GetInstance()->GetDataDir();
		if(fDataDir.at(fDataDir.size()-1)!='/') fDataDir += "/";
		
		if( access(fDataDir.c_str(), X_OK|W_OK|R_OK) )
		{
			exit(-1);
		}
	}
	else
	{
		if(!RSTPC_Options::GetInstance()->IsOutDirSet())
		{
			exit(-2);
		}
	}
	
	if( RSTPC_Options::GetInstance()->IsOutDirSet() )
	{
		fOutDir = RSTPC_Options::GetInstance()->GetOutDir();
		if( access(fOutDir.c_str(), X_OK|W_OK|R_OK) )
		{
			gSystem->mkdir(fOutDir.c_str(),true);
		}
		if( access(fOutDir.c_str(), X_OK|W_OK|R_OK) )
		{
			exit(-3);
		}
	}
	else
	{
		fOutDir = fDataDir + string("merged/");
		if( access(fOutDir.c_str(), X_OK|W_OK|R_OK) )
		{
			gSystem->mkdir(fOutDir.c_str(),true);
		}
		if( access(fOutDir.c_str(), X_OK|W_OK|R_OK) )
		{
			exit(-4);
		}
	}
	
	fTpcMan = new RSTPC_Analyser;
	fFebMan = new MppcTreeWrapper;
	fTimeDiffVec = new vector<Double_t>;
	fTpcTimes = new vector<Double_t>;
	fFebTimes = new vector<Double_t>;
	fFebEvs = new vector<Int_t>;
	
	fOutFile = NULL;
	fOutT1 = NULL;
	fT1wr = NULL;
	fOutT2 = NULL;
	
	fRunNumber = RunNumber;
	
	//fRunOpenFlag = InitT1proc(RunNumber);
	
	
	fProcT1 = false;
	fProcT2 = false;
	
	gEventData = NULL;//This is allocated in the processors
	
	gColWfsVect = NULL;
	gIndWfsVect = NULL;
	
	gColPulses = new vector<RSTPC_Pulse*>;
	gIndPulses = new vector<RSTPC_Pulse*>;
	gHits = new vector<RSTPC_Hit*>;
}


RSTPC_RunProcessor::~RSTPC_RunProcessor()
{
	if(fOutFile) delete fOutFile;
	
	if(fTpcMan) delete fTpcMan;
	if(fFebMan) delete fFebMan;
	
	if(fTimeDiffVec) delete fTimeDiffVec;
	if(fTpcTimes) delete fTpcTimes;
	if(fFebTimes) delete fFebTimes;
	if(fFebEvs) delete fFebEvs;
	
	if(!gEventData){
		delete gEventData;
		gEventData = NULL;
	}
	
	if(gColPulses) delete gColPulses;
	if(gIndPulses) delete gIndPulses;
	if(gHits) delete gHits;
}



Bool_t RSTPC_RunProcessor::InitT1proc(Int_t RunNumber)
{
	if(!RSTPC_Options::GetInstance())
	{
		return false;
	}
	
	if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		string datadir = RSTPC_Options::GetInstance()->GetDataDir();
		if(!datadir.size()) return false;
		
		fDataDir = datadir;
	}
	else
	{
		return false;
	}
	
	if( RunNumber<0 )
	{
		if(!RSTPC_Options::GetInstance()->IsRunNumberSet()) return false;
		
		fRunNumber = RSTPC_Options::GetInstance()->GetRun();
	}
	else
	{
		fRunNumber = RunNumber;
	}
	
	
	if(fDataDir.at(fDataDir.size()-1)!='/') fDataDir += "/";
	
	fTpcMan->fDataDir = fDataDir;
	fTpcMan->OpenRun(fRunNumber);
	if( !(fTpcMan->IsInit() && (fTpcMan->fTrigTree)) ) return false;
	
	fFebMan->fDataDir = fDataDir;
	stringstream feb_filename; feb_filename.str("");
	feb_filename << fDataDir + "July_run_" << fRunNumber << "_Cosmics.root";
	if( !fFebMan->Open(feb_filename.str()) ) return false;
	
	
	return true;
}


void RSTPC_RunProcessor::DescribeT1()
{
	fProcT1 = false;
	
	if(gEventData) delete gEventData;
	gEventData = new EventData;
	
	if( !(fTpcMan->IsInit() && fFebMan->IsInit() && (fTpcMan->fTrigTree)) ) return;
	
	/*
	string outdir;
	//string outdir = Options::Instance()->GetOutDir();
	if(outdir.at(outdir.size()-1)!='/') outdir += "/";
	
	if( access(outdir.c_str(), X_OK|W_OK|R_OK) )
	{
		gSystem->mkdir( outdir.c_str() );
	}
	
	if( access(outdir.c_str(), X_OK|W_OK|R_OK) ) return;
	*/
	
	if(fOutDir.at(fOutDir.size()-1)!='/') fOutDir += "/";
	
	stringstream outfilename; outfilename.str("");
	outfilename << fOutDir << "RSTPC_Run" << setfill('0') << setw(9) << fRunNumber << "_Merged.root";
	
	fOutFile = TFile::Open( outfilename.str().c_str(), "recreate" );
	
	fOutT1 = new TTree("T1","Tear 1 data of the merged tree for the Resistive Shell TPC data");
	
	fOutT1->Branch("TpcEvent", &gEventData->TpcEv, "TpcEvent/I");
	fOutT1->Branch("TpcTime", &gEventData->TpcTime, "TpcTime/D");
	fOutT1->Branch("RmsColWires", gEventData->RmsColWires, "RmsColWires[32]/D");
	fOutT1->Branch("RmsIndWires", gEventData->RmsIndWires, "RmsIndWires[32]/D");
	//fOutT1->Branch("ColHist", "TH2D", &gEventData->hCol);
	//fOutT1->Branch("IndHist", "TH2D", &gEventData->hInd);
	fOutT1->Branch("FebEvent", &gEventData->FebEv, "FebEvent/I");
	fOutT1->Branch("FebTime", &gEventData->FebTime, "FebTime/D");
	fOutT1->Branch("FebTopAmp", gEventData->FebTopAmp, "FebTopAmp[3]/s");
	fOutT1->Branch("FebTopTotAmp", &gEventData->FebTopTotAmp, "FebTopTotAmp/D");
	fOutT1->Branch("FebBotAmp", gEventData->FebBotAmp, "FebBotAmp[3]/s");
	fOutT1->Branch("FebBotTotAmp", &gEventData->FebBotTotAmp, "FebBotTotAmp/D");
	fOutT1->Branch("FebTotAmp", &gEventData->FebTotAmp, "FebTotAmp/D");
	
	fProcT1 = true;
}


void RSTPC_RunProcessor::T1process()
{
	if(!fProcT1)
	{
		cout << "\nERROR --> RSTPC_RunMerger::T1process(): \"fProcT1\" flag is false. Abort processing.\n" << endl;
		return;
	}
	
	Int_t nTpcEvs = fTpcMan->fTrigTree->GetEntries();
	Int_t nFebEvs = fFebMan->fChain->GetEntries();
	
	fTpcTimes->assign(nTpcEvs, -1.);
	fFebEvs->assign(nTpcEvs, -1);
	fFebTimes->assign(nFebEvs, -1.);
	
	set<Int_t> BadTpcEvs;
	stringstream ColHistName, IndHistName;
	
	TFile *TpcFile = fTpcMan->fInfile;
	
	//Before playing with the tree determine which FEB event corresponds to which TPC event
	for(Int_t iTrig=0; iTrig<nTpcEvs; iTrig++)
	{
		fTpcMan->fTrigTree->GetEntry(iTrig);
		fTpcTimes->at(iTrig) = ((Double_t)fTpcMan->fEventTime)/1e6;
		
		ColHistName.str(""); ColHistName << "Col_" << iTrig;
		IndHistName.str(""); IndHistName << "Ind_" << iTrig;
		
		TH2 *hCol = (TH2*)TpcFile->Get(ColHistName.str().c_str());
		TH2 *hInd = (TH2*)TpcFile->Get(IndHistName.str().c_str());
		
		if( !( hCol && hInd ) )
		{
			BadTpcEvs.insert(iTrig);
			cout << "\nWARNING --> The TPC event " << iTrig << " is marked as bad event. It will not be processed.\n" << endl;
		}
	}
	
	cout << "Debug --> RSTPC_RunMerger::T1process: " << BadTpcEvs.size() << " number of TPC bad events out of " << nTpcEvs <<  " total." << endl;
	
	set<Int_t> FebFreeEvs;
	for(int iEv=0; iEv<nFebEvs; iEv++)
	{
		fFebMan->fChain->GetEntry(iEv);
		fFebTimes->at(iEv) = fFebMan->GetEventTime();
		FebFreeEvs.insert(iEv);
	}
	
	
	
	//This object contains the correspondences between Tpc events and Feb events
	map<Int_t, Int_t> CoinMap; //Key=TpcEv, Val=FebEv
	
	for(Int_t iTrig=0; iTrig<nTpcEvs; iTrig++)
	{
		if( BadTpcEvs.find(iTrig)!=BadTpcEvs.end() ) continue;
		
		Double_t minabsdiff;
		
		set<Int_t>::iterator iEv, selEv;
		for(iEv=FebFreeEvs.begin(); iEv!=FebFreeEvs.end(); ++iEv)
		{
			Double_t diff = ( fTpcTimes->at(iTrig)-fFebTimes->at(*iEv) );
			if( iEv==FebFreeEvs.begin() )
			{
				minabsdiff = TMath::Abs(diff);
				selEv = iEv;
			}
			else if( TMath::Abs(diff)<minabsdiff )
			{
				minabsdiff = TMath::Abs(diff);
				selEv = iEv;
			}
		}//End cycle over the FEB events
		
		CoinMap[iTrig] = *selEv;
		FebFreeEvs.erase(selEv);
	}
	
	cout << "Debug --> RSTPC_RunMerger::T1process: " << CoinMap.size() << " number of Tpc events to be processed." << endl;
	
	//Now I can fill the T1 tree simply looping through the "CoinMap" object
	map<Int_t, Int_t>::iterator evIt;
	for( evIt=CoinMap.begin(); evIt!=CoinMap.end(); ++evIt )
	{
		Int_t iTpcEv = evIt->first;
		Int_t iFebEv = evIt->second;
		
		fTpcMan->LoadEvent(iTpcEv);
		fTpcMan->ApplyWiresMaps(fTpcMan->hC0,fTpcMan->hI0);
		fTpcMan->fTrigTree->GetEntry(iTpcEv);
		fFebMan->fChain->GetEntry(iFebEv);
		
		//Process the TPC event for CM rejection
		fTpcMan->BaselineCorr(fTpcMan->hC0, fTpcMan->fColRMS);
		fTpcMan->BaselineCorr(fTpcMan->hI0, fTpcMan->fIndRMS);
		
		fTpcMan->CMrej(fTpcMan->hC0, fTpcMan->fColRMS, fTpcMan->fCMrejIter, true);
		fTpcMan->CMrej(fTpcMan->hI0, fTpcMan->fIndRMS, fTpcMan->fCMrejIter, true);
		
		
		gEventData->TpcEv = iTpcEv;
		gEventData->TpcTime = ((Double_t)fTpcMan->fEventTime)/1e6;
		//gEventData->hCol = fTpcMan->hC0;
		//gEventData->hInd = fTpcMan->hI0;
		
		for(Int_t iWire=0; iWire<32; iWire++)
		{
			gEventData->RmsColWires[iWire] = fTpcMan->fColRMS->at(iWire);
			gEventData->RmsIndWires[iWire] = fTpcMan->fIndRMS->at(iWire);
		}
		
		
		gEventData->FebEv = iFebEv;
		gEventData->FebTime = fFebMan->GetEventTime();
		
		//The actually used FEB channels connected to the APDs are are hardcoded here below
		gEventData->FebTopAmp[0] = fFebMan->chg[1];
		gEventData->FebTopAmp[1] = fFebMan->chg[2];
		gEventData->FebTopAmp[2] = fFebMan->chg[3];
		gEventData->FebTopTotAmp = gEventData->FebTopAmp[0] + gEventData->FebTopAmp[1] + gEventData->FebTopAmp[2];
		
		gEventData->FebBotAmp[0] = fFebMan->chg[4];
		gEventData->FebBotAmp[1] = fFebMan->chg[6];
		gEventData->FebBotAmp[2] = fFebMan->chg[7];
		gEventData->FebBotTotAmp = gEventData->FebBotAmp[0] + gEventData->FebBotAmp[1] + gEventData->FebBotAmp[2];
		
		gEventData->FebTotAmp = gEventData->FebTopTotAmp + gEventData->FebBotTotAmp;
		
		fOutT1->Fill();
		
		fOutFile -> WriteTObject(fTpcMan->hC0,0,"overwrite");
		fOutFile -> WriteTObject(fTpcMan->hI0,0,"overwrite");
		
	}// End cycling over the events
	
	
	fOutFile -> WriteTObject(fOutT1,0,"overwrite");
	fOutT1 -> ResetBranchAddresses();
	
	if(gEventData)
	{
		delete gEventData;
		gEventData=NULL;
	}
}



Bool_t RSTPC_RunProcessor::InitT2proc()
{
	if(fProcT2)
	{//It might be that the initialisation was already performed
		
		if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): The class might already be initialised for the T2 processing." << endl;
		
		if( fOutFile && fOutFile->IsWritable() )
		{//check that the file name is the correct one
			string fname = fOutFile->GetName();
			
			if( RSTPC_Options::GetInstance()->IsOutDirSet() && ((RSTPC_Options::GetInstance()->GetOutFile()).size()>0) )
			{//Here I already give the out file that however is also the the input file where the T1 tree is stored
				string gfname;
				gfname = RSTPC_Options::GetInstance()->GetOutDir();
				if(gfname.at(gfname.size()-1)!='/') gfname += "/";
				gfname += RSTPC_Options::GetInstance()->GetOutFile();
				
				if(fname == gfname)
				{
					if(fOutT2)
					{
						if(fOutFile!=fOutT2->GetCurrentFile())
						{
							delete fOutT2;
							fOutT2 = NULL;
						}
					}
					fOutFile->cd();
					if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): No need to make an initialisation for the T2 processing." << endl;
					return true;
				}
				else
				{
					if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): The output file name do not match. Deleting the \"fOutFile\" pointer." << endl;
					delete fOutFile;
					fOutFile = NULL;
				}
			}
		}
	}
	fProcT2 = false;
	
	if(fOutT2)
	{
		if(fOutT2->GetCurrentFile()) delete fOutT2->GetCurrentFile();
		fOutFile = NULL;
		fOutT2 = NULL;
	}
	
	
	if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): Start the initialisation for T2 processing...." << endl;
	stringstream outfilename;
	
	if(fRunOpenFlag)
	{//I am making the T1 processing also check if the T1 is present
		if( (!fOutT1) || (fOutT1->GetEntries()==0) || (!fOutT1->GetCurrentFile()->IsWritable()) )
		{//There is an inconsistency problem!!!
			return false;
		}
		
		if(!fT1wr)
		{
			fT1wr = new RSTPC_T1wrapper( fOutT1 );
		}
		else
		{//Check if the right file is open and if the class is properly initialised
			if( !(fT1wr->fInfile==fOutT1->GetCurrentFile()) )
			{
				delete fT1wr;
				fT1wr = new RSTPC_T1wrapper( fOutT1 );
			}
			else if( !fT1wr->IsInit() )
			{
				fT1wr->Init( (TTree*)fOutFile->Get("T1") );
			}
		}
		
		
		//At this point either the T1wrapper class is initialised or better to delete it and return a false flag
		if(!fT1wr->IsInit())
		{
			delete fT1wr;
			fT1wr = NULL;
			return false;
		}
		
		return true;
	}
	else
	{//Here I will not "open the run" (the original dta files) as all the info are already in the outfile where also T1 is
		
		if(RSTPC_Options::GetInstance()->IsOutFileMode())
		{//Here I already give the out file that however is also the the input file where the T1 tree is stored
			outfilename.str("");
			outfilename << RSTPC_Options::GetInstance()->GetOutDir();
			if(outfilename.str().at(outfilename.str().size()-1)!='/') outfilename << "/";
			outfilename << RSTPC_Options::GetInstance()->GetOutFile();
			
			
			if( access(outfilename.str().c_str(), R_OK|W_OK )!=0 )
			{
				cout << "ERROR --> RSTPC_RunProcessor::InitT2proc(): Cannot access in read and write mode to file: <" << outfilename.str() << ">" << endl;
				return false;
			}
			
			if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): Output file name (complete path): <" << outfilename.str() << ">" << endl;
		}
		else
		{//Here I try to find the right file starting fromm the run number and either the data directory or an output directory where the file with proper naming and containing the T1 should already be
			if(!RSTPC_Options::GetInstance()->IsRunNumberSet())
			{//I don't know which file to open
				return false;
			}
			fRunNumber = RSTPC_Options::GetInstance()->GetRun();
			
			if(!RSTPC_Options::GetInstance()->IsOutDirSet())
			{//I do not know the output directory check if the data dir exists and if there the file with the T1 tree inside
				if(!RSTPC_Options::GetInstance()->IsDataDirSet())
				{//Here I cannot do much
					return false;
				}
				
				//I can compose the name of the output file
				fOutDir = RSTPC_Options::GetInstance()->GetDataDir() + string("/merged/");
			}
			else
			{
				fOutDir = RSTPC_Options::GetInstance()->GetOutDir();
			}
			
			outfilename.str("");
			outfilename << fOutDir << "RSTPC_Run" << setfill('0') << setw(9) << fRunNumber << "_Merged.root";
			if( access(outfilename.str().c_str(), R_OK|W_OK) )
			{
				return false;
			}
		}
		//Here I determined that the output file exists and that it is writable
		
		if(fgDebug)
		{
			cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): Determined output directory: <" << fOutDir << ">." << endl;
			cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): Determined output file: <" << outfilename.str() << ">." << endl;
			cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): Checking the T1 wrapper...." << endl;
		}
		
		if(!fOutFile)
		{
			fOutFile = TFile::Open(outfilename.str().c_str(), "update");
			if(!fOutFile)
			{
				cout << "ERROR --> RSTPC_RunProcessor::InitT2proc(): Cannot access in update mode to file: <" << outfilename.str() << ">" << endl;
				return false;
			}
		}
		
		if(!fT1wr)
		{
			cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): Allocating the T1 wrapper object with the root file in update mode." << endl;
			
			fT1wr = new RSTPC_T1wrapper( fOutFile );
			fT1wr -> SetFileOwner(false);
		}
		else
		{//Check if the right file is open and if the class is properly initialised
			if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): T1 wrapper already allocated. Checking if it is consistent for the T2 processing." << endl;
			
			if( !fT1wr->fInfile )
			{
				if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): The pointer to input file in the T1 wrapper is 0. Delete the object and reallocate it." << endl;
				delete fT1wr;
				fT1wr = new RSTPC_T1wrapper( fOutFile );
				fT1wr -> SetFileOwner(false);
			}
			else if( !fT1wr->fInfile->IsWritable() )
			{
				if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): The input file in the T1 wrapper is not writable. Delete the object and reallocate it." << endl;
				
				if( (fT1wr->IsFileOwner()) && (fT1wr->fInfile==fOutFile) ) fT1wr -> SetFileOwner(false);
				delete fT1wr;
				fT1wr = new RSTPC_T1wrapper( fOutFile );
			}
			else if( string(fT1wr->fInfile->GetName())!=string(fOutFile->GetName()) )
			{
				if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): The input file in the T1 wrapper is not the same as the one that should be processed for T2 production. Delete the object and reallocate it." << endl;
				if( (fT1wr->IsFileOwner()) && (fT1wr->fInfile==fOutFile) ) fT1wr -> SetFileOwner(false);
				delete fT1wr;
				fT1wr = new RSTPC_T1wrapper( fOutFile );
			}
			else if( !fT1wr->IsInit() )
			{
				if(fgDebug) cout << "Debug --> RSTPC_RunProcessor::InitT2proc(): The T1 wrapper is present but not initialised. Initialising it." << endl;
				TTree* T1 = (TTree*)fOutFile->Get("T1");
				if(!T1)
				{
					cout << "ERROR --> RSTPC_RunProcessor::InitT2proc(): Cannot find the T1 tree in file <" << fOutFile->GetName() << ">. Deleteing the T1 wrapper object." << endl;
					if( (fT1wr->IsFileOwner()) && (fT1wr->fInfile==fOutFile) ) fT1wr -> SetFileOwner(false);
					delete fT1wr;
					fT1wr = NULL;
					return false;
				}
				fT1wr -> SetFileOwner(false);
				fT1wr->Init( T1 );
			}
		}
		//At this point the T1 wrapper should be initialised
		
		if(!fT1wr->IsInit())
		{
			cout << "ERROR --> RSTPC_RunProcessor::InitT2proc(): Cannot initialise the T1 wrapper with the file <" << fOutFile->GetName() << ">. Deleteing the T1 wrapper object." << endl;
			if( (fT1wr->IsFileOwner()) && (fT1wr->fInfile==fOutFile) ) fT1wr -> SetFileOwner(false);
			delete fT1wr;
			fT1wr = NULL;
			return false;
		}
	}
	
	if(gEventData) delete gEventData;
	gEventData = new EventData;
	
	fProcT2 = true;
	return true;
}


void RSTPC_RunProcessor::DescribeT2()
{
	//fProcT2 = false;
	
	Bool_t procflag = false;
	if(!fProcT2)
	{
		if(!InitT2proc())
		{
			fProcT2 = false;
			return;
		}
	}
	
	if(!fOutFile) fOutFile = fT1wr->fInfile; //Maybe redundant
	
	fOutFile->cd(); //Maybe redundant
	
	fOutT2 = new TTree("T2","Tear 2 data of the merged tree for the Resistive Shell TPC data");
	
	if(!TClassTable::GetDict("RSTPC_Hits")) {
		gSystem->Load("RSTPC_Hits");
		gROOT->ProcessLine("#include \"RSTPC_Hits.hh\"");
	}
	
	fOutT2->Branch("GoodEvent",&gEventData->GoodEvent,"GoodEvent/O");
	//fOutT2->Branch("ColPulses", "TObjArray", &gEventData->ColPulses, 32000, 2);
	//fOutT2->Branch("IndPulses", "TObjArray", &gEventData->IndPulses, 32000, 2);
	//fOutT2->Branch("Hits", "TObjArray", &gEventData->Hits, 32000, 2);
	fOutT2->Branch("ColPulses", "TClonesArray", &gEventData->ColPulses, 32000, 2);
	fOutT2->Branch("IndPulses", "TClonesArray", &gEventData->IndPulses, 32000, 2);
	fOutT2->Branch("Hits", "TClonesArray", &gEventData->Hits, 32000, 2);
	
	fProcT2 = true;
}


void RSTPC_RunProcessor::T2Process()
{
	if(!fProcT2) return;
	
	
	Int_t nEvs = fT1wr->fChain->GetEntries();
	
	
	map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >* pulsesMap = NULL;
	
	if(fgDebug)
	{
		cout << "Debug --> RSTPC_RunProcessor::T2Process(): Number of events to process: " << nEvs << endl;
	}
	
	for(Int_t iEv=0; iEv<nEvs; iEv++)
	{
		//Clear and or reset all the structures
		
		if(fgDebug)
		{
			if(iEv%fgNevPrint==0)
			{
				cout << "\n====================================================================================" << endl;
				cout << "Debug --> RSTPC_RunProcessor::T2Process(): Processing event " << iEv << ". Resetting EventData class and counters." << endl;
			}
		}
		
		gEventData->Reset("C");
		RSTPC_Pulse::ResetCounter();
		RSTPC_Hit::ResetCounter();
		gColPulses->clear();
		gIndPulses->clear();
		gHits->clear();
		
		if(fgDebug)
		{
			if(iEv%fgNevPrint==0)
			{
				cout << "Debug --> RSTPC_RunProcessor::T2Process(): Reading event from the T1 wrapper." << endl;
			}
		}
		
		fT1wr->GetEntry(iEv);
		
		if( !(fT1wr->ColHist && fT1wr->IndHist) )
		{
			cout << "ERROR --> RSTPC_RunProcessor::T2Process(): Either the histogram for the collection or the induction wires is missing. Event " << iEv << " is broken." << endl;
			gEventData->GoodEvent = false;
		}
		
		//Find all the pulses of the collection and induction wires
		if(fgDebug)
		{
			if(iEv%100==0)
			{
				cout << "Debug --> RSTPC_RunProcessor::T2Process(): Finding pulses for the collection wires." << endl;
			}
		}
		FindPulses(fT1wr->ColHist, kCol, fgDebug&&(iEv%fgNevPrint==0) );
		
		
		if(fgDebug)
		{
			if(iEv%fgNevPrint==0)
			{
				cout << "Debug --> RSTPC_RunProcessor::T2Process(): Finding pulses for the induction wires." << endl;
			}
		}
		FindPulses(fT1wr->IndHist, kInd, fgDebug&&(iEv%fgNevPrint==0));
		
		
		//Now find which induction pulses correspond to each collection pulses
		if(fgDebug)
		{
			if(iEv%fgNevPrint==0)
			{
				cout << "Debug --> RSTPC_RunProcessor::T2Process(): Found " << gColPulses->size() << " collection pulses and " << gIndPulses->size() << " induction pulses. Combining them." << endl;
			}
		}
		pulsesMap = CombinePulses(gColPulses, gIndPulses, fgDebug&&(iEv%fgNevPrint==0)); //The key is a collection pulse the value is the vector of the corresponding induction pulse
		
		if( !pulsesMap )
		{
			cout << "ERROR --> RSTPC_RunProcessor::T2Process(): The \"pulsesMap\" pointer is 0. Event " << iEv << " is broken." << endl;
		}
		
		if(fgDebug)
		{
			if(iEv%fgNevPrint==0)
			{
				cout << "Debug --> RSTPC_RunProcessor::T2Process(): Finding the hits...." << endl;
			}
		}
		
		Int_t nHits = HitsFinder(pulsesMap, fgDebug&&(iEv%fgNevPrint==0));
		
		if(fgDebug && (iEv%fgNevPrint==0))
		{
			cout << "Debug --> RSTPC_RunProcessor::T2Process(): " << nHits <<" hits found." << endl;
		}
		
		
		
		//Fill the TClonesArray before saving the tree
		
		gEventData->ColPulses->ExpandCreate(gColPulses->size());
		TClonesArray &arr = *(gEventData->ColPulses);
		for(Int_t iPulse=0; iPulse<gColPulses->size(); iPulse++)
		{
			new(arr[iPulse]) RSTPC_Pulse(*gColPulses->at(iPulse));
		}
		
		gEventData->IndPulses->ExpandCreate(gIndPulses->size());
		arr = *(gEventData->IndPulses);
		for(Int_t iPulse=0; iPulse<gIndPulses->size(); iPulse++)
		{
			new(arr[iPulse]) RSTPC_Pulse(*gIndPulses->at(iPulse));
		}
		
		gEventData->Hits->ExpandCreate(gHits->size());
		arr = *(gEventData->Hits);
		for(Int_t iHit=0; iHit<gHits->size(); iHit++)
		{
			new(arr[iHit]) RSTPC_Hit(*gHits->at(iHit));
		}
		
		
		fOutT2->Fill();
		
		
		//Clean the pulses map in order to avoid memory leaks
		if(pulsesMap)
		{
			map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >::iterator mapIt;
			for(mapIt=pulsesMap->begin(); mapIt!=pulsesMap->end(); ++mapIt)
			{
				vector<RSTPC_Pulse*>::iterator vecIt;
				for(vecIt=mapIt->second->begin(); vecIt!=mapIt->second->end(); ++vecIt)
				{
					if(*vecIt) delete *vecIt;
				}
				if(mapIt->second)
				{
					delete mapIt->second;
				}
			}
			delete pulsesMap;
			pulsesMap = NULL;
		}
	}//End cycling over the events
	
	fOutFile->WriteTObject(fOutT2,0,"overwrite");
	fOutT2->ResetBranchAddresses();
	
	if(gEventData)
	{
		delete gEventData;
		gEventData=NULL;
	}
	
	if(gColWfsVect) delete gColWfsVect; gColWfsVect=NULL;
	if(gIndWfsVect) delete gIndWfsVect; gIndWfsVect=NULL;
}


void RSTPC_RunProcessor::FindPulses(TH2D* h, WireType type, Bool_t debug)
{
	if(!h) return;
	
	//Clear the array of the pulses
	if( (type!=kCol)&&(type!=kInd) )
	{
		return;
	}
	
	Int_t nChs = h->GetNbinsY(); //Number of channels
	Int_t nSamps = h->GetNbinsX(); //Number of time samples
	
	if(debug)
	{
		cout << "Debug -->  RSTPC_RunProcessor::FindPulses(...): From histogram \""<< h->GetName() << "\" determined " << nChs << " channels and " << nSamps << " time samples." << endl;
	}
	
	vector<Double_t> *flWf = new vector<Double_t>(nSamps);
	
	
	if(type==kCol)
	{
		if(!gColWfsVect)
		{
			gColWfsVect = new vector<vector<Double_t> >(nChs);
		}
		else if(gColWfsVect->size()!=nChs)
		{
			gColWfsVect->resize(nChs);
		}
		
		for(Int_t iCh=0; iCh<nChs; iCh++)
		{
			if((gColWfsVect->at(iCh)).size()!=nSamps) (gColWfsVect->at(iCh)).resize(nSamps);
		}
	}
	
	if(type==kInd)
	{
		if(!gIndWfsVect)
		{
			gIndWfsVect = new vector<vector<Double_t> >(nChs);
		}
		else if(gIndWfsVect->size()!=nChs)
		{
			gIndWfsVect->resize(nChs);
		}
		
		for(Int_t iCh=0; iCh<nChs; iCh++)
		{
			if((gIndWfsVect->at(iCh)).size()!=nSamps) (gIndWfsVect->at(iCh)).resize(nSamps);
		}
	}
	
	
	if(debug)
	{
		cout << "Debug -->  RSTPC_RunProcessor::FindPulses(...): Vector table of the waveforms initialised" << endl;
	}
	
	
	for(Int_t iCh=0; iCh<nChs; iCh++)
	{
		TH1D* hChWf = h->ProjectionX("hChWf",iCh+1,iCh+1);
		Double_t thr = fgSigmaThr*((fT1wr->RmsColWires)[iCh]);
		
		if( hChWf->GetNbinsX()!=nSamps )
		{
			cout << "ERROR --> RSTPC_RunProcessor::FindPulses(...): The number of samples in the histogram \"" << hChWf->GetName() << "\" are " << hChWf->GetNbinsX() << " while they should be " << nSamps << endl;
			return;
		}
		
		//Find the pulses on each channel
		if(type==kCol)
		{
			if( (gColWfsVect->at(iCh)).size()!=nSamps )
			{
				cout << "ERROR --> RSTPC_RunProcessor::FindPulses(...): The number of slots in the \"gColWfsVect\" table for channel" << iCh <<" are " << (gColWfsVect->at(iCh)).size() << " while they should be " << nSamps << endl;
				return;
			}
			
			//First fill up the waveform container
			if(flWf->size()!=nSamps) flWf->resize(nSamps);
			for(Int_t iSamp=0; iSamp<nSamps; iSamp++)
			{
				Double_t val = hChWf->GetBinContent(iSamp+1);
				(gColWfsVect->at(iCh)).at(iSamp) = val;
				if( val>(fgSigmaThr*((fT1wr->RmsColWires)[iCh])) )
				{
					flWf->at(iSamp) = hChWf->GetBinContent(iSamp+1);
				}else{
					flWf->at(iSamp) = 0;
				}
			}
			
			//Find the pulses on the channel iCh and fill the corresponding objects
			nSamps = flWf->size();
			for(Int_t iSamp=0; iSamp<nSamps; iSamp++)
			{
				if( (iSamp<nSamps) && flWf->at(iSamp)>thr )
				{
					Int_t ledge = iSamp;
					Int_t redge = ledge;
					Double_t max = flWf->at(iSamp);
					Int_t maxpos = iSamp;
					iSamp++;
					while( (iSamp<nSamps) && (flWf->at(iSamp)>0) )
					{
						if(max<flWf->at(iSamp))
						{
							max=flWf->at(iSamp);
							maxpos = iSamp;
						}
						iSamp++;
						
						if(iSamp<nSamps)
						{
							if(flWf->at(iSamp)==0) redge = iSamp-1;
						}
						else
						{
							redge = iSamp-1;
						}
					}
					if(redge>ledge+1)
					{//Make a pulse object and fill it up with the positive part of the pulse
						RSTPC_Pulse *pulse = new RSTPC_Pulse(kCol);
						pulse->fWireNum = iCh;
						pulse->fLedge = ledge;
						pulse->fRedge = redge;
						pulse->fMax = max;
						pulse->fMaxPos = maxpos;
						pulse->fMin = 0; //This might change if the negative part of the pulse is found
						pulse->fMinPos = -1; //This might change if the negative part of the pulse is found
						
						iSamp = pulse->fRedge+1; //The loop restarts from the sample after the end of the loop
						
						//Add the pulse to the array that will be saved into the tree
						//gEventData->ColPulses->Add(pulse);
						
						pulse->SetFWHM(flWf);
						pulse->SetFWTM(flWf);
						pulse->SetMeanTime(flWf);
						pulse->SetSigma(flWf);
						
						gColPulses->push_back(pulse);
					}
					//if(iSamp>=(flWf->size()-1)) break; //This exits also from the for loop
				}//End of the pulse condition (flWf->at(iSamp)>thr)
			}//End loop over the time samples (Collection wire)
		}
		
		if(type==kInd)
		{//Here the pulse has (ideally) positive and negative part
			if( (gIndWfsVect->at(iCh)).size()!=nSamps )
			{
				cout << "ERROR --> RSTPC_RunProcessor::FindPulses(...): The number of slots in the \"gIndWfsVect\" table for channel" << iCh <<" are " << (gIndWfsVect->at(iCh)).size() << " while they should be " << nSamps << endl;
				return;
			}
			
			//First fill up the waveform container
			if(flWf->size()!=nSamps) flWf->resize(nSamps);
			for(Int_t iSamp=0; iSamp<nSamps; iSamp++)
			{
				Double_t val = hChWf->GetBinContent(iSamp+1);
				(gIndWfsVect->at(iCh)).at(iSamp) = val;
				if( abs(val)>thr )
				{
					flWf->at(iSamp) = hChWf->GetBinContent(iSamp+1);
				}else{
					flWf->at(iSamp) = 0;
				}
			}
			
			
			//Find the pulses on the channel iCh and fill the corresponding objects
			nSamps = flWf->size();
			for(Int_t iSamp=0; iSamp<nSamps; iSamp++)
			{
				if( (iSamp<nSamps) && flWf->at(iSamp)>thr )
				{//The left edge is positive fluctuation
					Int_t ledge = iSamp;
					Int_t redge = ledge;
					Double_t max = flWf->at(iSamp);
					Int_t maxpos = iSamp;
					iSamp++;
					while( (iSamp<flWf->size()) && (flWf->at(iSamp)>0) )
					{
						if(max<flWf->at(iSamp))
						{
							max=flWf->at(iSamp);
							maxpos = iSamp;
						}
						
						iSamp++;
						if(iSamp<flWf->size())
						{
							if(flWf->at(iSamp)==0) redge = iSamp-1;
						}
						else
						{
							redge = iSamp-1;
						}
					}
					if(redge>ledge+1)
					{//A valid pulse has been found. Make a pulse object and fill it up with the positive part of the pulse
						RSTPC_Pulse *pulse = new RSTPC_Pulse(kInd);
						pulse->fWireNum = iCh;
						pulse->fLedge = ledge;
						pulse->fRedge = redge;
						pulse->fMax = max;
						pulse->fMaxPos = maxpos;
						pulse->fMin = 0; //This might change if the negative part of the pulse is found
						pulse->fMinPos = -1; //This might change if the negative part of the pulse is found
						
						//Try to find the negative part of the pulse
						Int_t searchWin = (Int_t)(GetPeakingTime()*GetSamplingFreq()+0.5); //Max extension where I look for the left edge of the negative part of the pulse
						Int_t negledge = redge;
						Int_t negredge = negledge;
						Double_t min = 0;
						Int_t minpos = -1;
						for(Int_t jSamp=redge; (jSamp<nSamps) && (jSamp<=redge+searchWin); jSamp++)
						{
							//if(jSamp<=redge+searchWin) break; //This terminate the loop
							
							if( (jSamp<nSamps) && (flWf->at(jSamp)<0) )
							{//Possible start of the negative part of the induction pulse
								negledge = jSamp;
								negredge = negledge;
								jSamp++;
								while( (jSamp<nSamps) && (flWf->at(jSamp)<0) )
								{
									if(min>flWf->at(jSamp))
									{
										min=flWf->at(jSamp);
										minpos = jSamp;
									}
									jSamp++;
									
									if(jSamp<nSamps)
									{
										if(flWf->at(jSamp)==0) negredge = jSamp-1;
									}
									else
									{
										negredge = jSamp-1;
									}
								}
								//Check that the negative part of the pulse is ok and then save it
								if(negredge > negledge+1)
								{
									pulse->fMin = min;
									pulse->fMinPos = minpos;
									pulse->fRedge = negredge;
									jSamp = nSamps; //This terminate the loop
								}
							}
						}//End of the loop to search for the negative part of the pulse
						iSamp = pulse->fRedge+1; //The loop restarts from the sample after the end of the loop
						//gEventData->IndPulses->Add(pulse);
						
						pulse->SetMeanTime(flWf);
						pulse->SetSigma(flWf);
						
						gIndPulses->push_back(pulse);
					}//Exit from the pulse making scope
					
				}//End of pulse serching
			}//End loop over the time samples (Induction wire)
		}
	}//End of the loop over the channels
	
	//Avoid memory leaks
	if(gROOT->FindObject("hChWf")) delete gROOT->FindObject("hChWf");
	if(flWf) delete flWf;
}


map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >* RSTPC_RunProcessor::CombinePulses(vector<RSTPC_Pulse*>* ColPulses,  vector<RSTPC_Pulse*>* IndPulses, Bool_t debug)
{
	if( !(ColPulses && IndPulses) ) return NULL;
	if( !((ColPulses->size()>0) && (IndPulses->size()>0)) ) return NULL; //They must be empty!!!
	
	Int_t nMatches = 0;
	
	map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* > *outmap = new map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >;
	
	Int_t nColPulses = ColPulses->size();
	Int_t nIndPulses = IndPulses->size();
	
	//This is the maximum shift time (in samples) that two pulses must have in order to be "coincident".
	//It should correspond to the electronics peaking time which also defines the time slices of the tpc;
	//Int_t timeWindow = (Int_t)(GetPeakingTime()*GetSamplingFreq()+0.5); //Trick to round to the closest integer
	
	
	vector<RSTPC_Pulse*>::iterator cVecIt, iVecIt;
	
	for(cVecIt=ColPulses->begin(); cVecIt!=ColPulses->end(); cVecIt++ )
	{
		vector<RSTPC_Pulse*> *IndPulsesVec = new vector<RSTPC_Pulse*>;
		
		RSTPC_Pulse *ColPulse = (*cVecIt);
		
		for(iVecIt=IndPulses->begin(); iVecIt!=IndPulses->end(); iVecIt++ )
		{
			RSTPC_Pulse *IndPulse = (*iVecIt);
			
			if( !( (ColPulse->fLedge>IndPulse->fRedge)||(ColPulse->fRedge<IndPulse->fLedge) ) )
			{//This is the overlapping condition
				IndPulsesVec->push_back(IndPulse);
			}
			
			if(IndPulsesVec->size()>0)
			{
				if(debug) nMatches += IndPulsesVec->size();
			}
			(*outmap)[ColPulse] = IndPulsesVec;
		}
	}
	
	if(debug) cout << "Debug --> RSTPC_RunProcessor::CombinePulses(...): Total matches found: " << nMatches << endl;
	
	return outmap;
}


Int_t RSTPC_RunProcessor::HitsFinder(map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >* pulseMap, Bool_t debug)
{
	gHits->clear();
	
	//Make time slices and find the hits
	Int_t timeSlice = (Int_t)(GetPeakingTime()*GetSamplingFreq()+0.5); //Trick to round to the closest integer
	Int_t nSlices = ceil( GetSamplingFreq()*GetDriftLenght()/GetDriftVel()/timeSlice );
	
	if(debug)
	{
		cout << "Debug --> RSTPC_RunProcessor::HitsFinder(...): Using " << nSlices << " to determine the hits." << endl;
		cout << "Debug --> RSTPC_RunProcessor::HitsFinder(...): Time slice: " << timeSlice << " samples = " << timeSlice/GetSamplingFreq() << " usecs = " << GetDriftVel()*timeSlice/GetSamplingFreq() << " mm." << endl;
	}
	
	vector<Int_t> TimeEdges(nSlices+1);
	for(Int_t iSl=0; iSl<=nSlices; iSl++)
	{
		TimeEdges.at(iSl) = iSl*timeSlice;
	}
	
	for(Int_t iSl=0; iSl<nSlices; iSl++)
	{
		Int_t ledge = TimeEdges.at(iSl);
		Int_t redge = TimeEdges.at(iSl+1);
		Int_t nsamps = redge-ledge+1;
		
		
		//if(debug)
		if(false)
		{
			cout << "Debug --> RSTPC_RunProcessor::HitsFinder(...): Slice " << iSl << ": ledge=" << ledge << "\tredge=" << redge << endl;
		}
		
		map<RSTPC_Pulse*, vector<RSTPC_Pulse*>* >::iterator mapIt;
		for(mapIt=pulseMap->begin(); mapIt!=pulseMap->end(); mapIt++)
		{
			RSTPC_Pulse *ColPulse = mapIt->first;
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
					sum2 += pow((gColWfsVect->at(ColPulse->fWireNum)).at(kSamp) ,2);
					meantime += kSamp* pow((gColWfsVect->at(ColPulse->fWireNum)).at(kSamp) ,2); 
					meanheight += (gColWfsVect->at(ColPulse->fWireNum)).at(kSamp)/nsamps;
				}
				
				meantime = meantime/sum2;
				
				if(meantime==0) meantime = ((Double_t)(redge-ledge))/2;
				
				hit->fMeanTime = meantime;
				hit->fMeanHeight = meanheight;
				
				
				gHits->push_back(hit);
				
				//gEventData->Hits->Add(hit);
				
				//if(debug)
				if(false)
				{
					cout << "           Hit ID: " << hit->fHitID << "->\tledge=" << hit->fLedge << "\tredge=" << hit->fRedge << "\tcentre time="<< hit->fCentreTime <<"\tmeantime=" << meantime << "\tmeanheight=" << meanheight << endl;
				}
				else
				{
					if(hit->fCentreTime <= 0)
					{
						cout << "WARNING ---> RSTPC_RunProcessor::HitsFinder(...): Non positive centre time!" << endl;
						cout << "           Hit ID: " << hit->fHitID << "->\tledge=" << hit->fLedge << "\tredge=" << hit->fRedge << "\tcentre time="<< hit->fCentreTime <<"\tmeantime=" << meantime << "\tmeanheight=" << meanheight << endl;
					}
					else if(meantime<=0)
					{
						cout << "WARNING ---> RSTPC_RunProcessor::HitsFinder(...): Non positive mean time!" << endl;
						cout << "           Hit ID: " << hit->fHitID << "->\tledge=" << hit->fLedge << "\tredge=" << hit->fRedge << "\tcentre time="<< hit->fCentreTime <<"\tmeantime=" << meantime << "\tmeanheight=" << meanheight << endl;
					}
				}
			}
			
		}
	}
	
	return gHits->size();
}


/*
void RSTPC_RunProcessor::CloseRun()
{
	if(fOutT1)
	{
		if(fOutT1->GetCurrentFile()) delete fOutT1->GetCurrentFile();
		fOutFile = NULL;
		fOutT1 = NULL;
	}
	
	if(fTpcMan) delete fTpcMan;
	fTpcMan = NULL;
	
	if(fFebMan) delete fFebMan;
	fFebMan = NULL;
}
*/


EventData::EventData()
{
	ColPulses = new TClonesArray("RSTPC_Pulse");
	//ColPulses->SetOwner(kTRUE);
	
	IndPulses = new TClonesArray("RSTPC_Pulse");
	//IndPulses->SetOwner(kTRUE);
	
	Hits = new TClonesArray("RSTPC_Hit");
	//Hits->SetOwner(kTRUE);
	
	GoodEvent = true;
}

EventData::~EventData()
{
	if(ColPulses) delete ColPulses;
	if(IndPulses) delete IndPulses;
	if(Hits) delete Hits;
}

void EventData::Reset(string opt)
{
	if(ColPulses) ColPulses->Clear(opt.c_str());
	if(IndPulses) IndPulses->Clear(opt.c_str());
	if(Hits) Hits->Clear(opt.c_str());
	
	GoodEvent = true;
}


#endif /* RSTPC_RUNPROCESSOR_CC */
