#ifndef RSTCP_ANALYSER_CC
#define RSTCP_ANALYSER_CC

#include "RSTPC_Analyser.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"

#include <string>
#include <sstream>
#include <fstream>
#include <map>



using namespace std;



RSTPC_Analyser::RSTPC_Analyser()
{
	fCanv = NULL;
	//char str[128];
	hI0 = NULL;
	hC0 = NULL;
	hI = NULL;
	hC = NULL;
	hII = NULL;
	hE = NULL;
	hU = NULL;
	hR = NULL;
	
	fRunNumber = -1; //Negative number is also a flag to simbolyse that no run is loaded
	fInfile = NULL;
	//FileStat_t stat;
	scaled = 0;
	pattern = 0;
	fTrigTree = NULL;
	adc_to_mV = 0.1512; //mV per ADC
	//#define SCALE 1
	
	
	if( !RSTPC_Options::GetInstance() )
	{
		fDataDir = "/home/francescop/data/ResistiveShell/";
		fOutDir = "/home/francescop/ArCube/analysis/ResistiveShellTPC/";
	}
	else if( RSTPC_Options::GetInstance()->IsDataDirSet() )
	{
		fDataDir = RSTPC_Options::GetInstance()->GetDataDir();
		if(fDataDir.size())
		{
			if(fDataDir.at(fDataDir.size()-1)!='/') fDataDir += "/";
		}
		else
		{
			fDataDir = "./";
		}
	}
	
	
	fCMrej = true;
	fCMrejIter = 1;
	
	//bgx = NULL;
	//bgy = NULL;
	
	fXlowBin = -1;
	fXupBin = -1;
	
	fPrintFlag = true;
	
	fCollectionMap = NULL;
	fInductionMap = NULL;
	
	
	fColRMS = new vector<Double_t>;
	fIndRMS = new vector<Double_t>;
	
	fColRelGains = new vector<Double_t>(32, 1.);
	fIndRelGains = new vector<Double_t>(32, 1.);
	
	fSigmaThr = 3.0;
	
	//Open the file with the data and read it
	//fLoadedRun = OpenRun(_RunNum);
	//if(!fLoadedRun) fRunNumber = -1;
}


RSTPC_Analyser::~RSTPC_Analyser()
{
	if(fTrigTree)
	{
		fTrigTree->ResetBranchAddresses();
		delete fTrigTree->GetCurrentFile();
	}
	else
	{
		if(fInfile) delete fInfile;
	}
	
	if(fCollectionMap) delete fCollectionMap;
	if(fInductionMap) delete fInductionMap;
	
	if(fColRMS) delete fColRMS;
	if(fIndRMS) delete fIndRMS;
	
	if(fColRelGains) delete fColRelGains;
	if(fIndRelGains) delete fIndRelGains;
	
	fTrigTree = NULL;
	fInfile = NULL;
	fCollectionMap = NULL;
	fInductionMap = NULL;
	fColRMS = NULL;
	fIndRMS = NULL;
}


void RSTPC_Analyser::LoadCollMap(string filename)
{
	ifstream CollMapFile(filename.c_str());
	if(!CollMapFile)
	{
		cout << "\nERROR: Could not find/open the file <" << filename << "> for loading the channel map for the collection wires\n" << endl;
		return;
	}
	
	if(!fCollectionMap)
	{
		fCollectionMap = new map<Int_t, Int_t>;
	}
	else
	{
		fCollectionMap->clear();
	}
	
	stringstream line;
	string str;
	while( getline(CollMapFile, str) )
	{
		line.clear(); line.str(""); line << str;
		
		line >> str;
		Int_t key = atoi(str.c_str()); //This is the actual wire number
		
		line >> str;
		Int_t val = atoi(str.c_str()); //This is the channel
		
		(*fCollectionMap)[key] = val;
		
	}
}

void RSTPC_Analyser::LoadIndcMap(string filename)
{
	ifstream IndcMapFile(filename.c_str());
	if(!IndcMapFile)
	{
		cout << "\nERROR: Could not find/open the file <" << filename << "> for loading the channel map for the induction wires\n" << endl;
		return;
	}
	
	if(!fInductionMap)
	{
		fInductionMap = new map<Int_t, Int_t>;
	}
	else
	{
		fInductionMap->clear();
	}
	
	stringstream line;
	string str;
	while( getline(IndcMapFile, str) )
	{
		line.clear(); line.str(""); line << str;
		
		line >> str;
		Int_t key = atoi(str.c_str()); //This is the actual wire number
		
		line >> str;
		Int_t val = atoi(str.c_str()); //This is the channel
		
		(*fInductionMap)[key] = val;
		
	}
}


Bool_t RSTPC_Analyser::OpenRun(Int_t RunNumber)
{
	stringstream ss_tmp;
	
	fLoadedRun = false;
	
	if(RunNumber<0) return false;
	
	ss_tmp.str(""); ss_tmp << fDataDir << "Run_" << setfill('0') << setw(9) << RunNumber << ".root";
	cout << "Opening file <" << ss_tmp.str() << ">" << endl;
	fInfile = TFile::Open( ss_tmp.str().c_str(), "read" );
	if( (!fInfile) || (!fInfile->IsOpen()) )
	{
		cout << "\nERROR: File <" << ss_tmp.str() << "> can't be opened.\n" << endl;
		fRunNumber = -1;
		return false;
	}
	
	fRunNumber=RunNumber;
	cout << "<" << fInfile->GetName() << "> opened.\n" << endl;
	
	if(fTrigTree)
	{
		delete fTrigTree;
		fTrigTree = NULL;
	}
	
	fTrigTree = (TTree*)fInfile->Get("Ev_Head");
	if(!fTrigTree)
	{
		cout << "\nWARNING: Could not find/load tree \"Ev_Head\" in file <" << fInfile->GetName() << ">. The event timestamps are not accessible for this run.\n" << endl;
	}
	else
	{
		fTrigTree->SetBranchAddress("EventTime", &fEventTime);
	}
	
	
	//Find the channel numbers for the collection and for the induction wires
	Int_t nChsCol = 0;
	TH2 *hCol = (TH2*)fInfile->Get("Col_0");
	if(hCol) nChsCol = hCol->GetNbinsY();
	
	Int_t nChsInd = 0;
	TH2 *hInd = (TH2*)fInfile->Get("Ind_0");
	if(hInd) nChsInd = hInd->GetNbinsY();
	
	if(nChsCol)
	{
		fColRMS->assign(nChsCol, 0.);
	}
	
	if(nChsInd)
	{
		fIndRMS->assign(nChsInd, 0.);
	}
	
	fLoadedRun = true;
	
	return true;
}


void RSTPC_Analyser::LoadEvent(Int_t evNumber)
{
	stringstream sstr;
	
	if( (!fLoadedRun) || (evNumber<0) || (!fInfile) || (!fInfile->IsOpen()) ) return;
	
	sstr.str(""); sstr << "Col_" << evNumber;
	if(hC0){delete hC0; hC0 = NULL;}
	hC0 = (TH2D*)fInfile->Get( sstr.str().c_str() );
	if(!hC0)
	{
		cout << "WARNING: Could not find object \"" << sstr.str() << "\" in file <" << fInfile->GetName() << ">" << endl;
	}
	
	sstr.str(""); sstr << "Ind_" << evNumber;
	if(hI0){delete hI0; hI0 = NULL;}
	hI0 = (TH2D*)fInfile->Get( sstr.str().c_str() );
	if(!hI0)
	{
		cout << "WARNING: Could not find object \"" << sstr.str() << "\" in file <" << fInfile->GetName() << ">" << endl;
	}
	
	return;
}


//Staff for proper scaling
TH2I *sI,*sC;
TF1 *F1;
TF1 *f;
Double_t xbins[8198];
Double_t ybins[64];
void RSTPC_Analyser::scale()
{
	Double_t McsPerBin=0.051;
	Double_t Length=0.15; //m
	Double_t Pitch=6.0/32.; //wire pitch
	Int_t BinAtAnod=120;
	Int_t BinAtCath=2000;

	for(int j=0;j<4100;j++) xbins[j]=j*McsPerBin;
	for(int j=0;j<=64;j++) ybins[j]=j*Pitch;

	scaled=1;

	cout << "Scale converted to cm: " << endl;
	cout << "Anode : Bin " << BinAtAnod << " Time " << xbins[BinAtAnod] << endl;
	cout << "Cathode : Bin " << BinAtCath << " Time " << xbins[BinAtCath] << endl;
}



//Adjusting the maps for the wires
void RSTPC_Analyser::ApplyWiresMaps(TH2 *hc, TH2 *hi)
{
	if( fCollectionMap && (fCollectionMap->size()==hc->GetNbinsY()) )
	{
		TH2 *hcopy = (TH2*)hc->Clone("hcopy");
		
		Int_t nXbins = hc->GetNbinsX();
		
		map<Int_t, Int_t>::iterator it;
		for(it=fCollectionMap->begin(); it!=fCollectionMap->end(); it++)
		{
			Int_t iWire = it->first;
			Int_t iCh = it->second;
			
			for(Int_t xBin=1; xBin<=nXbins; xBin++)
			{
				hc->SetBinContent(xBin, iWire, hcopy->GetBinContent(xBin, iCh) );
			}
		}
		if(hcopy) delete hcopy;
	}
	
	if(fInductionMap && (fInductionMap->size()==hi->GetNbinsY()) )
	{
		TH2 *hcopy = (TH2*)hi->Clone("hcopy");
		
		Int_t nXbins = hi->GetNbinsX();
		
		map<Int_t, Int_t>::iterator it;
		for(it=fInductionMap->begin(); it!=fInductionMap->end(); it++)
		{
			Int_t iWire = it->first;
			Int_t iCh = it->second;
			
			for(Int_t xBin=1; xBin<=nXbins; xBin++)
			{
				hi->SetBinContent(xBin, iWire, hcopy->GetBinContent(xBin, iCh) );
			}
		}
		if(hcopy) delete hcopy;
	}
}


//Stuff for common mode noise suppression
void RSTPC_Analyser::BaselineCorr(TH2* h, Int_t LowBinX, Int_t UpBinX, vector<Double_t>* RMS_vect)
{
	if(!h) return;
	
	//The projection is performed by summing over the integrated direction (here summing over the X variable)
	//Here I can find the baseline (or pedestal) for each channel
	if( (LowBinX<0) || (UpBinX<0) )
	{//In this case use the entire window to calculate the baseline
		LowBinX = 0;
		UpBinX = h->GetNbinsX();
	}
	
	if( (UpBinX>h->GetNbinsX()) || UpBinX<LowBinX ) UpBinX = h->GetNbinsX();
	
	Int_t nChs = h->GetNbinsY(); //Number of channels
	//Int_t nSamps = h->GetNbinsX()
	
	//if(bgy){delete bgy; bgy=NULL;}
	TH1D* bgy = h->ProjectionY("bgy", LowBinX, UpBinX);
	bgy->Scale(1./(UpBinX-LowBinX+1));
	
	
	//Shift the waveform of each channel to zero baseline
	for (int x=0; x<h->GetNbinsX(); x++)
	{
		for (int y=0; y<nChs; y++)
		{
			Double_t bsln = bgy->GetBinContent(y);
			h->SetBinContent(x,y,h->GetBinContent(x,y)-bsln);
		}
	}
	
	
	//After the baseline correction calculate the RMS
	//Be careful that if common mode is not cancelled the RMS will be affected by that
	if(!RMS_vect)
	{
		if(bgy) delete bgy;
		return;
	}
	
	if(RMS_vect->size()!=nChs) RMS_vect->assign(nChs,0);
	
	vector<Double_t> samples(UpBinX-LowBinX+1);
	for(Int_t iCh=0; iCh<nChs; iCh++)
	{
		for (Int_t xBin=LowBinX; xBin<=UpBinX; xBin++)
		{
			samples.at(xBin-LowBinX) = h->GetBinContent(xBin, iCh);
		}
		
		RMS_vect->at(iCh) = TMath::RMS(samples.size(), &samples.at(0) );
	}
	
	if(bgy) delete bgy;
	
	return;
}


TH1D* RSTPC_Analyser::CMrej(TH2* h, vector<Double_t>* RMS_vect, Int_t iterations, Bool_t bslncorr)
{
	stringstream sstr;
	iterations--;
	
	if(!h) return NULL;
	if(!RMS_vect)
	{
		cout << "\nERROR: Without the vector of the RMS values the CM noise reduction is not possible. No CM rejection will be applied for the histogram \"" << h->GetName() << "\"!\n" << endl;
		return NULL;
	}
	
	//With the channels shifted to zero baseline now I can compute and subtract the common mode noise
	//The projection is performed by summing over the direction integrated out (here summing over the Y variable)
	//With this I find the common mode induced oscillations (noise)
	sstr.str(""); sstr << h->GetName() << "_corr";
	TH1D *h_corr = (TH1D*)gROOT->FindObject( sstr.str().c_str() );
	if(!h_corr)
	{
		h_corr = (TH1D*)h->ProjectionX( sstr.str().c_str() );
		h_corr->Reset();
	}
	
	Int_t nChs = h->GetNbinsY(); //Number of channels
	
	
	for (Int_t xBin=1; xBin<h->GetNbinsX(); xBin++)
	{
		Int_t contribchs = 0;
		Double_t meanX = 0;
		for (Int_t iCh=1; iCh<=nChs; iCh++)
		{
			double val = h->GetBinContent(xBin,iCh);
			if( (RMS_vect) && (abs(val)<(fSigmaThr*RMS_vect->at(iCh-1))) )
			{
				meanX += val;
				contribchs++;
			}
		}
		
		if(contribchs>0) meanX /= contribchs;
		h_corr->SetBinContent(xBin, h_corr->GetBinContent(xBin)+meanX);
		
		//Now apply the correction to the mean
		for (Int_t iCh=1; iCh<=nChs; iCh++)
		{
			h->SetBinContent(xBin, iCh, h->GetBinContent(xBin,iCh)-meanX );
		}
	}
	
	
	if(iterations>0)
	{//Before reapplying the CM correction recalculate baselines and RMS values for each channel
		if(bslncorr) BaselineCorr(h, RMS_vect);
		h_corr = CMrej(h, RMS_vect, iterations);
	}
	
	return h_corr;	
}


void RSTPC_Analyser::PrintRMSvalues()
{
	if(fColRMS && (fColRMS->size()>0))
	{
		cout << "\nCollection wires. Noise RMS values:" << endl;
		for(unsigned iWr=0; iWr<fColRMS->size(); iWr++)
		{
			cout << "\tWire " << iWr << ": " << fColRMS->at(iWr) << endl;
		}
	}
	
	if(fIndRMS && (fIndRMS->size()>0))
	{
		cout << "\nInduction wires. Noise RMS values:" << endl;
		for(unsigned iWr=0; iWr<fIndRMS->size(); iWr++)
		{
			cout << "\tWire " << iWr << ": " << fIndRMS->at(iWr) << endl;
		}
	}
	cout << endl;
}


void RSTPC_Analyser::Display(int Event0, int Event1, Bool_t calib)
{
	stringstream sstr;
	
	
	sstr.str(""); sstr << fOutDir << "Event_Display/Run_" << setfill('0') << setw(9) << fRunNumber;
	if( fPrintFlag && access( sstr.str().c_str(), R_OK|W_OK|X_OK )  )
	{
		string outdir = sstr.str();
		cout << "Making output directory <" << outdir << ">." << endl;
		sstr.str(""); sstr << "mkdir -p " << outdir;
		gSystem->Exec( sstr.str().c_str() );
		if( access( outdir.c_str(), R_OK|W_OK|X_OK ) )
		{
			cout << "\nERROR: Output directory <" << outdir << "> cannot be created or accessed! Switching off the flag for canvas printing." << endl;
			fPrintFlag = false;
		}
	}
	
	
	if(!fCanv) 
	{
		if(calib) scale();
		//if(SCALE==1)   scale();
		//gStyle->SetPalette(1);
		if(!fCanv) fCanv = new TCanvas("fCanv","RSTPC event display",1020,640);
		fCanv->SetFillColor(kWhite);
		fCanv->Divide(1,2);
	}
	
	
	if(hC) { delete hC; hC=NULL;}
	if(hI) { delete hI; hI=NULL;}
	
	if(Event0<0) Event0=0;
	if(Event1<Event0) Event1=Event0;
	
	for(Int_t iEv=Event0; iEv<=Event1; iEv++)
	{
		cout << "Processing event " << iEv << "..." << endl;
		
		if(!fInfile)
		{
			cout << "\nERROR: Pointer to fInfile: " << fInfile << ". Abort execution!\n" << endl;
			return;
		}
		/*
		sstr.str(""); sstr << "Ind_" << iEv;
		if(hI0){delete hI0; hI0 = NULL;}
		hI0 = (TH2D*)fInfile->Get( sstr.str().c_str() );
		if(!hI0) continue;
		
		sstr.str(""); sstr << "Col_" << iEv;
		if(hC0){delete hC0; hC0 = NULL;}
		hC0 = (TH2D*)fInfile->Get( sstr.str().c_str() );
		if(!hC0) continue;
		*/
		
		LoadEvent(iEv);
		if( !(hC0 && hI0) ) continue;
		
		//hI0->GetXaxis()->SetTitle("tZ, mcs");
		//hC0->GetXaxis()->SetTitle("tZ, mcs");
		
		sstr.str(""); sstr << "RSTPC Induction view, Event " << iEv << "; Time sample; X Wire index";
		hI0->SetTitle( sstr.str().c_str() );
		
		sstr.str(""); sstr << "RSTPC Collection view, Event " << iEv << "; Time sample; Y Wire index";
		hC0->SetTitle( sstr.str().c_str() );
		
		ApplyWiresMaps(hC0, hI0);
		
		//Find the baselines again and fill the RMS of each channel
		BaselineCorr(hC0, fColRMS);
		BaselineCorr(hI0, fIndRMS);
		
		//Before adding the histograms together clean for their common mode bg
		if(fCMrej)
		{
			CMrej(hC0, fColRMS, fCMrejIter);
			CMrej(hI0, fIndRMS, fCMrejIter);
			
			BaselineCorr(hC0, fColRMS);
			BaselineCorr(hI0, fIndRMS);
		}
		
		
		
		
		if(!hC) { hC=(TH2D*)(hC0->Clone()); hC->Reset();}
		hC->Add(hC0);
		
		if(!hI) { hI=(TH2D*)(hI0->Clone()); hI->Reset();}
		hI->Add(hI0);
		
		
		//sstr.str(""); sstr << "RSTPC Induction view, Run " << fRunNumber << ", Event " << iEv;
		//hI->SetTitle( sstr.str().c_str() );
		
		//sstr.str(""); sstr << "RSTPC Collection view, Run " << fRunNumber << ", Event " << iEv;
		//hC->SetTitle( sstr.str().c_str() );
	}
	
	if(Event1>Event0)
	{
		sstr.str(""); sstr << "RSTPC Induction view, Events " << Event0 << " - " << Event1 << "; Time sample; X Wire index";
		hI->SetTitle( sstr.str().c_str() );
		
		sstr.str(""); sstr << "RSTPC Collection view, Events " << Event0 << " - " << Event1 << "; Time sample; Y Wire index";
		hC->SetTitle( sstr.str().c_str() );
	}
	
	
	if(scaled)
	{
		hI->GetXaxis()->Set(4096,xbins);
		hC->GetXaxis()->Set(4096,xbins);
		hI->GetXaxis()->SetTitle("Z, mcs");
		hC->GetXaxis()->SetTitle("Z, mcs");
	
		hI->GetYaxis()->Set(32,ybins);
		hC->GetYaxis()->Set(32,ybins);
		hI->GetYaxis()->SetTitle("X, cm");
		hC->GetYaxis()->SetTitle("Y, cm");
		hI->Scale(adc_to_mV);
		hC->Scale(adc_to_mV);
		hI->GetZaxis()->SetTitle("mV");
		hC->GetZaxis()->SetTitle("mV");
	}

	if(scaled)
	{
		hI->SetMinimum(0*sqrt(Event1-Event0+1)); hI->SetMaximum(adc_to_mV*500*sqrt(Event1-Event0+1));
		hC->SetMinimum(0*sqrt(Event1-Event0+1)); hC->SetMaximum(adc_to_mV*500*sqrt(Event1-Event0+1));
	}
	else
	{
		//hI->SetMinimum(0*sqrt(Event1-Event0+1)); hI->SetMaximum(500*sqrt(Event1-Event0+1));
		//hC->SetMinimum(0*sqrt(Event1-Event0+1)); hC->SetMaximum(500*sqrt(Event1-Event0+1));
	}

	//hC->GetXaxis()->SetRangeUser(minbin,maxbin);
	//hI->GetXaxis()->SetRangeUser(minbin,maxbin);
	
	
	fCanv->cd(1);
	hC->Draw("colz");
	
	fCanv->cd(2);
	hI->Draw("colz");
	
	//c->SetFillColor(kWhite);
	//c->Draw();

	sstr.str(""); sstr << setfill('0') << setw(9);
	sstr << fOutDir << "Event_Display/Run_" << setfill('0') << setw(9) << fRunNumber << "/R" << fRunNumber << "E" << setfill('0') << setw(5);
	sstr << setfill('0') << setw(5);
	sstr << Event0;
	
	if(!fPrintFlag) return;
	
	if(Event0==Event1)
	{
		sstr << ".png";
	}
	else
	{
		sstr << "-" << setfill('0') << setw(5) << Event1 << ".png";
	}
	
	cout << "Saving canvas as <" << sstr.str() << ">" << endl;
	fCanv->SaveAs( sstr.str().c_str() );
	
}


void RSTPC_Analyser::preamp(TH2I * h) //droop compensation is here
{
	Double_t Igl=0;
	Double_t Val=0;
	for (int y=0; y<h->GetNbinsY(); y++)
	{
		Igl=0;
		//   bgx=h->ProjectionX("bgx",y,y);
		for (int x=0; x<h->GetNbinsX(); x++) 
		{
			Val= h->GetBinContent(x,y);
			h->SetBinContent(x,y,Val+Igl*0.7e-4);
			Igl=Igl+Val;
		}
	}
}


void RSTPC_Analyser::RunDisplay(ULong_t rnumber, Int_t evt0, Int_t evt1, Float_t minbin, Float_t maxbin, Bool_t calib)
{
	OpenRun(rnumber);
	//if(evt1>evt0) mdisplay(evt0,evt1, minbin,maxbin,calib);
	//else display(evt0,evt1, minbin,maxbin,calib);
	Display(evt0, evt1);
}


#endif /* RSTCP_ANALYSER_CC */