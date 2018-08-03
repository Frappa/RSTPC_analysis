#ifndef RSTPC_ANALYSER_HH
#define RSTPC_ANALYSER_HH

#include "RSTPC_Globals.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"

#include <string>
#include <sstream>
#include <fstream>
#include <map>

using namespace std;

//This class is specifically designed to analyse the data acquired for the Resistive Shell TPC, and hardly can be adapted to other experiments.
//The data is found in the rootfile of the run as 2 histograms for each event.
//One histogram is for the "collection" wires (32 in total) and the other for the "induction" wires (32 in total).
class RSTPC_Analyser
{
public:
	
	TCanvas *fCanv = NULL;
	
	TH2D *hI, *hC; //These are the only non-temporary histograms
	
	//All these histograms are overwritten at each event
	TH2D *hI0, *hC0, *hII;
	TH1F *hE, *hU, *hR;

	Int_t fRunNumber;
	TFile* fInfile;
	FileStat_t stat;
	Int_t scaled;
	UChar_t pattern;
	TTree* fTrigTree;
	Float_t adc_to_mV; //mV per ADC
	//#define SCALE 1
	TColor vvv;
	
	string fDataDir;
	string fOutDir;
	
	Int_t fXlowBin, fXupBin;
	bool fCMrej;
	Int_t fCMrejIter;
	
	Bool_t fPrintFlag;
	Bool_t fLoadedRun; //It is true when the rootfile of the run is found and opened
	
	vector<Double_t> *fColRelGains, *fIndRelGains;
	
	map<Int_t, Int_t> *fCollectionMap, *fInductionMap;//The pairs are: (wire, channel)
	
	vector<Double_t> *fColRMS, *fIndRMS; //They store the noise level for each channel after the common mode suppression
	
	Double_t fSigmaThr; //Threshold above which the fluctuation is considered a signal
	
	ULong64_t fEventTime;
	
	RSTPC_Analyser();
	
	virtual ~RSTPC_Analyser();
	
	
	Bool_t IsInit()
	{
		return fLoadedRun && fInfile && fInfile->IsOpen() && fTrigTree;
	}
	
	void Set_CMnoiseRej(bool flag, Int_t iterations=1)
	{
		fCMrej = flag;
		if(flag) fCMrejIter = iterations;
	};
	
	void SetPrintFlag(bool flag)
	{
		fPrintFlag = flag;
	};
	
	void SetBaselineROI(Int_t _XlowBin, Int_t _XupBin)
	{
		//Set the ROI region (in samples units) where the baseline for each channel is calculated
		fXlowBin = _XlowBin;
		fXupBin = _XupBin;
	};
	
	void SetSigmaThr(Double_t thr){fSigmaThr = thr;};
	
	virtual Bool_t OpenRun(Int_t RunNumber);
	void LoadCollMap(string filename);
	void LoadIndcMap(string filename);
	void ApplyWiresMaps(TH2 *hc, TH2 *hi);
	void BaselineCorr(TH2 * h, Int_t _XlowBin, Int_t _XupBin, vector<Double_t>* RMS_vect=NULL);
	void BaselineCorr(TH2 * h, vector<Double_t>* RMS_vect=NULL)
	{
		BaselineCorr(h, fXlowBin, fXupBin, RMS_vect);
	};
	
	void Display(int Event0, int Event1, Bool_t calib=false);
	void RunDisplay(ULong_t rnumber, Int_t evt0, Int_t evt1=0, Float_t minbin=0, Float_t maxbin=5., Bool_t calib=false);
	
	void preamp(TH2I* h); //droop compensation is here
	void scale();
	
	//The common mode noise rejection methods return the correction histogram "bgx"
	TH1D* CMrej(TH2* h, vector<Double_t>* RMS_vect, Int_t iterations, Bool_t bslncorr=true);
	TH1D* CMrej(TH2* h, Int_t iterations, Bool_t bslncorr=true)
	{
		return CMrej(h, NULL, iterations, bslncorr);
	};
	TH1D* CMrej(TH2 * h, Bool_t bslncorr)
	{
		return CMrej(h, NULL, 1, bslncorr);
	};
	TH1D* CMrej(TH2 * h, vector<Double_t>* RMS_vect, Bool_t bslncorr)
	{
		return CMrej(h, RMS_vect, 1, bslncorr);
	};
	
	void LoadEvent(Int_t evNumber);
	
	void FindGaines(){;};
	
	void PrintRMSvalues();
};


//#if defined(__CLING__)
//#include "RSTPC_Analyser.cc"
//#endif

#endif /* RSTPC_ANALYSER_HH */