#ifndef RSTPCGLOBALS_HH
#define RSTPCGLOBALS_HH

#include "TRint.h"
#include "TROOT.h"

#include <string>
#include <vector>

using namespace std;


//Global functions for general use
int pow(int base, int exponent);
bool fexists(const char *filename);
void usage(){return;};

//Use this class instead of extern declared variables

//Better if this class is a singleton class (only one instance for the entire run)
class RSTPC_Options
{//Singletorn class
private:

	static RSTPC_Options* globalThis;

	bool parseArgs(const vector<string> &args);
	
	RSTPC_Options();
	RSTPC_Options(int argc, char** argv);

	~RSTPC_Options(){;};

	string fDataDir, fXmlProcFile, fOutDir, fManOutFile;
	
	Int_t fRunNumber;
	
	bool fRunSetFlag, fDataDirFlag, fOutDirModeFlag, fXmlProcFlag, fT1ProcFlag, fT2ProcFlag, fManOutFileFlag;
	
		
public:
	static RSTPC_Options* GetInstance(int argc, char** argv);
	static RSTPC_Options* GetInstance();
	
	void SetRunNumber(Int_t RunNumber){ fRunSetFlag=(RunNumber>=0); fRunNumber=RunNumber;};
	void SetDataDir(string dirpath){fDataDirFlag=true; fDataDir=dirpath;};
	void SetOutDir(string dirpath){fOutDirModeFlag=true; fOutDir=dirpath;};
	void SetOutFile(string outfilename);
	
	Int_t GetRun() const {return fRunNumber;};
	string GetDataDir(){return fDataDir;};
	//string GetXmlProcFile() const {return fXmlProcFile;};
	string GetOutDir() const {return fOutDir;};
	string GetOutFile() const {return fManOutFile;};
	
	Bool_t IsDataDirSet() const {return fDataDirFlag;};
	bool IsOutFileMode() const {return fManOutFileFlag;};
	bool IsRunNumberSet() const {return fRunSetFlag;};
	bool IsT1proc() const {return fT1ProcFlag;};
	bool IsT2proc() const {return fT2ProcFlag;};
	bool IsOutDirSet() const {return fOutDirModeFlag;};
	//bool IsEvList() const {return fEvListFlag;};

};



#endif /* RSTPCGLOBALS_HH */