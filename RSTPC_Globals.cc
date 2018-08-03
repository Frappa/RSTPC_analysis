#include "RSTPC_Globals.hh"

#include <iostream>


using namespace std;

RSTPC_Options* RSTPC_Options::globalThis = NULL;



RSTPC_Options* RSTPC_Options::GetInstance()
{
	if(!globalThis) globalThis = new RSTPC_Options();
	return globalThis;
}

RSTPC_Options* RSTPC_Options::GetInstance(int argc, char** argv)
{
	
	if(!globalThis) globalThis = new RSTPC_Options(argc, argv);
	
	return globalThis;
}

RSTPC_Options::RSTPC_Options()
{
	fRunNumber = -1;
	
	fRunSetFlag = false;
	//fXmlProcFlag = false;
	fT1ProcFlag = true; //This is always true unless a T2 process is defined
	//fT2ProcFlag = false;
	fOutDirModeFlag = false;
	fManOutFileFlag = false;
}


RSTPC_Options::RSTPC_Options(int argc, char** argv)
{
	fRunNumber = -1;
	
	fRunSetFlag = false;
	//fXmlProcFlag = false;
	fT1ProcFlag = true; //This is always true unless a T2 process is defined
	//fT2ProcFlag = false;
	fOutDirModeFlag = false;
	fManOutFileFlag = false;

	vector<string> args(argc-1);

	for(int iArg=0; iArg<(argc-1); iArg++){
		args.at(iArg) = argv[iArg+1];
	}

	if(!parseArgs(args))
	{
		usage();
		exit(1);
	}
	
	
}


bool RSTPC_Options::parseArgs(const vector<string> &args)
{
	string tmpstr;
	for(unsigned iArg=0; iArg<args.size(); iArg++)
	{
		if(args.at(iArg)==string("-r"))
		{
			tmpstr = args.at(iArg+1);
			if(tmpstr.at(0)=='-') return false;
			fRunNumber = atoi(tmpstr.c_str());
			if(fRunNumber>=0) fRunSetFlag = true;
			iArg++;
		}
		else if(args.at(iArg)==string("-f"))
		{
			tmpstr = args.at(iArg+1);
			//cout << "\nOutput file (manual): <" << outfilename << ">" << endl;
			if(tmpstr.at(0)=='-') return false;
			SetOutFile(tmpstr);
			iArg++;
		}
		else if(args.at(iArg)==string("-O"))
		{
			if(!fManOutFileFlag)
			{
				tmpstr = args.at(iArg+1);
				if(tmpstr.at(0)=='-') return false;
				SetOutDir(tmpstr);
				//cout << "\nOutput main directory (manual): <" << gXuOptions.outdir << ">" << endl;
			}
			iArg++;
		}
		else
		{
			cerr << "\nERROR --> RSTPC_Options::parseArgs: Unknown option \"" << args.at(iArg) << "\"\n" << endl;
			return false;
		}
		
		/*
		switch(args.at(iArg))
		{
			case string("-p"):
			fXmlProcFile = args(iArg+1);
			//cout << "\nSingle file to process: <" << rootfilename << ">" << endl;
			if(fXmlProcFile.at(0)=='-') return false;
			iArg++;
			fXmlProcFlag = true;
			break;
			
			case string("-T1"):
			iArg++;
			fT1ProcFlag = true; //Only T1 processing
			break;
			
			case string("-T2"): //Only T2 processing
			iArg++;
			fT2ProcFlag = true;
			break;
		}
		*/
	}
	
	return true;
}


void RSTPC_Options::SetOutFile(string outfilename)
{
	fManOutFileFlag = true;
	
	Bool_t path = false;
	
	if(outfilename.find('/')!=string::npos)
	{
		path = true;
	}
	
	if( !path )
	{
		fManOutFile=outfilename;
	}
	else
	{
		size_t strsize = outfilename.size();
		size_t find1 = outfilename.find_last_of("/");
		string dir = outfilename.substr( 0, find1);
		SetOutDir(dir);
		fManOutFile = outfilename.substr( find1+1, strsize-find1);
	}
}