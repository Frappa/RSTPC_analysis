#include "RSTPC_Globals.hh"
#include "RSTPC_RunProcessor.hh"
#include "RSTPC_Analyser.hh"
#include "MppcTreeWrapper.hh"


//#if defined(__CLING__) && defined(__ROOTCLING__)
RSTPC_RunProcessor *T2processor = NULL;
void RSTPC_T2proc()
{
	RSTPC_RunProcessor::SetDebug(true, 100);
	
	RSTPC_Options::GetInstance()->SetOutFile("/home/francescop/data/ResistiveShell/merged/test/RSTPC_Run000002032_Merged.root");
	RSTPC_Options::GetInstance()->SetDataDir( RSTPC_Options::GetInstance()->GetOutDir() ); //I need this instruction to avoid that the RSTPC_RunProcessor ctor crashes
	
	T2processor = new RSTPC_RunProcessor;
	
	RSTPC_RunProcessor::SetPitchSize( 52.5/31 ); //Measured
	
	RSTPC_RunProcessor::SetPeakingTime(1.0);
	
	RSTPC_RunProcessor::SetDriftVel(2.06); //At 1 kV/cm
	//RSTPC_RunProcessor::SetDriftVel(2.1); //At 1.33 kV/cm
	
	if(T2processor->InitT2proc())
	{
		T2processor->DescribeT2();
		T2processor->T2Process();
	}
	
	
	return;
}

//#endif

/*
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell/"); RSTPC_Options::GetInstance()->SetRunNumber(2031); RSTPC_RunMerger *merger = new RSTPC_RunMerger; merger->LoadColMap("/home/francescop/ArCube/analysis/ResistiveShellTPC/CollWireMap.txt"); merger->LoadIndMap("/home/francescop/ArCube/analysis/ResistiveShellTPC/InducWireMap.txt"); merger->DescribeT1();
*/