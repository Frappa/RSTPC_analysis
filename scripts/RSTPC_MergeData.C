#include "RSTPC_Globals.hh"
#include "RSTPC_RunProcessor.hh"
#include "RSTPC_Analyser.hh"
#include "MppcTreeWrapper.hh"


//#if defined(__CLING__) && defined(__ROOTCLING__)
RSTPC_RunProcessor *merger = NULL;
void RSTPC_MergeData(Int_t RunNumber)
{
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell/");
	//RSTPC_Options::GetInstance()->SetRunNumber(RunNumber);
	
	merger = new RSTPC_RunProcessor;
	
	merger->GetTpcManager()->Set_CMnoiseRej(true, 5);
	merger->GetTpcManager()->SetBaselineROI(2000, 3000);
	merger->GetTpcManager()->SetSigmaThr(3.0);
	merger->GetTpcManager()->SetPrintFlag(false);
	
	if( !merger->InitT1proc(RunNumber) ) return;
	
	merger->LoadColMap("/home/francescop/ArCube/analysis/ResistiveShellTPC/CollWireMap.txt");
	merger->LoadIndMap("/home/francescop/ArCube/analysis/ResistiveShellTPC/InducWireMap.txt");
	
	merger->DescribeT1();
	
	merger->T1process();
	
	//merger->CloseRun();
	
	//delete merger;
	
	return;
}

//#endif

/*
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell/"); RSTPC_Options::GetInstance()->SetRunNumber(2031); RSTPC_RunMerger *merger = new RSTPC_RunMerger; merger->LoadColMap("/home/francescop/ArCube/analysis/ResistiveShellTPC/CollWireMap.txt"); merger->LoadIndMap("/home/francescop/ArCube/analysis/ResistiveShellTPC/InducWireMap.txt"); merger->DescribeT1();
*/