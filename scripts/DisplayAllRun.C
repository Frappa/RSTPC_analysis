#include "RSTPC_Globals.hh"
#include "RSTPC_Analyser.hh"
#include "HistoManipulators.hh"
#include "DigitalFilters.hh"


RSTPC_Analyser *an = NULL;


void DisplayAllRun(Int_t RunNumber, Bool_t CMrej=true, Bool_t bAllEvs=true)
{
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell");
	
	if(!an)
	{
		an = new RSTPC_Analyser;
		an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
		an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
	}
	
	an->OpenRun(RunNumber);
//an->OpenRun(2026); //Peaking time 2 usec, 15 kV
//an->OpenRun(2027); //Peaking time 2 usec, 3 kV
//an->OpenRun(2028); //Peaking time 2 usec, 3 kV
//an->OpenRun(2031); //Peaking time 1 usec, 15 kV
//an->OpenRun(2032); //Peaking time 1 usec, 20 kV
//an->OpenRun(2033); //Peaking time 1 usec, 30 kV
	
	an->SetPrintFlag(false);
	if(CMrej) an->Set_CMnoiseRej(true, 5);
	
	an->SetBaselineROI(3000, 4000);
	an->SetSigmaThr(3.0);
	
	if(bAllEvs)
	{
		Int_t nEvs = an->fTrigTree->GetEntries();
		an->Display(0, nEvs-1);
		return;
	}
	
	an->Display(0, 0);
	return;
}
