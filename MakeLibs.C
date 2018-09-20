void MakeLibs()
{
	gSystem->SetBuildDir(gSystem->GetBuildDir(), true);
	
	gInterpreter->AddIncludePath("/home/francescop/ArCube/analysis/code");
	gInterpreter->AddIncludePath("/home/francescop/analysis/include");
	
	//gSystem->Unload("RSTPC_Globals");
	//gSystem->Unload("RSTPC_Analyser");
	//gSystem->Unload("MppcTreeWrapper");
	//gSystem->Unload("RSTPC_RunProcessor");
	//gSystem->Unload("RSTPC_T1wrapper");
	
	gSystem->CompileMacro("RSTPC_Globals.cc", "gkf", "RSTPC_Globals");
	gSystem->CompileMacro("RSTPC_Analyser.cc", "gkf", "RSTPC_Analyser");
	gSystem->CompileMacro("MppcTreeWrapper.cc", "gkf", "MppcTreeWrapper");
	gSystem->CompileMacro("RSTPC_T1wrapper.cc", "gkf", "RSTPC_T1wrapper");
	gSystem->CompileMacro("RSTPC_Hits.cc", "gkf", "RSTPC_Hits");
	gSystem->CompileMacro("RSTPC_RunProcessor.cc", "gkf", "RSTPC_RunProcessor");
}