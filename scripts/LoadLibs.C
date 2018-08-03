{
gInterpreter->AddIncludePath("/home/francescop/ArCube/analysis/code");
gInterpreter->AddIncludePath("/home/francescop/analysis/code");

gSystem->AddDynamicPath("/home/francescop/ACLiC_build");

gSystem->Load("DigitalFilters.so");
gSystem->Load("HistoManipulators.so");
gSystem->Load("RSTPC_Globals.so");
gSystem->Load("RSTPC_Analyser.so");
gSystem->Load("MppcTreeWrapper.so");
gSystem->Load("RSTPC_T1wrapper.so");
gSystem->Load("RSTPC_Hits.so");
gSystem->Load("RSTPC_RunProcessor.so");
}
