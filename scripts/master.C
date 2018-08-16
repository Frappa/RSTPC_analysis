void master() {
  gROOT->ProcessLine(".x LoadLibs.C");
  gROOT->ProcessLine(".L RSTPC_T2wrapper.cc");
  //gROOT->ProcessLine(".x StraightLines.C");

}
