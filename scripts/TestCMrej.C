#include "RSTPC_Analyser.hh"
#include "HistoManipulators.hh"
#include "DigitalFilters.hh"


RSTPC_Analyser *an = NULL;

void TestCMrej(Int_t RunNum, Int_t EvNnum=28){

RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell/");
//R__LOAD_LIBRARY(RSTPC_Analyser.so)
//R__LOAD_LIBRARY(HistoManipulators.so)
//R__LOAD_LIBRARY(DigitalFilters.so)


	//gROOT->ProcessLine(".L display_vRS.C+");
	//gSystem->Load("RSTPC_Analyser.so");
	
	//gROOT->ProcessLine("RSTPC_Analyser *ed = new RSTPC_Analyser;");
	
	if(!an){
		an = new RSTPC_Analyser;
		an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
		an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
		an->OpenRun(RunNum);
		an->Set_CMnoiseRej(false);
		//an->Set_CMnoiseRej(true, 2);
		an->SetPrintFlag(false);
		an->SetBaselineROI(2000, 3000);
		an->SetSigmaThr(3.0);
	}
	
	//
	
	
	an->LoadEvent(EvNnum);
	an->ApplyWiresMaps(an->hC0,an->hI0);
	
	
	TCanvas* c0 = new TCanvas("c0","Before any correction", 1020,640);
	c0->Divide(1,2);
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH2D* hC0_0 = (TH2D*)an->hC0->Clone("hC0_0");
	hC0_0->GetXaxis()->SetTitle("Drift time [samples]");
	hC0_0->GetYaxis()->SetTitle("Wire number");
	
	TH2D* hI0_0 = (TH2D*)an->hI0->Clone("hI0_0");
	hI0_0->GetXaxis()->SetTitle("Drift time [samples]");
	hI0_0->GetYaxis()->SetTitle("Wire number");
	
	
	c0->cd(1);
	hC0_0->Draw("colz");
	
	c0->cd(2);
	hI0_0->Draw("colz");
	
	
	
	//First round of corrections
	TCanvas* c1 = new TCanvas("c1","After first round of bsln shift and corrections", 1020,640);
	c1->Divide(1,2);
	
	an->CMrej(an->hC0, an->fColRMS, false);
	an->CMrej(an->hI0, an->fIndRMS, false);
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH2D* hC0_1 = (TH2D*)an->hC0->Clone("hC0_1");
	TH2D* hI0_1 = (TH2D*)an->hI0->Clone("hI0_1");
	
	c1->cd(1);
	hC0_1->Draw("colz");
	
	c1->cd(2);
	hI0_1->Draw("colz");
	
	
	
	//Second round of corrections
	TCanvas* c2 = new TCanvas("c2","After the second round of bsln shift and corrections", 1020,640);
	c2->Divide(1,2);
	
	an->CMrej(an->hC0, an->fColRMS, false);
	an->CMrej(an->hI0, an->fIndRMS, false);
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH2D* hC0_2 = (TH2D*)an->hC0->Clone("hC0_2");
	TH2D* hI0_2 = (TH2D*)an->hI0->Clone("hI0_2");
	
	c2->cd(1);
	hC0_2->Draw("colz");
	
	c2->cd(2);
	hI0_2->Draw("colz");
	
}


void RMSplot(Int_t RunNum, Int_t EvNnum)
{
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell/");
	
	if(!an){
		an = new RSTPC_Analyser;
		an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
		an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
		an->OpenRun(RunNum);
		an->SetPrintFlag(false);
		an->SetSigmaThr(3.0);
		an->SetBaselineROI(2000, 3000);
	}
	
	
	an->LoadEvent(EvNnum);
	an->ApplyWiresMaps(an->hC0,an->hI0);
	
	
	
	//Before corrections
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH1D* hC_RMS_0 = VectorToHisto(an->fColRMS, -0.5, 31.5);
	if(hC_RMS_0) hC_RMS_0->SetName("hC_RMS_0"); hC_RMS_0->SetLineColor(kBlue); hC_RMS_0->SetTitle("Collection wires;Wire number;RMS amplitude [AU]");
	TH1D* hI_RMS_0 = VectorToHisto(an->fIndRMS, -0.5, 31.5);
	if(hI_RMS_0) hI_RMS_0->SetName("hI_RMS_0"); hI_RMS_0->SetLineColor(kBlue); hI_RMS_0->SetTitle("Induction wires;Wire number;RMS amplitude [AU]");
	
	
	
	//First round of CM rejection
	
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH1D* hC_RMS_1 = VectorToHisto(an->fColRMS, -0.5, 31.5);
	if(hC_RMS_1) hC_RMS_1->SetName("hC_RMS_1"); hC_RMS_1->SetLineColor(kMagenta);
	TH1D* hI_RMS_1 = VectorToHisto(an->fIndRMS, -0.5, 31.5);
	if(hI_RMS_1) hI_RMS_1->SetName("hI_RMS_1"); hI_RMS_1->SetLineColor(kMagenta);
	
	
	
	//Second round of CM rejection
	TH1D* hC0_2_corr = an->CMrej(an->hC0, an->fColRMS, false);
	if(hC0_2_corr) hC0_2_corr = (TH1D*)hC0_2_corr->Clone("hC0_2_corr");
	
	TH1D* hI0_2_corr = an->CMrej(an->hI0, an->fIndRMS, false);
	if(hI0_2_corr) hI0_2_corr = (TH1D*)hI0_2_corr->Clone("hI0_2_corr");
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH1D* hC_RMS_2 = VectorToHisto(an->fColRMS, -0.5, 31.5);
	if(hC_RMS_2) hC_RMS_2->SetName("hC_RMS_2"); hC_RMS_2->SetLineColor(kRed);
	TH1D* hI_RMS_2 = VectorToHisto(an->fIndRMS, -0.5, 31.5);
	if(hI_RMS_2) hI_RMS_2->SetName("hI_RMS_2"); hI_RMS_2->SetLineColor(kRed);
	
	
	//Third round of CM rejection
	TH1D* hC0_3_corr = an->CMrej(an->hC0, an->fColRMS, false);
	if(hC0_3_corr) hC0_3_corr = (TH1D*)hC0_3_corr->Clone("hC0_3_corr");
	
	TH1D* hI0_3_corr = an->CMrej(an->hI0, an->fIndRMS, false);
	if(hI0_3_corr) hI0_3_corr = (TH1D*)hI0_3_corr->Clone("hI0_3_corr");
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH1D* hC_RMS_3 = VectorToHisto(an->fColRMS, -0.5, 31.5);
	if(hC_RMS_3) hC_RMS_3->SetName("hC_RMS_3"); hC_RMS_3->SetLineColor(kOrange);
	TH1D* hI_RMS_3 = VectorToHisto(an->fIndRMS, -0.5, 31.5);
	if(hI_RMS_3) hI_RMS_3->SetName("hI_RMS_3"); hI_RMS_3->SetLineColor(kOrange);
	
	
	//Fourth round of CM rejection
	TH1D* hC0_4_corr = an->CMrej(an->hC0, an->fColRMS, false);
	if(hC0_4_corr) hC0_4_corr = (TH1D*)hC0_4_corr->Clone("hC0_4corr");
	
	TH1D* hI0_4_corr = an->CMrej(an->hI0, an->fIndRMS, false);
	if(hI0_4_corr) hI0_4_corr = (TH1D*)hI0_4_corr->Clone("hI0_4_corr");
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH1D* hC_RMS_4 = VectorToHisto(an->fColRMS, -0.5, 31.5);
	if(hC_RMS_4) hC_RMS_4->SetName("hC_RMS_4"); hC_RMS_4->SetLineColor(kGreen+1);
	TH1D* hI_RMS_4 = VectorToHisto(an->fIndRMS, -0.5, 31.5);
	if(hI_RMS_4) hI_RMS_4->SetName("hI_RMS_4"); hI_RMS_4->SetLineColor(kGreen+1);
	
	
	//Fifth round of CM rejection
	TH1D* hC0_5_corr = an->CMrej(an->hC0, an->fColRMS, false);
	if(hC0_5_corr) hC0_5_corr = (TH1D*)hC0_5_corr->Clone("hC0_5corr");
	
	TH1D* hI0_5_corr = an->CMrej(an->hI0, an->fIndRMS, false);
	if(hI0_5_corr) hI0_5_corr = (TH1D*)hI0_5_corr->Clone("hI0_5_corr");
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	TH1D* hC_RMS_5 = VectorToHisto(an->fColRMS, -0.5, 31.5);
	if(hC_RMS_5) hC_RMS_5->SetName("hC_RMS_5"); hC_RMS_5->SetLineColor(kCyan+1);
	TH1D* hI_RMS_5 = VectorToHisto(an->fIndRMS, -0.5, 31.5);
	if(hI_RMS_5) hI_RMS_5->SetName("hI_RMS_5"); hI_RMS_5->SetLineColor(kCyan+1);
	
	
	TCanvas *c1 = (TCanvas*)gROOT->FindObject("canvRMSplot1");
	if(!c1)
	{
		c1 = new TCanvas("canvRMSplot1","RMS evolution",800,600);
		c1->Divide(1,2);
	}
	
	c1->cd(1);
	hC_RMS_0->Draw();
	hC_RMS_1->Draw("same");
	hC_RMS_2->Draw("same");
	hC_RMS_3->Draw("same");
	hC_RMS_4->Draw("same");
	hC_RMS_5->Draw("same");
	
	TLegend *leg = new TLegend(0.7,0.9,0.9,0.65);
	leg->AddEntry(hC_RMS_0, "Before CM rejection");
	leg->AddEntry(hC_RMS_1, "CM rejection: 1 iteration", "L");
	leg->AddEntry(hC_RMS_2, "CM rejection: 2 iterations", "L");
	leg->AddEntry(hC_RMS_3, "CM rejection: 3 iterations", "L");
	leg->AddEntry(hC_RMS_4, "CM rejection: 4 iterations", "L");
	leg->AddEntry(hC_RMS_5, "CM rejection: 5 iterations", "L");
	leg->Draw();
	
	c1->cd(2);
	hI_RMS_0->Draw();
	hI_RMS_1->Draw("same");
	hI_RMS_2->Draw("same");
	hI_RMS_3->Draw("same");
	hI_RMS_4->Draw("same");
	hI_RMS_5->Draw("same");
}


void CorrPlot(Int_t RunNum, Int_t EvNnum)
{
	RSTPC_Options::GetInstance()->SetDataDir("/home/francescop/data/ResistiveShell/");
	
	if(!an){
		an = new RSTPC_Analyser;
		an->LoadCollMap("../ResistiveShellTPC/CollWireMap.txt");
		an->LoadIndcMap("../ResistiveShellTPC/InducWireMap.txt");
		an->OpenRun(RunNum);
		an->Set_CMnoiseRej(true, 1);
		an->SetPrintFlag(false);
		an->SetBaselineROI(2000, 3000);
		an->SetSigmaThr(3.0);
	}
	
	
	an->LoadEvent(EvNnum);
	an->ApplyWiresMaps(an->hC0,an->hI0);
	
	
	
	//Before corrections
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	
	//First round of corrections
	TH1D* hC0_1_corr = an->CMrej(an->hC0, an->fColRMS, false);
	if(hC0_1_corr) hC0_1_corr = (TH1D*)hC0_1_corr->Clone("hC0_1_corr"); hC0_1_corr->SetLineColor(kBlue);
	
	TH1D* hI0_1_corr = an->CMrej(an->hI0, an->fIndRMS, false);
	if(hI0_1_corr) hI0_1_corr = (TH1D*)hI0_1_corr->Clone("hI0_1_corr"); hI0_1_corr->SetLineColor(kBlue);
	
	TH1D* hC0_1_corr_ampl = MakeAmplitudeHisto("hC0_1_corr_ampl", hC0_1_corr, 100, 1.10, kBlue);
	TH1D* hI0_1_corr_ampl = MakeAmplitudeHisto("hI0_1_corr_ampl", hI0_1_corr, 100, 1.10, kBlue);
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	
	//Second round of corrections
	TH1D* hC0_2_corr = an->CMrej(an->hC0, an->fColRMS, false);
	if(hC0_2_corr) hC0_2_corr = (TH1D*)hC0_2_corr->Clone("hC0_2_corr"); hC0_2_corr->SetLineColor(kGreen+1);
	
	TH1D* hI0_2_corr = an->CMrej(an->hI0, an->fIndRMS, false);
	if(hI0_2_corr) hI0_2_corr = (TH1D*)hI0_2_corr->Clone("hI0_2_corr"); hI0_2_corr->SetLineColor(kGreen+1);
	
	TH1D* hC0_2_corr_ampl = MakeAmplitudeHisto("hC0_2_corr_ampl", hC0_2_corr, 100, 1.10, kGreen+1);
	TH1D* hI0_2_corr_ampl = MakeAmplitudeHisto("hI0_2_corr_ampl", hI0_2_corr, 100, 1.10, kGreen+1);
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	
	//Third round of corrections
	TH1D* hC0_3_corr = an->CMrej(an->hC0, an->fColRMS, false);
	if(hC0_3_corr) hC0_3_corr = (TH1D*)hC0_3_corr->Clone("hC0_3_corr"); hC0_3_corr->SetLineColor(kRed);
	
	TH1D* hI0_3_corr = an->CMrej(an->hI0, an->fIndRMS, false);
	if(hI0_3_corr) hI0_3_corr = (TH1D*)hI0_3_corr->Clone("hI0_3_corr"); hI0_3_corr->SetLineColor(kRed);
	
	TH1D* hC0_3_corr_ampl = MakeAmplitudeHisto("hC0_3_corr_ampl", hC0_3_corr, 100, 1.10, kRed);
	TH1D* hI0_3_corr_ampl = MakeAmplitudeHisto("hI0_3_corr_ampl", hI0_3_corr, 100, 1.10, kRed);
	
	an->BaselineCorr(an->hC0, an->fColRMS);
	an->BaselineCorr(an->hI0, an->fIndRMS);
	
	
	TCanvas *c3 = (TCanvas*)gROOT->FindObject("c3");
	if(!c3)
	{
		c3 = new TCanvas("c3","Correction histograms",1280,640);
		c3->Divide(1,2);
	}
	
	c3->cd(1);
	hC0_1_corr->Draw();
	hC0_2_corr->Draw("same");
	hC0_3_corr->Draw("same");
	
	c3->cd(2);
	hI0_1_corr->Draw();
	hI0_2_corr->Draw("same");
	hI0_3_corr->Draw("same");
	
	
	TCanvas *c4 = (TCanvas*)gROOT->FindObject("c4");
	if(!c4)
	{
		c4 = new TCanvas("c4","Correction histograms",1280,640);
		c4->Divide(2,1);
	}
	
	c4->cd(1);
	hC0_1_corr_ampl->Draw();
	hC0_2_corr_ampl->Draw("same");
	hC0_3_corr_ampl->Draw("same");
	
	c4->cd(2);
	hI0_1_corr_ampl->Draw();
	hI0_2_corr_ampl->Draw("same");
	hI0_3_corr_ampl->Draw("same");
	
}