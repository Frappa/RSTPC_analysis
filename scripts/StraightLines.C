#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TObjArray.h"
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <algorithm>

//#include "RSTPC_T1wrapper.hh"
//#include "RSTPC_T2wrapper.hh"
//#include "TestPulses.C"
//#include "../RSTPC_T1wrapper.hh" // for the RSTPC_T1wrapper class; DOESN'T WORK
//#include "../RSTPC_Hits.hh" // for the RSTPC_Hit and RSTPC_Pulse classes

//#include "PCA.C"



// ***************************************************
void StraightLines() {
// ***************************************************

    double e_drift_velocity = 0.1648; // cm/us @ 0.5 kV/cm

    //gROOT->ProcessLine(".x LoadLibs.C");
    //gROOT->ProcessLine(".L RSTPC_T2wrapper.cc");



    // Define data files
    // ===================================================================
    TFile * tfile_run_2032 = new TFile("/home/francescop/data/ResistiveShell/merged/test/RSTPC_Run000002032_Merged.root");
									 // /home/rberner   /data/ResistiveShell/merged/test/RSTPC_Run000002032_Merged.root");


    // Access the T1 data
    // ===================================================================
    // Get T1 tree and link variables
    TTree * T1 = (TTree*)tfile_run_2032->Get("T1");

	Int_t       TpcEvent;
    Double_t    TpcTime;
    Double_t    RmsColWires[32];
    Double_t    RmsIndWires[32];
    Int_t       FebEvent;
    Double_t    FebTime;
    UShort_t    FebTopAmp[3];
    Double_t    FebTopTotAmp;
    UShort_t    FebBotAmp[3];
    Double_t    FebBotTotAmp;
    Double_t    FebTotAmp;

    T1->SetBranchAddress("TpcEvent", &TpcEvent);
    T1->SetBranchAddress("TpcTime", &TpcTime);
    T1->SetBranchAddress("RmsColWires", &RmsColWires);
    T1->SetBranchAddress("RmsIndWires", &RmsIndWires);
    T1->SetBranchAddress("FebEvent", &FebEvent);
    T1->SetBranchAddress("FebTime", &FebTime);
    T1->SetBranchAddress("FebTopAmp", &FebTopAmp);
    T1->SetBranchAddress("FebTopTotAmp", &FebTopTotAmp);
    T1->SetBranchAddress("FebBotAmp", &FebBotAmp);
    T1->SetBranchAddress("FebBotTotAmp", &FebBotTotAmp);
    T1->SetBranchAddress("FebTotAmp", &FebTotAmp);

    std::cout << " ---------------------------------- " << std::endl;
    std::cout << " Processing T1 with " << T1->GetEntries() << " entries: " << std::endl;
    std::cout << " ---------------------------------- " << std::endl;
    //std::cout << " Tpc \t Tpc \t\t Rms \t Rms \t\t Feb \t Feb \t\t FebT \t FebT \t FebB \t FebB \t Feb " << std::endl;
    //std::cout << " Evnt: \t Time: \t\t ColW: \t IndW: \t\t Evnt: \t Time: \t\t Amp: \t TotAmp: Amp: \t TotAmp: TotAmp: " << std::endl << std::endl;
    for(int entry=2; entry<3; entry++) {
      	T1->GetEntry(entry);
      	/*std::cout 	<< " " << TpcEvent	<< " \t " << TpcTime 		<< " \t " 		<< RmsColWires[0] 	<< " \t "
                	<< RmsIndWires[5] 	<< " \t " << FebEvent 		<< " \t " 		<< FebTime 			<< " \t "
                	<< FebTopAmp[2] 	<< " \t " << FebTopTotAmp 	<< " \t " 		<< FebBotAmp[0] 	<< " \t "
                	<< FebBotTotAmp 	<< " \t " << FebTotAmp 		<< std::endl;*/
    	//T1->Show(entry);
    }
    std::cout << std::endl;


    // Access the T2 data
    // ===================================================================
    // Get T2 tree and link variables
    TTree * T2 = (TTree*)tfile_run_2032->Get("T2");
    TTree * testtree = (TTree*)tfile_run_2032->Get("T2");

    TClonesArray * tclonesarray_T2_ColPulses = new TClonesArray("RSTPC_Pulse");
    TClonesArray * tclonesarray_T2_IndPulses = new TClonesArray("RSTPC_Pulse");
    TClonesArray * tclonesarray_T2_Hits 	 = new TClonesArray("RSTPC_Hit");

    //TClonesArray &ColPulses_array = *tclonesarray_T2_ColPulses;
    //TClonesArray &IndPulses_array = *tclonesarray_T2_IndPulses;

   	T2->SetBranchAddress("ColPulses", &tclonesarray_T2_ColPulses);
	T2->SetBranchAddress("IndPulses", &tclonesarray_T2_IndPulses);
    T2->SetBranchAddress("Hits", 	  &tclonesarray_T2_Hits);

    std::cout << " ---------------------------------- " << std::endl;
    std::cout << " Processing T2 with " << T2->GetEntries() << " entries: " << std::endl;
    std::cout << " ---------------------------------- " << std::endl;

	// Check that T1 as well as T2 have the same number of entries
	if( T1->GetEntries() != T2->GetEntries() ) {
		std::cout << " WARNING: T1 AND T2 HAVE DIFFERENT NUMBER OF ENTRIES !! " << std::endl;
	}

	// Loop over all events in T2
    for(int event=4; event<5; event++) { // 2: reference track // 4: muon + delta // 7: muon only
		T2->GetEntry(event);
		/*T2->Show(event);
		std::cout << " ---------------------------------------------------------------------------------------- " << std::endl;*/
		std::cout << "======> EVENT:" << event << std::endl;

		// Get the number of ColPulses, IndPulses and Hits in this event
		UInt_t NColPulses = tclonesarray_T2_ColPulses->GetEntries();
		UInt_t NIndPulses = tclonesarray_T2_IndPulses->GetEntries();
		UInt_t NHits = tclonesarray_T2_Hits->GetEntries();
		std::cout << " NColPulses: " << NColPulses << " \tNIndPulses: " << NIndPulses << " \tNHits: " << NHits << std::endl;

		// Loop over all ColPulses and put the bad pulses into a set
    	std::set<ULong_t> badColIDs;
    	std::set<ULong_t> badIndIDs;
        Int_t width_ColPulse = 0;
		std::cout << " Good ColPulseIDs: " << std::endl << " --------------- " << std::endl;
		for(UInt_t pulse=0; pulse<NColPulses; pulse++) {
			RSTPC_Pulse	* t2event_ColPulse	= (RSTPC_Pulse *)tclonesarray_T2_ColPulses->At(pulse);
			width_ColPulse = t2event_ColPulse->fRedge - t2event_ColPulse->fLedge;
            //std::cout << " fColPulseID: " << t2event_ColPulse->fPulseID << " \tfColCoinNum: " << t2event_ColPulse->fColCoinNum << " \twidth_ColPulse: " << width_ColPulse << " \tfSigma: " << t2event_ColPulse->fSigma << " \tfMax: " << t2event_ColPulse->fMax << std::endl;
            if( t2event_ColPulse->fMax<300 ) {
            	badColIDs.insert( t2event_ColPulse->fPulseID );
        	}
        	else std::cout << " fColPulseID: " << t2event_ColPulse->fPulseID << " \tfColCoinNum: " << t2event_ColPulse->fColCoinNum << "  \tfMax: " << t2event_ColPulse->fMax << "  \twidth_ColPulse: " << width_ColPulse << " \tfSigma: " << t2event_ColPulse->fSigma << "  \tfLedge: " << t2event_ColPulse->fLedge << " \tfRedge: " << t2event_ColPulse->fRedge << std::endl;
		}
		std::cout << std::endl;

		// Loop over all IndPulses and put the bad pulses into a set
        Int_t width_IndPulse = 0;
		Double_t max_Ind_ampl = 0.;
		std::cout << " Good IndPulseIDs: " << std::endl << " --------------- " << std::endl;
		for(UInt_t pulse=0; pulse<NIndPulses; pulse++) {
			RSTPC_Pulse	* t2event_IndPulse	= (RSTPC_Pulse *)tclonesarray_T2_IndPulses->At(pulse);
			width_IndPulse = t2event_IndPulse->fRedge - t2event_IndPulse->fLedge;
            //std::cout << " fIndPulseID: " << t2event_IndPulse->fPulseID << " \tfIndCoinNum: " << t2event_IndPulse->fIndCoinNum << " \twidth_IndPulse: " << width_IndPulse << " \tfSigma: " << t2event_IndPulse->fSigma << " \tfMax: " << t2event_IndPulse->fMax << std::endl;
            if( t2event_IndPulse->fMax<150 ) {
            	badIndIDs.insert( t2event_IndPulse->fPulseID );
        	}
            else std::cout << " fIndPulseID: " << t2event_IndPulse->fPulseID << " \tfIndCoinNum: " << t2event_IndPulse->fIndCoinNum << "  \tfMax: " << t2event_IndPulse->fMax << "  \twidth_IndPulse: " << width_IndPulse << " \tfSigma: " << t2event_IndPulse->fSigma << "  \tfLedge: " << t2event_IndPulse->fLedge << " \tfRedge: " << t2event_IndPulse->fRedge << std::endl;
		}
		std::cout << std::endl;

		// Clear file with 3Dhits, PrincipalComponents and barycentre for offline event display
		ofstream clear_file;
		clear_file.open("3Dhits.txt");
		if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }
		clear_file.open("PrincipalComponents.txt");
		if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }
		clear_file.open("Barycentre.txt");
		if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }
		clear_file.open("Hits_and_Residuals.txt");
		if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }

		// Loop over all hits and look for those with good ColPulseIDs and good IndPulseIDs
		std::vector<float> hits_x, hits_y, hits_z;
        std::vector<float> weights;
		
		for(int hit=0; hit<NHits; hit++) {
			RSTPC_Hit * t2event_Hit = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit);
		    UInt_t ColPulseID = t2event_Hit->fColPulseID;
		    if( badColIDs.find(ColPulseID) == badColIDs.end() ) { // if true: good ColPulse is found
		        UInt_t IndPulseID = t2event_Hit->fIndPulseID;
		        if( badIndIDs.find(IndPulseID) == badIndIDs.end() ) { // if true: good IndPulse is found
					//std::cout << " Good hit. HitID: " << t2event_Hit->fHitID << " \tfCentreTime: " << t2event_Hit->fCentreTime << std::endl;
                    //t2event_Hit->fX // ->fY // ->fZ // ->fMeanTime // ->fCentreTime // ->fColWireNum // ->fIndWireNum // ->fColPulseID // ->fIndPulseID // ->fMeanHeight // ->fLedge // ->fRedge

		            // Only proceed if fX, fY and fMeanTime are > than certain value (or are not -nan)
		            if( t2event_Hit->fX>-999. && t2event_Hit->fY>-999. && t2event_Hit->fCentreTime>-999. ) {
		                // Write spatial coordinates to file for offline event display
		                ofstream textfile;
		                textfile.open("3Dhits.txt", ios::out | ios::app);
		                if(textfile.is_open()) {
		                    textfile << t2event_Hit->fIndWireNum << " " << t2event_Hit->fColWireNum << " " << t2event_Hit->fCentreTime/20 * e_drift_velocity * 10 << " " << t2event_Hit->fMeanHeight << std::endl; // [mm]. Factor 1/20 comes from the fact that 1 sample corresponds to 50 ns, so 1/20 us = 50 ns. // Factor 10: converting cm to mm // WHAT ABOUT USING fMeanTime??
		                    textfile.close();
		                }

		                // Save the 3D hit coordinates (IN UNITS OF WIRE PITCHES FOR X AND Y COORDINATES)
		                hits_x.push_back(t2event_Hit->fIndWireNum);
		                hits_y.push_back(t2event_Hit->fColWireNum);
		                hits_z.push_back(t2event_Hit->fCentreTime/20 * e_drift_velocity * 10); // [mm]. Factor 1/20 comes from the fact that 1 sample corresponds to 50 ns, so 1/20 us = 50 ns. // Factor 10: converting cm to mm // WHAT ABOUT USING fMeanTime??
                        //map<UInt_t, RSTPC_Pulse*> ColPulsesMap;
                        //RSTPC_Pulse* ColPulse = ColPulsesMap[t2event_Hit->fColPulseID];
                        weights.push_back(1.); //ColPulse->fMax);
		            }
		        }
		    }
		}

		std::cout << " Found " << hits_x.size() << " good hits." << std::endl;

		// Access the information from the Hits
		//for(int k=0; k<NHits; k++) {
			//RSTPC_Hit 	* t2event_Hit 		= (RSTPC_Hit *)tclonesarray_T2_Hits	  	  	->At(event+k);
        	//std::cout << " Hit: " << " \tfCentreTime: " << t2event_Hit->fCentreTime << " \tfMeanTime: " << t2event_Hit->fMeanTime << std::endl;
		//}

        for(int entry=0; entry<hits_x.size(); entry++) {
            std::cout << " x: " << hits_x[entry] << " \tweight: " << weights[entry] << std::endl;
        }

		// Erase the sets
		badColIDs.erase( badColIDs.begin(), badColIDs.end() );
		badIndIDs.erase( badIndIDs.begin(), badIndIDs.end() );


		// Perform Principal Components Analysis (PCA)
		// ===========================================
    	std::cout << " ---------------------------------- " << std::endl;
    	std::cout << " Principal Components Analysis: " << std::endl;
    	std::cout << " ---------------------------------- " << std::endl;

		std::vector<double> Lambda_PC_and_barycentre = principal_components(hits_x,hits_y,hits_z,weights);

		std::vector<double> eigenvalues(3,0.);
		std::vector<double> eigenvector1(3,0.);
		std::vector<double> eigenvector2(3,0.);
		std::vector<double> eigenvector3(3,0.);
		std::vector<double> barycentre(3,0.);

		for(int i=0; i<3; i++) {
			eigenvalues[i] 	= Lambda_PC_and_barycentre[i];
			eigenvector1[i] = Lambda_PC_and_barycentre[i+3];
			eigenvector2[i] = Lambda_PC_and_barycentre[i+6];
			eigenvector3[i] = Lambda_PC_and_barycentre[i+9];
			barycentre[i] 	= Lambda_PC_and_barycentre[i+12];
		}

		std::cout << " Eigenvalues:             (" << eigenvalues[0]  << " , " << eigenvalues[1]  << " , " << eigenvalues[2]  << ")" << std::endl;
		std::cout << " 1st eigenvector (x,y,z): (" << eigenvector1[0] << " , " << eigenvector1[1] << " , " << eigenvector1[2] << ")" << std::endl;
		std::cout << " 2nd eigenvector (x,y,z): (" << eigenvector2[0] << " , " << eigenvector2[1] << " , " << eigenvector2[2] << ")" << std::endl;
		std::cout << " 3rd eigenvector (x,y,z): (" << eigenvector3[0] << " , " << eigenvector3[1] << " , " << eigenvector3[2] << ")" << std::endl;
		std::cout << " Barycentre (x,y,z):      (" << barycentre[0]   << " , " << barycentre[1]   << " , " << barycentre[2]   << ")" << std::endl;

		ofstream textfile;
		textfile.open("PrincipalComponents.txt", ios::out | ios::app);
		if(textfile.is_open()) {
			// ONLY DRAW THIRD PC (MUON TRAVELLING IN Z DIRECTION)
		    //textfile << eigenvector1[0] << " " << eigenvector1[1] << " " << eigenvector1[2] << std::endl;
		    //textfile << eigenvector2[0] << " " << eigenvector2[1] << " " << eigenvector2[2] << std::endl;
		    textfile << eigenvector3[0] << " " << eigenvector3[1] << " " << eigenvector3[2] << std::endl;
		    textfile.close();
		}

		textfile.open("Barycentre.txt", ios::out | ios::app);
		if(textfile.is_open()) {
		    textfile << barycentre[0] << " " << barycentre[1] << " " << barycentre[2] << std::endl;
		    textfile.close();
		}


		// Calculate residual vector for each 3D hit to the straight line (the principal component)
		// ========================================================================================
    	std::cout << " ---------------------------------- " << std::endl;
    	std::cout << " Computing the residuals: " << std::endl;
    	std::cout << " ---------------------------------- " << std::endl;

		// Initialize a vector (length = number of hits) of vectors (with length = 3 for the residuals in each spatial direction)
		std::vector<std::vector<double>> residuals(hits_x.size(),std::vector<double>(3,0.));

		// Calculate the residuals for all good hits in the event
		// Note that the residual vector points from the hit to the closest point on the principal component's line
		residuals = calculate_residuals(hits_x,hits_y,hits_z,eigenvector3,barycentre);

		// Write Hits and residuals in .txt file
		textfile.open("Hits_and_Residuals.txt", ios::out | ios::app);
		if(textfile.is_open()) {
			for(int hit=0; hit<hits_x.size(); hit++) {
				//std::cout << " Hit " << hit << " \t \tr_x: " << residuals[hit][0] << " \tr_y: " << residuals[hit][1] << " \tr_z: "  << residuals[hit][2] << std::endl;
				textfile << hits_x[hit] << " " << hits_y[hit] << " " << hits_z[hit] <<  " " << residuals[hit][0] << " " << residuals[hit][1] << " " << residuals[hit][2] << std::endl;
			}
			std::cout << " Hits and residuals have been written to file 'Hits_and_Residuals.txt' " << std::endl;
			textfile.close();
		}


    } // End loop over all events in T2
















	// =============================================================================================
	// =============================================================================================
	// =============================================================================================
	// =============================================================================================


    // Access the T1 and T2 data (via RSTPC_T1wrapper and RSTPC_T2wrapper)
    // ===================================================================
    /*RSTPC_T1wrapper * t1w = new RSTPC_T1wrapper("/home/rberner/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root");
    RSTPC_T2wrapper * t2w = new RSTPC_T2wrapper( t1w ); //Alternatively: RSTPC_T2wrapper * t2w = new RSTPC_T2wrapper( (TTree*)(t1w->fInfile->Get("T2")) );

    // Check if t1w is initialised:
    if( !t1w->IsInit() )  {
      std::cout << " ERROR: t1w is not initialized: " << std::endl;
      return;
    }

    // Read T1 (note: t1w->fChain is a tree)
    std::cout << " ---------------------------------- " << std::endl;
    std::cout << " Total entries in T1: " << t1w->fChain->GetEntries() << std::endl;
    std::cout << " ---------------------------------- " << std::endl;
    std::cout << " Tpc \t Tpc \t\t Rms \t Rms \t\t Feb \t Feb \t\t FebT \t FebT \t FebB \t FebB \t Feb " << std::endl;
    std::cout << " Evnt: \t Time: \t\t ColW: \t IndW: \t\t Evnt: \t Time: \t\t Amp: \t TotAmp: Amp: \t TotAmp: TotAmp: " << std::endl << std::endl;
        for(int entry=0; entry<5; entry++) { // }t1w->fChain->GetEntries(); entry++) {
      t1w->GetEntry(entry);
      //t1w->Show(entry);
      std::cout << " " << t1w->TpcEvent << " \t " << t1w->TpcTime << " \t " << t1w->RmsColWires[0] << " \t "
                << t1w->RmsIndWires[5] << " \t " << t1w->FebEvent << " \t " << t1w->FebTime << " \t "
                << t1w->FebTopAmp[2] << " \t " << t1w->FebTopTotAmp << " \t " << t1w->FebBotAmp[0] << " \t "
                << t1w->FebBotTotAmp << " \t " << t1w->FebTotAmp << std::endl;
    }
    std::cout << std::endl;

    // Read T2
    std::cout << " ---------------------------------- " << std::endl;
    std::cout << " Total entries in T2: " << t2w->fChain->GetEntries() << std::endl;
    std::cout << " ---------------------------------- " << std::endl;
    for(int entry=0; entry<5; entry++) { // }t2w->fChain->GetEntries(); entry++) {
        t2w->GetEntry(entry);
        //t2w->Show(entry);
        std::cout << " Entry: " << entry << " \t GoodEvent: " << t2w->GoodEvent << " \t First ColPulse: " << t2w->ColPulses->First() << " \t n ColPulses: " << t2w->ColPulses->GetEntries() << std::endl;
    }
    std::cout << std::endl;

    // For test purposes: Get t2w entry 2
    std::cout << " ---------------------------------- " << std::endl;
    std::cout << " t2w->entry 2 of " << t2w->fChain->GetEntries() << ": " << std::endl;
    std::cout << " ---------------------------------- " << std::endl;
    t2w->GetEntry(2);
    int NColPulses = t2w->ColPulses->GetEntries();
    int NIndPulses = t2w->IndPulses->GetEntries();
    int NHits      = t2w->Hits->GetEntries();

    // Loop over all hits in t2w entry 2 and print the number of Col- and IndPulses
    //for(int hit=0; hit<5; hit++) { // hit<NHits; hit++) {
    //    std::cout << " hit: " << hit << " \t x: " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fX  << " \t y: " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fY << " \t NColPulses: " << NColPulses << " \t NIndPulses " << NIndPulses << std::endl;
    //}
    //std::cout << std::endl;


    // Put the bad collection and induction pulses into a set
    std::set<ULong_t> badColIDs;
    std::set<ULong_t> badIndIDs;

    for(int counter=0; counter<NColPulses; counter++) {
        Int_t width_Col = ((RSTPC_Pulse *)t2w->ColPulses->At(counter))->fRedge - ((RSTPC_Pulse *)t2w->ColPulses->At(counter))->fLedge;
        Double_t max_Col_ampl = ((RSTPC_Pulse *)t2w->ColPulses->At(counter))->fMax;
        if(max_Col_ampl<300 || width_Col>700) { // HARD CUT, VALUE CHOSEN WITH EVENT 2 OF RUN 2032
            badColIDs.insert( ((RSTPC_Pulse *)t2w->ColPulses->At(counter))->fPulseID  );
        }
        else std::cout << " GOOD COL PULSE ID: " << ((RSTPC_Pulse *)t2w->ColPulses->At(counter))->fPulseID << std::endl;
    }
    for(int counter=0; counter<NIndPulses; counter++) {
        Int_t width_Ind = ((RSTPC_Pulse *)t2w->IndPulses->At(counter))->fRedge - ((RSTPC_Pulse *)t2w->IndPulses->At(counter))->fLedge;
        Double_t max_Ind_ampl = ((RSTPC_Pulse *)t2w->IndPulses->At(counter))->fMax;
        if(max_Ind_ampl<100 || width_Ind>700) { // HARD CUT, VALUE CHOSEN WITH EVENT 1 OF RUN 2032
            badIndIDs.insert( ((RSTPC_Pulse *)t2w->IndPulses->At(counter))->fPulseID   );
        }
        else std::cout << " GOOD IND PULSE ID: " << ((RSTPC_Pulse *)t2w->IndPulses->At(counter))->fPulseID << std::endl;
    }


    // Clear file with 3Dhits for offline event display
    ofstream clear_file;
    clear_file.open("3Dhits.txt");
    if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }

    // Clear file with PrincipalComponents for offline event display
    clear_file.open("PrincipalComponents.txt");
    if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }

    // Clear file with barycentre for offline event display
    clear_file.open("Barycentre.txt");
    if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }


    // Loop over all hits and look for the good hits
    std::vector<double> hits_x, hits_y, hits_z;
    
    for(int hit=0; hit<NHits; hit++) {
        UInt_t ColPulseID = ((RSTPC_Hit *)t2w->Hits->At(hit))->fColPulseID;

        if( badColIDs.find(ColPulseID) == badColIDs.end() ) { // if true: good ColPulse is found
            UInt_t IndPulseID = ((RSTPC_Hit *)t2w->Hits->At(hit))->fIndPulseID;
            if( badIndIDs.find(IndPulseID) == badIndIDs.end() ) { // if true: good IndPulse is found
                //std::cout   << " Good hit: x: "     << ((RSTPC_Hit *)t2w->Hits->At(hit))->fX
                //            << " \t y: "            << ((RSTPC_Hit *)t2w->Hits->At(hit))->fY
                //            << " \t Mean time: "    << ((RSTPC_Hit *)t2w->Hits->At(hit))->fMeanTime
                //            << " \t Centre time: "  << ((RSTPC_Hit *)t2w->Hits->At(hit))->fCentreTime
                //            << " \t ColWireNum: "   << ((RSTPC_Hit *)t2w->Hits->At(hit))->fColWireNum
                //            << " \t IndWireNum: "   << ((RSTPC_Hit *)t2w->Hits->At(hit))->fIndWireNum << std::endl;

                // Only proceed if fX, fY and fMeanTime are > than certain value (or are not -nan)
                if( ((RSTPC_Hit *)t2w->Hits->At(hit))->fX>-999. && ((RSTPC_Hit *)t2w->Hits->At(hit))->fY>-999. && ((RSTPC_Hit *)t2w->Hits->At(hit))->fMeanTime>-999.  ) {
                    // Write spatial coordinates to file for offline event display
                    ofstream textfile;
                    textfile.open("3Dhits.txt", ios::out | ios::app);
                    if(textfile.is_open()) {
                        textfile << ((RSTPC_Hit *)t2w->Hits->At(hit))->fX << " " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fY << " " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fMeanTime/20 * e_drift_velocity << std::endl;
                        textfile.close();
                    }

                    // Save the 3D hit coordinates
                    std::cout << " fColPulseID: " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fColPulseID << " \tfIndPulseID: " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fIndPulseID << " \tfX: " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fX << " \tfY: " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fY << " \tfMeanTime: " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fMeanTime << " \tfCentreTime: " << ((RSTPC_Hit *)t2w->Hits->At(hit))->fCentreTime << std::endl;
                    hits_x.push_back(((RSTPC_Hit *)t2w->Hits->At(hit))->fX);
                    hits_y.push_back(((RSTPC_Hit *)t2w->Hits->At(hit))->fY);
                    hits_z.push_back(((RSTPC_Hit *)t2w->Hits->At(hit))->fMeanTime/20 * e_drift_velocity); // [cm]. Factor 1/20 comes from the fact that 1 sample corresponds to 50 ns, so 1/20 us = 50 ns.
                }   
            }
        }
    }
    */


    /*std::vector<double> a = {1.+5.,2.+5.,3.+5.,4.+5.,5.+5.,6.+5.,7.+5.,8.+5.,9.+5.,10.+5.};
    std::vector<double> b = {2.+5.,4.+5.,6.+5.,8.+5.,10.+5.,12.+5.,14.+5.,16.+5.,18.+5.,20.+5.};
    std::vector<double> c = {3.+5.,6.+5.,9.+5.,12.+5.,15.+5.,18.+5.,21.+5.,24.+5.,27.+5.,30.+5.};

    ofstream textfile;
    textfile.open("3Dhits.txt");
    if(textfile.is_open()) {
        textfile << "";
        textfile.close();
    }

    textfile.open("3Dhits.txt", ios::out | ios::app);
    if(textfile.is_open()) {
        for(int i=0; i<a.size(); i++) {
            textfile << a[i] << " " << b[i] << " " << c[i] << std::endl;
        }
        textfile.close();
    }
    */


    // Perform Principal Components Analysis (PCA)
    // ===========================================
    // HERE: ALSO RETURN THREE EIGENVALUS AND THREE EIGENVECTORS!
    /*
    std::vector<double> PC_and_barycentre = principal_components(hits_x,hits_y,hits_z);
    std::cout << " The principal component vector is calculated to be: " << std::endl;
    std::cout << PC_and_barycentre[0] << std::endl << PC_and_barycentre[1] << std::endl << PC_and_barycentre[2] << std::endl;
    std::cout << " with barycentre at: " << std::endl;
    std::cout << PC_and_barycentre[3] << std::endl << PC_and_barycentre[4] << std::endl << PC_and_barycentre[5] << std::endl;

    ofstream textfile;
    textfile.open("PrincipalComponents.txt", ios::out | ios::app);
    if(textfile.is_open()) {
        textfile << PC_and_barycentre[0] << " " << PC_and_barycentre[1] << " " << PC_and_barycentre[2] << std::endl;
        textfile.close();
    }

    textfile.open("Barycentre.txt", ios::out | ios::app);
    if(textfile.is_open()) {
        textfile << PC_and_barycentre[3] << " " << PC_and_barycentre[4] << " " << PC_and_barycentre[5] << std::endl;
        textfile.close();
    }
    */

}
