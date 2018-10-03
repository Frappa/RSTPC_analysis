#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TObjArray.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <TTree.h>
#include <TFile.h>
#include <vector>

#include "HistoManipulators.hh"
#include "DigitalFilters.hh"
#include "RSTPC_RunProcessor.hh"
#include "RSTPC_Analyser.hh"
#include "MppcTreeWrapper.hh"
#include "RSTPC_T1wrapper.hh"
#include "RSTPC_T2wrapper.hh"
#include "RSTPC_Hits.hh"
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TClassTable.h"

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
                                     
    int run_number = 2032;


	// Produce some plots for initial data inspection
	// (-> which ADC cut values should be chosen, etc.)
	// Functions defined in TestPulses.C
	// ===================================================================
    //PlotColPulses(); // Amplitudes of all ColPulses from the entire run plotted vs. width (= LeftEdge - RightEdge)
    //Plot_ColPulses_Ampl_vs_Width(); // Amplitudes of each single Col Wire plotted vs. width
    //Plot_ColPulses_Ampl_vs_FWHM(); // Amplitudes of each single Col Wire plotted vs. FWHM
    //Plot_ColPulses_Ampl_vs_FWTM(); // Amplitudes of each single Col Wire plotted vs. FWTM
    //Plot_ColPulses_Ampl_vs_sigma(); // Amplitudes of each single Col Wire plotted vs. sigma
    //Plot_IndPulses_Ampl_vs_sigma(); // Amplitudes of each single Ind Wire plotted vs. sigma



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

    T1->SetBranchAddress("TpcEvent",     &TpcEvent);
    T1->SetBranchAddress("TpcTime",      &TpcTime);
    T1->SetBranchAddress("RmsColWires",  &RmsColWires);
    T1->SetBranchAddress("RmsIndWires",  &RmsIndWires);
    T1->SetBranchAddress("FebEvent",     &FebEvent);
    T1->SetBranchAddress("FebTime",      &FebTime);
    T1->SetBranchAddress("FebTopAmp",    &FebTopAmp);
    T1->SetBranchAddress("FebTopTotAmp", &FebTopTotAmp);
    T1->SetBranchAddress("FebBotAmp",    &FebBotAmp);
    T1->SetBranchAddress("FebBotTotAmp", &FebBotTotAmp);
    T1->SetBranchAddress("FebTotAmp",    &FebTotAmp);

    /*
    std::cout << " ---------------------------------------- "            << std::endl;
    std::cout << " Reading T1 with " << T1->GetEntries() << " entries: " << std::endl;
    std::cout << " ---------------------------------------- "            << std::endl;
    //std::cout << " Tpc \t Tpc \t\t Rms \t Rms \t\t Feb \t Feb \t\t FebT \t FebT \t FebB \t FebB \t Feb " << std::endl;
    //std::cout << " Evnt: \t Time: \t\t ColW: \t IndW: \t\t Evnt: \t Time: \t\t Amp: \t TotAmp: Amp: \t TotAmp: TotAmp: " << std::endl << std::endl;
    for(int entry=2; entry<3; entry++) {
        T1->GetEntry(entry);
      	std::cout << " " << TpcEvent       << " \t " << TpcTime      << " \t "     << RmsColWires[0] << " \t "
                         << RmsIndWires[5] << " \t " << FebEvent     << " \t "     << FebTime        << " \t "
                         << FebTopAmp[2]   << " \t " << FebTopTotAmp << " \t "     << FebBotAmp[0]   << " \t "
                         << FebBotTotAmp   << " \t " << FebTotAmp    << std::endl;
        //T1->Show(entry);
    }
    std::cout << std::endl;
    */



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

    std::cout << " ---------------------------------------- "            << std::endl;
    std::cout << " Reading T2 with " << T2->GetEntries() << " entries: " << std::endl;
    std::cout << " ---------------------------------------- "            << std::endl;



    // Check that T1 as well as T2 have the same number of entries
    if( T1->GetEntries() != T2->GetEntries() ) {
        std::cout << " WARNING: T1 AND T2 HAVE DIFFERENT NUMBER OF ENTRIES !! " << std::endl;
    }


    bool single_event = true;

    // Clear file with 3Dhits and residuals
    ofstream clear_file;
    clear_file.open("Hits_and_Residuals.txt");
    if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }


    // Loop over all events in T2
    // -------------------------------------------------------------------
    for(int event=0; event<T2->GetEntries(); event++) { //T2->GetEntries(); event++) { // 2: reference track // 4: muon + delta // 7: muon only // 931: long muon track
        T2->GetEntry(event);
        //T2->Show(event);
        std::cout << "======> EVENT:" << event << std::endl;



        // Get the number of ColPulses, IndPulses and Hits in this event
        UInt_t NColPulses = tclonesarray_T2_ColPulses->GetEntries();
        UInt_t NIndPulses = tclonesarray_T2_IndPulses->GetEntries();
        UInt_t NHits = tclonesarray_T2_Hits->GetEntries();
        std::cout << " NColPulses: " << NColPulses << " \tNIndPulses: " << NIndPulses << " \tNHits: " << NHits << std::endl;



        // Loop over all ColPulses and put those with a small ADC amplitude in a set 'badColPulseIDs' (containing bad ColPulseIDs)
        Double_t ColPulse_ADC_threshold = 100.;
        std::set<ULong_t> badColPulseIDs;

        std::set<ULong_t> ColPulsesWithLowerThreshold;
        //ColPulsesWithLowerThreshold.insert( (ULong_t)0 );
        //ColPulsesWithLowerThreshold.insert( (ULong_t)27 );
        //ColPulsesWithLowerThreshold.insert( (ULong_t)29 );
        //ColPulsesWithLowerThreshold.insert( (ULong_t)30 );

        //std::cout << " ColPulses with amplitudes larger than " << ColPulse_ADC_threshold << " ADC:" << std::endl;
        //std::cout << " ----------------------------------------------- " << std::endl;
        for(UInt_t pulse=0; pulse<NColPulses; pulse++) {
            RSTPC_Pulse * t2event_ColPulse = (RSTPC_Pulse *)tclonesarray_T2_ColPulses->At(pulse);

            ////if ( ColPulsesWithLowerThreshold.find(t2event_ColPulse->fWireNum) == ColPulsesWithLowerThreshold.end()) {
                if( t2event_ColPulse->fMax<ColPulse_ADC_threshold ) {
                    badColPulseIDs.insert( t2event_ColPulse->fPulseID );
                }
            ////}
            ////else {
            ////    if( (t2event_ColPulse->fMax < 40) ) {
            ////        badColPulseIDs.insert( t2event_ColPulse->fPulseID );
            ////    }
            ////}
            
            /*else {
                std::cout << "  fColPulseID: "      << t2event_ColPulse->fPulseID
                          << "  \tfColCoinNum: "    << t2event_ColPulse->fColCoinNum
                          << "  \tfMax: "           << t2event_ColPulse->fMax
                          << "  \twidth_ColPulse: " << t2event_ColPulse->fRedge - t2event_ColPulse->fLedge;
                          << "  \tfSigma: "         << t2event_ColPulse->fSigma
                          << "  \tfLedge: "         << t2event_ColPulse->fLedge
                          << "  \tfRedge: "         << t2event_ColPulse->fRedge      << std::endl;
            }*/
        }
        //std::cout << std::endl;


        // After the threshold cut, only a few hits should be remaining.
        // Loop over all remaining (those are the good ones) ColPulses and put those occurring at the same time (within +/- 1 us) in a set 'CoincidentColPulseIDs'
        // First, have to produce vector-map with the ColPulses (the index corresponds to the ColPulseID)
        vector<RSTPC_Pulse*> ColPulsesVec(NColPulses);
        map<UInt_t, RSTPC_Pulse*> ColPulsesMap;
        for(UInt_t pulse=0; pulse<NColPulses; pulse++) {
            RSTPC_Pulse * t2event_ColPulse = (RSTPC_Pulse *)tclonesarray_T2_ColPulses->At(pulse);
            ColPulsesVec.at(pulse) = t2event_ColPulse;
            ColPulsesMap[t2event_ColPulse->fPulseID] = t2event_ColPulse;
        }

        std::set<ULong_t> CoincidentColPulseIDs;
        for(int hit1=0; hit1<NHits; hit1++) {
            RSTPC_Hit * t2event_Hit1 = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit1);
            RSTPC_Pulse* CollPulse1 = ColPulsesMap[t2event_Hit1->fColPulseID];
            if( ! (badColPulseIDs.find(t2event_Hit1->fColPulseID) == badColPulseIDs.end() ) ) { // if true: Pulse is already in set 'badColPulseIDs'
                continue;
            }
            for(int hit2=0; hit2<NHits; hit2++ ) {
                if(hit1==hit2) continue;
                RSTPC_Hit * t2event_Hit2 = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit2);
                RSTPC_Pulse* CollPulse2 = ColPulsesMap[t2event_Hit2->fColPulseID];
                if( ! (badColPulseIDs.find(t2event_Hit2->fColPulseID) == badColPulseIDs.end() ) ) { // if true: Pulse is already in set 'badColPulseIDs'
                    continue;
                }
                if( abs(CollPulse1->fMaxPos-CollPulse2->fMaxPos)<50 ) { // 50 time samples = 1000 ns = 1 us
                    CoincidentColPulseIDs.insert( t2event_Hit1->fColPulseID );
                    CoincidentColPulseIDs.insert( t2event_Hit2->fColPulseID );
                }
            }
        }


        // Loop over all IndPulses and put those with a small ADC amplitude in a set 'badIndPulseIDs' (containing bad IndPulseIDs)
        Double_t IndPulse_ADC_threshold = 150.;
        std::set<ULong_t> badIndPulseIDs;
        
        //std::cout << " IndPulses with amplitudes larger than " << IndPulse_ADC_threshold << " ADC:" << std::endl;
        //std::cout << " ----------------------------------------------- " << std::endl;
        for(UInt_t pulse=0; pulse<NIndPulses; pulse++) {
            RSTPC_Pulse	* t2event_IndPulse	= (RSTPC_Pulse *)tclonesarray_T2_IndPulses->At(pulse);
            if( t2event_IndPulse->fMax<IndPulse_ADC_threshold ) {
                badIndPulseIDs.insert( t2event_IndPulse->fPulseID );
            }
            /*else {
                std::cout   << "  fIndPulseID: "      << t2event_IndPulse->fPulseID
                            << "  \tfIndCoinNum: "    << t2event_IndPulse->fIndCoinNum
                            << "  \tfMax: "           << t2event_IndPulse->fMax
                            << "  \twidth_IndPulse: " << t2event_IndPulse->fRedge - t2event_IndPulse->fLedge;
                            << "  \tfSigma: "         << t2event_IndPulse->fSigma
                            << "  \tfLedge: "         << t2event_IndPulse->fLedge
                            << "  \tfRedge: "         << t2event_IndPulse->fRedge      << std::endl;
            }*/
        }
        //std::cout << std::endl;


        // Clear files with 3Dhits, Principal components and barycentre
        clear_file.open("3Dhits.txt");
        if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }
        clear_file.open("PrincipalComponents.txt");
        if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }
        clear_file.open("Barycentre.txt");
        if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }


        // Loop over all hits and look for those with good ColPulseIDs and good IndPulseIDs
        std::vector<float> hits_x, hits_y, hits_z, weights;
        for(int hit=0; hit<NHits; hit++) {
            RSTPC_Hit * t2event_Hit = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit);
            UInt_t ColPulseID = t2event_Hit->fColPulseID;
            if( badColPulseIDs.find(ColPulseID) == badColPulseIDs.end() && CoincidentColPulseIDs.find(ColPulseID) == CoincidentColPulseIDs.end()) { // if true: good ColPulse is found
                UInt_t IndPulseID = t2event_Hit->fIndPulseID;
                if( badIndPulseIDs.find(IndPulseID) == badIndPulseIDs.end() ) { // if true: good IndPulse is found
                    //std::cout << " Good hit. HitID: " << t2event_Hit->fHitID
                    //          << " \tfCentreTime: "   << t2event_Hit->fCentreTime << std::endl;
                    //t2event_Hit->fX,->fY,->fZ,
                    //           ->fMeanTime,->fCentreTime,
                    //           ->fColWireNum,->fIndWireNum,->fColPulseID,->fIndPulseID,
                    //           ->fMeanHeight,->fLedge,->fRedge

                    // Only proceed if fX, fY and fMeanTime are > than certain value (or are not -nan)
                    if( t2event_Hit->fX>-999. && t2event_Hit->fY>-999. && t2event_Hit->fCentreTime>-999. ) {
                        // Save the 3D hit coordinates (IN UNITS OF WIRE PITCHES FOR X AND Y COORDINATES)
                        hits_x.push_back(t2event_Hit->fIndWireNum);
                        hits_y.push_back(t2event_Hit->fColWireNum);
                        hits_z.push_back(t2event_Hit->fCentreTime/20 * e_drift_velocity * 10); // [mm]. Factor 1/20 because: 1 sample is 50 ns, so 1/20 us = 50 ns. Factor 10: cm to mm conversion

                        // Give each hit a weight (corresponding to the ADC value, corresponding to ColPulse->fMax)
                        RSTPC_Pulse* CollPulse = ColPulsesMap[t2event_Hit->fColPulseID];
                        weights.push_back(CollPulse->fMax);

                        // Write spatial coordinates to file for offline event display
                        ofstream textfile;
                        textfile.open("3Dhits.txt", ios::out | ios::app);
                        if(textfile.is_open()) {
		                    textfile << t2event_Hit->fIndWireNum                            << " "
                                     << t2event_Hit->fColWireNum                            << " "
                                     << t2event_Hit->fCentreTime/20 * e_drift_velocity * 10 << " "
                                     << CollPulse->fMax                                     << std::endl;
                            textfile.close();
                        }
                    }
                }
            }
        }

        std::cout << " Found " << hits_x.size() << " good hits, rejecting those having > 1 ColPulses within 2.5 us (50samples)." << std::endl;


        // Only consider events with at least 6 hits // ADJUST THIS VALUE!
        if( hits_x.size()<=6 || hits_y.size()<=6 || hits_z.size()<=6 ) continue;


        // Erase the sets
        badColPulseIDs.erase( badColPulseIDs.begin(), badColPulseIDs.end() );
        badIndPulseIDs.erase( badIndPulseIDs.begin(), badIndPulseIDs.end() );


        // Initial Principal Components Analysis (PCA)
        // ===========================================
        std::cout << " ----------------------------------------------- " << std::endl;
        std::cout << " Initial Principal Components Analysis: "          << std::endl;
        std::cout << " ----------------------------------------------- " << std::endl;

        std::vector<double> Lambda_PC_and_barycentre = principal_components(hits_x,hits_y,hits_z,weights);

        std::vector<double> eigenvalues(3,0.);
        std::vector<double> eigenvector1(3,0.);
        std::vector<double> eigenvector2(3,0.);
        std::vector<double> eigenvector3(3,0.);
        std::vector<double> barycentre(3,0.);

        for(int i=0; i<3; i++) {
            eigenvalues[i]  = Lambda_PC_and_barycentre[i];
            eigenvector1[i] = Lambda_PC_and_barycentre[i+3];
            eigenvector2[i] = Lambda_PC_and_barycentre[i+6];
            eigenvector3[i] = Lambda_PC_and_barycentre[i+9];
            barycentre[i]   = Lambda_PC_and_barycentre[i+12];
        }

        std::cout << " Eigenvalues:             (" << eigenvalues[0]  << " , " << eigenvalues[1]  << " , " << eigenvalues[2]  << ")" << std::endl;
        std::cout << " 1st eigenvector (x,y,z): (" << eigenvector1[0] << " , " << eigenvector1[1] << " , " << eigenvector1[2] << ")" << std::endl;
        std::cout << " 2nd eigenvector (x,y,z): (" << eigenvector2[0] << " , " << eigenvector2[1] << " , " << eigenvector2[2] << ")" << std::endl;
        std::cout << " 3rd eigenvector (x,y,z): (" << eigenvector3[0] << " , " << eigenvector3[1] << " , " << eigenvector3[2] << ")" << std::endl;
        std::cout << " Barycentre (x,y,z):      (" << barycentre[0]   << " , " << barycentre[1]   << " , " << barycentre[2]   << ")" << std::endl;

        /*ofstream textfile2;
        textfile2.open("PrincipalComponents.txt", ios::out | ios::app);
        if(textfile2.is_open()) {
            // Only draw the main principal component (which is the third axis: mu travelling in z direction)
            textfile2 << eigenvector3[0] << " " << eigenvector3[1] << " " << eigenvector3[2] << std::endl;
            textfile2.close();
        }

        textfile2.open("Barycentre.txt", ios::out | ios::app);
        if(textfile2.is_open()) {
            textfile2 << barycentre[0] << " " << barycentre[1] << " " << barycentre[2] << std::endl;
            textfile2.close();
        }*/



        // Inserting back those hits which have ColPulse coincidences (stored in the set 'CoincidentColPulseIDs')
        // ======================================================================================================
        // If hit has a distance of < 1 wire pitch (ADJUST VALUE!), reinsert the hit
        std::vector<double> dist_hit_to_PC(3,-999.);
        for(int hit=0; hit<NHits; hit++) {
            RSTPC_Hit * t2event_Hit = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit);
            UInt_t ColPulseID = t2event_Hit->fColPulseID;
            if( CoincidentColPulseIDs.find(ColPulseID) != CoincidentColPulseIDs.end()) { // if true: Hit with ColPulse element set 'CoincidentColPulseIDs' is found
                if( t2event_Hit->fX>-999. && t2event_Hit->fY>-999. && t2event_Hit->fCentreTime>-999. ) {
                    dist_hit_to_PC = calculate_distance(t2event_Hit->fIndWireNum,t2event_Hit->fColWireNum,t2event_Hit->fCentreTime/20*e_drift_velocity*10,eigenvector3,barycentre);
                    //std::cout << " Res_x: " << dist_hit_to_PC[0] << " Res_y: " << dist_hit_to_PC[1] << " Res_z: " << dist_hit_to_PC[2] << std::endl;

                    // Reinsert the hits which have residuals smaller than 1 wire pitch (=52.5/31mm) in all three coordinates (x,y,z)
                    if( fabs(dist_hit_to_PC[0])<1 && fabs(dist_hit_to_PC[1])<1 && fabs(dist_hit_to_PC[2])<52.5/31. ) {
                        hits_x.push_back(t2event_Hit->fIndWireNum);
                        hits_y.push_back(t2event_Hit->fColWireNum);
                        hits_z.push_back(t2event_Hit->fCentreTime/20*e_drift_velocity*10);
                        RSTPC_Pulse* CollPulse = ColPulsesMap[t2event_Hit->fColPulseID];
                        weights.push_back(CollPulse->fMax);
                        std::cout << " Found another good hit: " << hit << std::endl;

                        // Write spatial coordinates to file for offline event display
                        ofstream textfile;
                        textfile.open("3Dhits.txt", ios::out | ios::app);
                        if(textfile.is_open()) {
		                	textfile << t2event_Hit->fIndWireNum                        << " "
                                     << t2event_Hit->fColWireNum                        << " "
                                     << t2event_Hit->fCentreTime/20*e_drift_velocity*10 << " "
                                     << CollPulse->fMax                                 << std::endl;
                            textfile.close();
                        }
                    }
                }
            }
        }


        // Second Principal Components Analysis (PCA)
        // ===========================================
        std::cout << " ----------------------------------------------- " << std::endl;
        std::cout << " 2nd Principal Components Analysis: "              << std::endl;
        std::cout << " ----------------------------------------------- " << std::endl;

        Lambda_PC_and_barycentre = principal_components(hits_x,hits_y,hits_z,weights);

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
            // Only draw the main principal component (which is the third axis: mu travelling in z direction)
            textfile << eigenvector3[0] << " "
                     << eigenvector3[1] << " "
                     << eigenvector3[2] << std::endl;
            textfile.close();
        }

        textfile.open("Barycentre.txt", ios::out | ios::app);
        if(textfile.is_open()) {
            textfile << barycentre[0] << " " << barycentre[1] << " " << barycentre[2] << std::endl;
            textfile.close();
        }



        // Computing residual vector for each 3D hit to the straight line (the principal component)
        // ========================================================================================
        std::cout << " ----------------------------------------------- " << std::endl;
        std::cout << " Computing the residuals: "                        << std::endl;
        std::cout << " ----------------------------------------------- " << std::endl;

        // Initialize a vector (length = number of hits) of vectors (with length = 3 for the residuals in each spatial direction)
        std::vector<std::vector<double>> residuals(hits_x.size(),std::vector<double>(3,0.));

        // Calculate the residuals for all good hits in the event
        // Note that the residual vector points from the hit to the closest point on the principal component's line
        residuals = calculate_residuals(hits_x,hits_y,hits_z,eigenvector3,barycentre); // output: in units of wire pitches (=52.5/31mm)

        // Write Hits and residuals in .txt file
        textfile.open("Hits_and_Residuals.txt", ios::out | ios::app);
        if(textfile.is_open()) {
            for(int hit=0; hit<hits_x.size(); hit++) {
                //std::cout << " Hit "       << hit
                //          << " \t \tr_x: " << residuals[hit][0]
                //          << " \tr_y: "    << residuals[hit][1]
                //          << " \tr_z: "    << residuals[hit][2] << std::endl;
                textfile << hits_x[hit]       << " "
                         << hits_y[hit]       << " "
                         << hits_z[hit]       << " "
                         << residuals[hit][0] << " "
                         << residuals[hit][1] << " "
                         << residuals[hit][2] << std::endl;
            }
            std::cout << " Hits and residuals have been written to file 'Hits_and_Residuals.txt' " << std::endl;
            textfile.close();
        }
        
        
        // Make 2D histograms to show the hits
        // ========================================================================================
        //plot_good_hits_of_event(hits_x,hits_y,hits_z,weights,run_number,event); // Function defined in TestPulses.C

    } // End loop over all events in T2














    // Access the T1 and T2 data (via RSTPC_T1wrapper and RSTPC_T2wrapper)
    // ===================================================================
    /*
    RSTPC_T1wrapper * t1w = new RSTPC_T1wrapper("/home/rberner/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root");
    RSTPC_T2wrapper * t2w = new RSTPC_T2wrapper( t1w );
    //Alternatively: RSTPC_T2wrapper * t2w = new RSTPC_T2wrapper( (TTree*)(t1w->fInfile->Get("T2")) );

    // Check if t1w is initialised:
    if( !t1w->IsInit() )  {
        std::cout << " ERROR: t1w is not initialized: " << std::endl;
        return;
    }
	
    // Read T1 (t1w->fChain is tree)
    std::cout << " Total entries in T1: " << t1w->fChain->GetEntries() << std::endl;
    std::cout << " Tpc \t Tpc \t\t Rms \t Rms \t\t Feb \t Feb \t\t FebT \t FebT \t FebB \t FebB \t Feb " << std::endl;
    std::cout << " Evnt: \t Time: \t\t ColW: \t IndW: \t\t Evnt: \t Time: \t\t Amp: \t TotAmp: Amp: \t TotAmp: TotAmp: " << std::endl << std::endl;
    for(int entry=0; entry<5; entry++) { // }t1w->fChain->GetEntries(); entry++) {
    t1w->GetEntry(entry);
    //t1w->Show(entry);
    std::cout << " " << t1w->TpcEvent       << " \t " << t1w->TpcTime      << " \t " << t1w->RmsColWires[0] << " \t "
                     << t1w->RmsIndWires[5] << " \t " << t1w->FebEvent     << " \t " << t1w->FebTime        << " \t "
                     << t1w->FebTopAmp[2]   << " \t " << t1w->FebTopTotAmp << " \t " << t1w->FebBotAmp[0]   << " \t "
                     << t1w->FebBotTotAmp   << " \t " << t1w->FebTotAmp    << std::endl;
    }

    // Read T2
    std::cout << " Total entries in T2: " << t2w->fChain->GetEntries() << std::endl;
    for(int entry=0; entry<t2w->fChain->GetEntries(); entry++) {
        t2w->GetEntry(entry);
        //t2w->Show(entry);
        std::cout   << " Entry: "             << entry
                    << " \t GoodEvent: "      << t2w->GoodEvent
                    << " \t First ColPulse: " << t2w->ColPulses->First()
                    << " \t n ColPulses: "    << t2w->ColPulses->GetEntries()
                    << std::endl;
    }
    */


}
