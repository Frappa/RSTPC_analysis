#include <iostream>
#include <fstream>
#include <string>
//using namespace std;


// ***************************************************
void eventDisplay() {
// ***************************************************

    // Get number of lines/hits in the file '3Dhits.txt'
    ULong_t lines = 0;
    std::string line;
    ifstream inputfile;
    inputfile.open("3Dhits.txt");
    if(inputfile.is_open()) {
        while(getline(inputfile,line)) {
          lines++;
        }
        inputfile.close();
    }
    std::cout << " Number of hits: " << lines << std::endl;


    // Create histogram
    TH3F * hist = new TH3F("hist","hist",32,0,32,32,0,32,75,0,150);

    
    // Get the min and max weights for the normalization
    double x,y,z,weight;
    double min_value = 9999.9;
    double max_value = -9999.9;
    inputfile.open("3Dhits.txt");
    if(inputfile.is_open()) {
        while( std::getline(inputfile,line) ) {
            istringstream strStream(line);
            strStream >> x >> y >> z >> weight;
            if ( weight > max_value ) {
                max_value = weight;
            }
            if( weight < min_value ) {
                min_value = weight;
            }
        }
        inputfile.close();
    }
    else { std::cout << " Could not open file '3Dhits.txt'" << std::endl; }

    
    // Fill histogram with 3Dhits
    inputfile.open("3Dhits.txt");
    if(inputfile.is_open()) {
        while( std::getline(inputfile,line) ) {
            istringstream strStream(line);
            strStream >> x >> y >> z >> weight;
            std::cout << " " << x << " " << y << " " << z << " " << weight << std::endl;
            hist->Fill(x,y,z,(weight-min_value+1)); // Shift the smallest value to value 1; z values given in mm
        }
        inputfile.close();
    }
    else std::cout << " Could not open file '3Dhits.txt'" << std::endl;


    // Get barycentre
    std::vector<double> barycentre(3,0.);
    inputfile.open("Barycentre.txt");
    if(inputfile.is_open()) {
        while( std::getline(inputfile,line)  ) {
            istringstream strStream(line);
            strStream >> barycentre[0] >> barycentre[1] >> barycentre[2];
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'Barycentre.txt'" << std::endl;

    
    // Fill principal component line
    TH3F * hist_line = new TH3F("hist_line","hist_line",210,0,32,210,0,32,300,0,150);
    int n = 1000; // number of points to draw for the line
    double PC_x, PC_y, PC_z;
    inputfile.open("PrincipalComponents.txt");
    if(inputfile.is_open()) {
        while( std::getline(inputfile,line) ) {
            istringstream strStream(line);
            strStream >> PC_x >> PC_y >> PC_z;
            //std::cout << " PrincipalComponent: " << PC_x << " \t " << PC_y << " \t " << PC_z << std::endl;
            for(Int_t i=-n; i<n; i++) {
                hist_line->Fill(barycentre[0]+(double)i*0.5*PC_x, barycentre[1]+(double)i*0.5*PC_y, barycentre[2]+(double)i*0.5*PC_z);
            }
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'PrincipalComponents.txt'" << std::endl;


    // Create two new 2D histograms with projection of 3D hist onto xt and yt planes // COLOR PALETTE IS WRONG! PROBABLY HAVE TO TAKE PROJECTION OF TPROFILE?
    //TProfile2D * projXT = hist->Project3DProfile("xz");
    //TProfile2D * projYT = hist->Project3DProfile("yz");


    // Plot 3D histogram
    TCanvas * canvas = new TCanvas("canvas","canvas");
    gStyle->SetOptStat(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetZaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitle("IndWireNum [-]"); //"x [mm]");
    hist->GetYaxis()->SetTitle("ColWireNum [-]"); //"y [mm]");
    hist->GetZaxis()->SetTitle("z [mm] (fMeanTime/20 * drift_vel)");
    hist->SetMarkerColor(kBlue);
    hist_line->SetMarkerColor(kRed);
    hist->Draw("BOX1");
    hist_line->Draw("BOX1 same");
    gPad->RedrawAxis();

    /*
    // Plot projections
    TCanvas * projection1 = new TCanvas("projection1","projection1");
    gStyle->SetOptStat(0);
    projXT->SetTitle("Projection on IndWire-Time-Plane");
    projXT->GetXaxis()->SetTitleOffset(1.2);
    projXT->GetYaxis()->SetTitleOffset(1.2);
    projXT->GetXaxis()->SetTitle("z [mm] (fMeanTime/20 * drift_vel)");
    projXT->GetYaxis()->SetTitle("IndWireNum [-]");
    projXT->Draw("COLZ");

    TCanvas * projection2 = new TCanvas("projection2","projection2");
    gStyle->SetOptStat(0);
    projYT->SetTitle("Projection on ColWire-Time-Plane");
    projYT->GetXaxis()->SetTitleOffset(1.2);
    projYT->GetYaxis()->SetTitleOffset(1.2);
    projYT->GetXaxis()->SetTitle("z [mm] (fMeanTime/20 * drift_vel)");
    projYT->GetYaxis()->SetTitle("ColWireNum [-]");
    projYT->Draw("COLZ");
    */

    return;
}
