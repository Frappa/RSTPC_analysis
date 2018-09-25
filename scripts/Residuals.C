#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include <iostream>
#include <algorithm>



// **************************************************************
std::vector<std::vector<double>> calculate_residuals(std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z, std::vector<double> principal_component, std::vector<double> barycentre) {
// **************************************************************

	// Initialize a vector (length = number of hits) of vectors (with length = 3 for the residuals in each spatial direction)
	std::vector<std::vector<double>> residuals(hits_x.size(),std::vector<double>(3,-1.));


    // First, check that all input vectors have the same size
    //double error[hit_x.size()][3];
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() ) {
        std::cout << " ERROR: size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return residuals;
    }


	// For each hit, calculate the shortest vector from the hit to the principal component's line
	double lambda;

	for(int hit=0; hit<hits_x.size(); hit++) {
		lambda = ( ( 	(hits_x[hit]-barycentre[0])*principal_component[0] +
						(hits_y[hit]-barycentre[1])*principal_component[1] +
						(hits_z[hit]-barycentre[2])*principal_component[2] ) / (sqrt( principal_component[0]*principal_component[0] + principal_component[1]*principal_component[1] + principal_component[2]*principal_component[2] )) );

		residuals[hit][0] = barycentre[0] + lambda * principal_component[0] - hits_x[hit];
		residuals[hit][1] = barycentre[1] + lambda * principal_component[1] - hits_y[hit];
		residuals[hit][2] = barycentre[2] + lambda * principal_component[2] - hits_z[hit];
	}

	return residuals;
}
