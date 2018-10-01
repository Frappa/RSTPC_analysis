#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include <iostream>
#include <algorithm>

#include "Eigen/Dense"
//#include "Eigen/Eigenvalues"
//using namespace Eigen;



// **************************************************************
std::vector<double> principal_components(std::vector<float> hit_x, std::vector<float> hit_y, std::vector<float> hit_z, std::vector<float> weights) {
// **************************************************************

    // First, check that all input vectors have the same size
    std::vector<double> error(3,0.);
    if( hit_x.size()!=hit_y.size() || hit_x.size()!=hit_z.size() || hit_x.size()!=weights.size() || hit_y.size()!=hit_z.size() || hit_y.size()!=weights.size() || hit_z.size()!=weights.size() ) {
        std::cout << " ERROR: size of hit_x, hit_y, hit_z and weights do not match !!" << std::endl;
        return error;
    }


	// Check that no negative weights are present
	for(int entry=0; entry<weights.size(); entry++) {
		if(weights[entry]<0.) std::cout << " WARNING: NEGATIVE WEIGHTS IN PCA MODULE !! " << std::endl;
	}


    // Create vector with the principal component
    std::vector<double> principal_components_vector(9,0.);


    // Get the barycentre of all hits
    std::vector<double> barycentre;
    double sum_x=0, sum_y=0, sum_z=0;
    for(int entry=0; entry<hit_x.size(); entry++) {
        // Reject values which are -nan
        //if( hit_x[entry]>-999. && hit_y[entry]>-999. && hit_z[entry]>-999.   ) {
        sum_x = sum_x + hit_x[entry];
        sum_y = sum_y + hit_y[entry];
        sum_z = sum_z + hit_z[entry];
        //}
    }

    barycentre = { sum_x/hit_x.size(), sum_y/hit_y.size(), sum_z/hit_z.size()  };
    //std::cout << " Barycentre at (x,y,z) = (" << barycentre[0] << ", " << barycentre[1] << ", " << barycentre[2] << ")" << std::endl;


    double sum_of_weights = 0.;

    // Shift the origin of spatial coordinates to the barycentre
    for(int entry=0; entry<hit_x.size(); entry++) {
        // Reject values which are -nan
        //if( hit_x[entry]>-999. && hit_y[entry]>-999. && hit_z[entry]>-999. ) {
            hit_x[entry] = (hit_x[entry] - barycentre[0]) * weights[entry];
            hit_y[entry] = (hit_y[entry] - barycentre[1]) * weights[entry];
            hit_z[entry] = (hit_z[entry] - barycentre[2]) * weights[entry];

            sum_of_weights += weights[entry] * weights[entry];
        //}
    }


    // Create the covariance matrix
    double xixi = 0.;
    double xiyi = 0.;
    double xizi = 0.;
    double yiyi = 0.;
    double yizi = 0.;
    double zizi = 0.;

    for(int entry=0; entry<hit_x.size(); entry++) {
        //if( hit_x[entry]>-999. && hit_y[entry]>-999. && hit_z[entry]>-999.  ) {
            xixi += hit_x[entry] * hit_x[entry];
            xiyi += hit_x[entry] * hit_y[entry];
            xizi += hit_x[entry] * hit_z[entry];
            yiyi += hit_y[entry] * hit_y[entry];
            yizi += hit_y[entry] * hit_z[entry];
            zizi += hit_z[entry] * hit_z[entry];
        //} 
    }

    Eigen::Matrix3f covMatrix;
    covMatrix   <<  xixi, xiyi, xizi,
                    xiyi, yiyi, yizi,
                    xizi, yizi, zizi;

    covMatrix *= 1. / sum_of_weights;

    //std::cout << " The covariance matrix has the form: " << std::endl;
    //std::cout << covMatrix << std::endl;


    // Find eigenvalues (at maximum three) and save them in eigenValColVec
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(covMatrix);

    if(eigenMat.info() == Eigen::ComputationInfo::Success) {
        using eigenValColPair = std::pair<float,size_t>;
        std::vector<eigenValColPair> eigenValColVec = {
            eigenValColPair(eigenMat.eigenvalues()(0),0),
            eigenValColPair(eigenMat.eigenvalues()(1),1),
            eigenValColPair(eigenMat.eigenvalues()(2),2)
        };

        ////std::sort(eigenValColVec.begin(), eigenValColVec.end(),[](const eigenValColPair& left, const eigenValColPair& right){return left.first > right.first;}); // DO NOT SORT BECAUSE EIGENVECTORS AREN'T SORT EITHER
        //std::cout << " with eigenvalues of: " << eigenMat.eigenvalues()(0) << " \t " << eigenMat.eigenvalues()(1) << " \t " << eigenMat.eigenvalues()(2) << std::endl;

        // Find eigenvectors and save them in eigenvectors
        std::vector<std::vector<double>> eigenvectors;
        Eigen::Matrix3f eigenVecs(eigenMat.eigenvectors());
        for(const auto &pair : eigenValColVec) {
            std::vector<double> tempVec = { eigenVecs(0,pair.second), eigenVecs(1,pair.second), eigenVecs(2,pair.second) };
            eigenvectors.push_back(tempVec);
        }

        for(int i=0; i<3; i++) {
            principal_components_vector[i] = eigenvectors[0][i];
            principal_components_vector[i+3] = eigenvectors[1][i];
            principal_components_vector[i+6] = eigenvectors[2][i];
        }

        //std::cout << " 1st eigenvector: " << eigenvectors[0][0] << " " << eigenvectors[0][1] << " " << eigenvectors[0][2] << std::endl;
        //std::cout << " 2nd eigenvector: " << eigenvectors[1][0] << " " << eigenvectors[1][1] << " " << eigenvectors[1][2] << std::endl;
        //std::cout << " 3rd eigenvector: " << eigenvectors[2][0] << " " << eigenvectors[2][1] << " " << eigenvectors[2][2] << std::endl;
    }
    else {
        std::cerr << " WARNING: PCA decompose failure! " << std::endl;
        return error;
    }


    //std::cout << " ======================================================================================= " << std::endl;
    //Eigen::SelfAdjointEigenSolver<Matrix3f> eigensolver(covMatrix);
    //if (eigensolver.info() != Success) abort();
    //std::cout << " The eigenvalues of the covariance matrix are: " << std::endl << eigensolver.eigenvalues() << std::endl;
    //std::cout << " Matrix whose columns are eigenvectors of of the covariance matrix: " << std::endl << std::endl << eigensolver.eigenvectors() << endl;

    

    // Put principal component and barycentre into one vector to return
    std::vector<double> Lambda_PC_and_barycentre(15,0.);
    for(int i=0; i<3; i++) {
        Lambda_PC_and_barycentre[i] = eigenMat.eigenvalues()(i);
        Lambda_PC_and_barycentre[i+3] = principal_components_vector[i];
        Lambda_PC_and_barycentre[i+6] = principal_components_vector[i+3];
        Lambda_PC_and_barycentre[i+9] = principal_components_vector[i+6];
        Lambda_PC_and_barycentre[i+12] = barycentre[i];
    }


    return Lambda_PC_and_barycentre;
} 
