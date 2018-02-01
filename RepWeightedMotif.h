/*
 * RepWeightedMotif.h
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#ifndef REPWEIGHTEDMOTIF_H_
#define REPWEIGHTEDMOTIF_H_

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <math.h>       /* fabs */

#include "Eigen/Dense"

extern int minNumberOfSubOpt;

using namespace std;

class RepWeightedMotif {

    Eigen::ArrayXXf matALLMotifWeigthed;
    Eigen::ArrayXXf matALLMotif;

    Eigen::RowVectorXi subOptIdx;
    map<string,float> allStructureFeatureForWeigthOfMotifs;
    vector< map<string,float> > nmotifsForEachStructure;
    vector<string> headers,structures,sequences;
    string headerOfRepresentative;
    string sequenceOfRepresentative;
    string structureOfRepresentative;
    Eigen::RowVectorXf distToCentroid;
    Eigen::RowVectorXf allWeight;

public:

	RepWeightedMotif();
    RepWeightedMotif(map<string,float>, vector<map<string, float> >,vector<string>,vector<string>,vector<string>,int,int);
	virtual ~RepWeightedMotif();

    Eigen::ArrayXXf filterMotifs(Eigen::RowVectorXi,int);
    void computeRepresentative(Eigen::ArrayXXf,vector<string>,vector<string>,vector<string> );

	//Getters
    const map<string, float>& getAllStructureFeatureForWeigthOfMotifs() const {
		return allStructureFeatureForWeigthOfMotifs;
	}

    const Eigen::ArrayXXf& getMatAllMotifWeigthed() const {
		return matALLMotifWeigthed;
	}

    const Eigen::RowVectorXi& getSubOptIdx() const {
		return subOptIdx;
	}

    const string& getHeaderOfRepresentative() const {
		return headerOfRepresentative;
	}

    const string& getSequenceOfRepresentative() const {
		return sequenceOfRepresentative;
	}

    const string& getStructureOfRepresentative() const {
		return structureOfRepresentative;
	}

    const Eigen::RowVectorXf& getDistToCentroid() const {
    return distToCentroid;
	}

    const Eigen::RowVectorXf& getAllWeight() const {
    return allWeight;
	}

    const Eigen::ArrayXXf& getMatALLMotif() const {
    return matALLMotif;
	}


};


#endif /* REPWEIGHTEDMOTIF_H_ */
